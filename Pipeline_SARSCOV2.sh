#!/bin/bash

#-----------------------------------------------------------------------#
#------------------------ATALHOS PARA O TERMINAL------------------------#
#Para rodar mais de uma amostra
#cd amostras
#criar lista 
#for f in *fastq.gz; do echo ${f%*R*}; done | uniq > ../Lista.txt
#Vai executar o pipe para cada nome presente na lista como $AMOSTRA
#cat Lista.txt | xargs -P 10 -I {} sh -c "bash ./charlie.sh {}"

#-----------------------------------------------------------------------#
#-------------------------NOMEAÇÃO DAS VARIÁVEIS------------------------#
#Cria padrão para nomes de arquivos e pastas de acordo com a amostra analisada

AMOSTRA=$1; #Nome da amostra

R1=${AMOSTRA}R1_001.fastq.gz #Forward paired-end
R2=${AMOSTRA}R2_001.fastq.gz #Reverse paired-end

#---Criando pastas---#
mkdir Trimagem
mkdir Qualidade
mkdir Quali_Trim
mkdir Alinhamento 
mkdir Estatisticas
mkdir Consenso


#-----------------------------------------------------------------------#
#-------------------------------TRIMAGEM--------------------------------#

#Limpar reads  
cd Trimagem
trim_galore -q 33 --illumina --paired ../amostras/${R1} ../amostras/${R2}
cd ../

#-----------------------------------------------------------------------#
#-----------------------VERIFICAÇÃO DA QUALIDADE------------------------#

#Criar arquivo de qualidade das amostras brutas usando FastQC 
cd amostras
fastqc amostras/${R1} -o Qualidade
fastqc amostras/${R2} -o Qualidade
cd ../

#Criar arquivo de qualidade das amostras trimadas usando FastQC
cd Trimagem
fastqc Trimagem/*.fq.gz -o Quali_Trim
cd ../
#-----------------------------------------------------------------------#
#------------------------------ALINHAMENTO------------------------------#

#Criar index de acordo com sequencia referência 
cd Reference_Genomes
bwa index Ref_Wuhan.fasta

cd ../Alinhamento

#Alinha sequências trimadas com referência 
bwa mem -t 1 ../Reference_Genomes/Ref_Wuhan.fasta ../Trimagem/${AMOSTRA}R1_001_val_1.fq.gz ../Trimagem/${AMOSTRA}R2_001_val_2.fq.gz > ${AMOSTRA}.sam 

ls -lh #Check se alinhamento deu certo

#-----------------------------------------------------------------------#
#---------PREPARANDO ARQUIVO PARA CHAMADA E ANÁLISE DE VARIANTES--------#

#Converte arquivo .SAM em . BAM 
samtools view -bS -T ../Reference_Genomes/Ref_Wuhan.fasta ${AMOSTRA}.sam > ${AMOSTRA}.bam

#Classifica o arquivo .BAM do alinhamento 
samtools sort -n ${AMOSTRA}.bam -o ${AMOSTRA}.sorted.bam
#cd Alinhamento

#Cria coordenadas genômicas no arquivo sorted
samtools fixmate ${AMOSTRA}.sorted.bam ${AMOSTRA}.sorted.fixmate.bam

#Reclassifica arquivo 
samtools sort ${AMOSTRA}.sorted.fixmate.bam -o ${AMOSTRA}.sorted.fixmate.position.bam

#indexando o arquivo BAM
samtools index ${AMOSTRA}.sorted.fixmate.position.bam

#-----------------------------------------------------------------------#
#-------------------------CHAMADA DE VARIANTES--------------------------#

#Criar arquivo vcf.gz
samtools mpileup -Ou -f ../Reference_Genomes/Ref_Wuhan.fasta ${AMOSTRA}.sorted.fixmate.position.bam | bcftools call -mv -Oz -o ${AMOSTRA}calls.vcf.gz 

#Indexando arquivo vcf
bcftools index ${AMOSTRA}calls.vcf.gz

#-----------------------------------------------------------------------#
#--------------------------MONTAGEM DO GENOMA---------------------------#

#Faz arquivo consenso entre a amostra e referência, anotando variantes encontradas
cd ../
cd Consenso
cat ../Reference_Genomes/Ref_Wuhan.fasta|bcftools consensus ../Alinhamento/${AMOSTRA}calls.vcf.gz > ${AMOSTRA}.consensus.fa

#Alterar reader do fasta consenso para nome da amostra
sed -i "1s/^.*$/>${AMOSTRA}/" ${AMOSTRA}.consensus.fa

#Concatena os arquivos consenso em um multifasta
cat *.consensus.fa>multifasta-SARS-CoV-2.fa

cd ../

#-----------------------------------------------------------------------#
#-------------------------DETECÇÃO DE VARIANTE--------------------------#

#Para a detecção das linhagens foi uzado o repositório pangolin (https://github.com/cov-lineages/pangolin)
#source usado para conseguir ativar o ambiente pangolin dentro do script
source ~/anaconda3/etc/profile.d/conda.sh

# ativa ambiente conda
conda activate pangolin

#Verifica linhagens do SARS-COV-2 presentes nas amostras
pangolin multifasta-SARS-CoV-2.fa

conda deactivate 

#-----------------------------------------------------------------------#
#---------------------OBTENDO AS ESTATÍSTICAS---------------------------#
#AMOSTRAL: reads das sequências, reads mapeadas e porcentagem de mapeamento
#GERAL: média de reads, cobertura média, mediana da cobertura


cd Estatisticas

#Arquivo geral de estatisticas da amostra
samtools flagstat ../Alinhamento/${AMOSTRA}.sorted.fixmate.position.bam > ${AMOSTRA}mappingstats.txt

#Total de reads da amostra
samtools view -c ../Alinhamento/${AMOSTRA}.sorted.fixmate.position.bam > ${AMOSTRA}.ReadCount

#Reads mapeadas
samtools view -c -F 260 ../Alinhamento/${AMOSTRA}.sorted.fixmate.position.bam > ${AMOSTRA}.ReadsMapped;


x=$(cat ${AMOSTRA}.ReadsMapped); y=$(cat ${AMOSTRA}.ReadCount); python -c "print(round(float(${x}/${y}*100), 2))" > ${AMOSTRA}.PercentMapped;

#--------------------------------------------------------------------------------------------#
#Profundidade
samtools depth -a  ../Alinhamento/${AMOSTRA}.sorted.fixmate.position.bam  |  awk '{sum+=$3} END {print sum/NR}' > ${AMOSTRA}.MeanDepth;

samtools depth -a  ../Alinhamento/${AMOSTRA}.sorted.fixmate.position.bam  |  awk '{print $3}' | sort -n | awk 'NF{a[NR]=$1;c++}END {print (c%2==0)?(a[int(c/2)+1]+a[int(c/2)])/2:a[int(c/2)+1]}' > ${AMOSTRA}.MedianDepth;

samtools depth -a  ../Alinhamento/${AMOSTRA}.sorted.fixmate.position.bam  |  awk '{print $3 >= 10}' | grep '1' | wc -l > ${AMOSTRA}.Depth10;

samtools depth -a  ../Alinhamento/${AMOSTRA}.sorted.fixmate.position.bam  |  awk '{print $3 >= 25}' | grep '1' | wc -l > ${AMOSTRA}.Depth25;

#Coverage
blastn -query ../Consenso/${AMOSTRA}.consensus.fa -db ../Reference_Genomes/Ref_Wuhan -outfmt 6 > ${AMOSTRA}.blastn;

cat ${AMOSTRA}.blastn | awk '{x = $8-$7; print x < 0 ? -x+1 : x+1}' | awk '{sum+=$1} END {coverage = sum/29903 * 100"%"; printf "%0.2f\n", coverage}' > ${AMOSTRA}.CoverageBlastn;

#Contagem de Ns
seqtk comp ../Consenso/${AMOSTRA}.consensus.fa | awk '{x+=$9}END{print x}' > ${AMOSTRA}.CountNs;

#--------------------------------------------------------------------------------------------
#Printar as estatísticas
cd ../Consenso
ls ${AMOSTRA}.consensus.fa > ../Estatisticas/${AMOSTRA}.GenomeName

cd ../Estatisticas
printf "Amostra\tN_Reads\tReads_Mapeados\tPorcentagem_Mapeada\tMean_depth\tMedian_depth\tNpos_Depth>=10\tNpos_Depth>=25\tCoverage\tNumber_of_Ns\n" > ${AMOSTRA}.Statistics;
paste -d "\t" *${AMOSTRA}.GenomeName* *${AMOSTRA}.ReadC* *${AMOSTRA}.ReadsM* *${AMOSTRA}.Percent* *${AMOSTRA}.Mean* *${AMOSTRA}.Median* *${AMOSTRA}.Depth10 *${AMOSTRA}.Depth25 *${AMOSTRA}.CoverageBlastn *${AMOSTRA}.CountN* >> ${AMOSTRA}.Statistics;

#Concatena em arquivo csv estatísticas de todas as amostras

cat ../Estatisticas/*.Statistics >> ALL_GENOMES_STATISTICS.csv;

echo "*************************" 
echo "*${AMOSTRA} terminou*"
echo "*************************"



