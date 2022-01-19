#!/bin/sh

##########################################################
### Pipeline for CNV-Seq data from TE PGT-A
###
### Author: Zhenyu Liu
### Date:   27 September 2021
### Email:  liuzhenyu@pku.edu.cn
###
##########################################################

cell=$1


###### Softwares by order
bwa=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/bwa/bwa
samtools=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/samtools/bin/samtools

database_dir=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/database
ref_fa=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/database/hg38/BWA/Homo_sapiens.GRCh38.dna.primary_assembly.fa

#java
JAVA_HOME=/apps/usr/java/jdk1.8.0_112
export JAVA_HOME
export PATH=$JAVA_HOME/bin:$JAVA_HOME/jre/bin:$PATH
export JRE_HOME=/apps/usr/java/jdk1.8.0_112/jre
JAVA_OPTS=-server
export CLASSPATH=.:$JAVA_HOME/lib/dt.jar:$JAVA_HOME/lib/tools.jar:$JAVA_HOME/jre/lib/rt.jar:$JAVA_HOME/lib/comm.jar:/usr/share/ant/lib/ant-launcher.jar
java=/apps/usr/java/jdk1.8.0_112/bin/java
picard=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/picard_2.19.0/picard.jar

mkdir -p ${cell}
cd ${cell}

###### bwa mapping
$bwa aln ${ref_fa} ../data/${cell}.fastq \
    > ${cell}.sai

$bwa samse ${ref_fa} ${cell}.sai ../data/${cell}.fastq | \
     $samtools view -bS -q 30 - > ${cell}_mapQ30.bam

###### samtools sort and index
$samtools sort ${cell}_mapQ30.bam -o ${cell}_mapQ30_sort.bam -m 8G
$samtools index ${cell}_mapQ30_sort.bam

###### picard remove duplication
$java -Xmx20g                                                              \
      -jar $picard MarkDuplicates                                          \
      I=${cell}_mapQ30_sort.bam                                 \
      O=${cell}_mapQ30_sort_rmdup.bam                           \
      M=${cell}_dup.metrics                                     \
      CREATE_INDEX=true                                                    \
      ASSUME_SORTED=true                                                   \
      VALIDATION_STRINGENCY=SILENT                                         \
      REMOVE_DUPLICATES=true

rm ${cell}_mapQ30_sort.bam ${cell}_mapQ30_sort.bam.bai

###### bam to bed
/lustre1/tangfuchou_pkuhpc/Software/bedtools2/bin/bamToBed -i ${cell}_mapQ30_sort_rmdup.bam |gzip > ${cell}_mapQ30_sort_rmdup.bed.gz
