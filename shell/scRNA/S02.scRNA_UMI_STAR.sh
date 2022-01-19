#!/bin/sh

############################################
### Tanglab Single-Cell RNA-Seq Pipeline ###
### Author: Zorro Dong                   ###
### Date: 2018-07-25                     ###
############################################

sample=$1
species=$2
Nmax=$3

######Get Parameters######
#source ~/.bashrc
echo ${RNA_PIPELINE}
CONFIG=${RNA_PIPELINE}/${species}_config
while read line
do
        eval "$line"
done < $CONFIG

######Step3 Mapping With STAR######
$STAR --runThreadN 8 \
     --genomeDir $genomeDir \
     --readFilesIn ${sample}.R1.clean.fq.gz \
     --readFilesCommand zcat \
     --outFilterMultimapNmax $Nmax \
     --outFileNamePrefix ${sample}. \
     --outSAMtype BAM SortedByCoordinate

######Step4 Add feature using featureCounts of subread#######
$featureCounts -a $gtf \
	-o gene_assigned \
	-R BAM ${sample}.Aligned.sortedByCoord.out.bam \
	-T 8

#######Step5 sort and index bam file######
$samtools sort \
	-m 15000000000 \
	${sample}.Aligned.sortedByCoord.out.bam.featureCounts.bam \
	-o ${sample}.assigned_sorted.bam

$samtools index \
	${sample}.assigned_sorted.bam

rm ${sample}.Aligned.sortedByCoord.out.bam.featureCounts.bam \
	${sample}.Aligned.sortedByCoord.out.bam
