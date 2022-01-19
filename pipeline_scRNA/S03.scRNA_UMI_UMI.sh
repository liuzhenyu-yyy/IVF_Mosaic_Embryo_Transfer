#!/bin/sh

############################################
### Tanglab Single-Cell RNA-Seq Pipeline ###
### Author: Zorro Dong                   ###
### Date: 2018-07-25                     ###
############################################

sample=$1
species=$2

######Get Parameters######
#source ~/.bashrc
echo ${RNA_PIPELINE}
CONFIG=${RNA_PIPELINE}/${species}_config
while read line
do
        eval "$line"
done < $CONFIG

######Step6 get umi count table######
$umi_tools count \
	--per-gene --gene-tag=XT \
	--per-cell --wide-format-cell-counts \
	-I ${sample}.assigned_sorted.bam \
	-S ${sample}.UMI_counts.tsv

#######Step7 Normalization umi count table######
$Rscript $norm_script $sample $species
