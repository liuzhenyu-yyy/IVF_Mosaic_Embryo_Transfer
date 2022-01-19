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

if [ -z $RNA_PIPELINE ];then
	echo "Could not find RNA_PIPELINE !
Please define the RNA_PIPELINE environment variable."
	exit 1
fi

echo "RNA_PIPELINE=${RNA_PIPELINE}" > conf
CONFIG=${RNA_PIPELINE}/${species}_config
echo "CONFIG=${RNA_PIPELINE}/${species}_config" >> conf

while read line
do
	eval "$line"
done < $CONFIG

######Step1 Extract Barcode and UMI######
$umi_tools extract --bc-pattern=CCCCCCCCNNNNNNNN \
                  --stdin ${sample}*2.f*q.gz \
                  --stdout ${sample}.R1.extracted.fq.gz \
                  --read2-stdout \
                  --read2-in ${sample}*1.f*q.gz \
                  --filter-cell-barcode \
                  --whitelist=$barcode

######Step2 Trim TSO and PolyA######
$perl $trim_script ${sample}.R1.extracted.fq.gz ${sample}.R1.trim.fq.gz 0
$seqtk trimfq ${sample}.R1.trim.fq.gz | gzip - > ${sample}.R1.clean.fq.gz

rm ${sample}.R1.extracted.fq.gz \
        ${sample}.R1.trim.fq.gz
