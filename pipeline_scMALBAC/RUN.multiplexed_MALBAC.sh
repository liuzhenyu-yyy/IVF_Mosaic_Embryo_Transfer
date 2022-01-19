#!/bin/sh

##########################################################
### Pipeline for multiplexed MALBAC
###
### Author: Shuhui Bian
### Date:   Mon Apr  1 13:23:25 CST 2019
### Email:  fondofbiology@163.com
###
### Modified by Zhenyu Liu 2020-9-22 21:41 CST
##########################################################


sample=$1
ref=$2

###### Reference and scripts######
ref_fa=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/database/$ref/BWA/*.fa
bin_dir=/gpfs1/tangfuchou_pkuhpc/Pipeline/01.scmMALBAC_bsh
sft_dir=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software
barcode=$bin_dir/ZY_barcode_withGAT_whitelist

main_dir=`pwd`
bwa_dir=$main_dir/01.BWA
if [ ! -d $bwa_dir ];then mkdir -p $bwa_dir;fi
cnv_dir=$main_dir/02.CNV
if [ ! -d $cnv_dir ];then mkdir -p $cnv_dir;fi


###### Softwares by order
umi_tools=$sft_dir/miniconda3/bin/umi_tools
seqtk=$sft_dir/seqtk/seqtk 
fastp=$sft_dir/miniconda3/bin/fastp 
bwa=$sft_dir/bwa/bwa
samtools=$sft_dir/samtools/bin/samtools
#java
JAVA_HOME=/apps/usr/java/jdk1.8.0_112
export JAVA_HOME
export PATH=$JAVA_HOME/bin:$JAVA_HOME/jre/bin:$PATH
export JRE_HOME=/apps/usr/java/jdk1.8.0_112/jre
JAVA_OPTS=-server
export CLASSPATH=.:$JAVA_HOME/lib/dt.jar:$JAVA_HOME/lib/tools.jar:$JAVA_HOME/jre/lib/rt.jar:$JAVA_HOME/lib/comm.jar:/usr/share/ant/lib/ant-launcher.jar
java=/apps/usr/java/jdk1.8.0_112/bin/java
picard=$sft_dir/picard_2.19.0/picard.jar
sinto=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/miniconda3/bin/sinto

# Read1:
# Read2: barcode + CAT + N5 + 3C


#==========================================================================
#      Step1 Extract Barcode and UMI
#==========================================================================
$umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCNNNNNNNN           \
                   --stdin ${sample}_*2.f*q.gz                                      \
                   --stdout $sample.R2.extracted.fq.gz                                \
                   --read2-in ${sample}_*1.f*q.gz                                   \
                   --read2-out $sample.R1.extracted.fq.gz                             \
                   --filter-cell-barcode                                              \
                   --whitelist=$barcode


##==========================================================================
##      Step1: trim reads 2 43bp (barcode + CAT + N5 + 3C)
##==========================================================================
#$seqtk trimfq -b 43 ${sample}_R2.fastq.gz > $sample.R2.trim.fq && gzip -f $sample.R2.trim.fq


#==========================================================================
#      Step2: QC for DNA (paired end mode)
#==========================================================================
$fastp -i $sample.R1.extracted.fq.gz                                       \
       -I $sample.R2.extracted.fq.gz                                       \
       -o $sample.R1.clean.fq.gz                                           \
       -O $sample.R2.clean.fq.gz
rm -f $sample.R1.extracted.fq.gz $sample.R2.extracted.fq.gz

#==========================================================================
#      Step3: BWA align and rmdup
#==========================================================================

#+++++++++++++++++++++
# BWA Mapping
#+++++++++++++++++++++
echo -e "["$(date)"]\tStart $sample aligning.."

# attention: t means cpu number
$bwa mem -M -t 15                                                           \
         $ref_fa                                                           \
         $sample.R1.clean.fq.gz                                            \
         $sample.R2.clean.fq.gz                                          | \
         $samtools view -bS -q 30 - > $bwa_dir/${sample}_mapQ30.bam
echo -e "["$(date)"]\t$sample BWA done!"


#+++++++++++++++++++++
# sort and index bam
#+++++++++++++++++++++
echo -e "["$(date)"]\tBegin samtools sort of $sample.."
$samtools sort $bwa_dir/${sample}_mapQ30.bam -o $bwa_dir/${sample}_mapQ30_sort.bam -m 20G
$samtools index $bwa_dir/${sample}_mapQ30_sort.bam
echo -e "["$(date)"]\tFinish samtools sorting and indexing of $sample.."

#+++++++++++++++++++++
# Picard remove Duplicates(after sorting)
#+++++++++++++++++++++
echo -e "["$(date)"]\tMarking duplicates of $sample.."
$java -Xmx20g                                                              \
      -jar $picard MarkDuplicates                                          \
      I=$bwa_dir/${sample}_mapQ30_sort.bam                                 \
      O=$bwa_dir/${sample}_mapQ30_sort_rmdup.bam                           \
      M=$bwa_dir/${sample}_dup.metrics                                     \
      CREATE_INDEX=true                                                    \
      ASSUME_SORTED=true                                                   \
      VALIDATION_STRINGENCY=SILENT                                         \
      REMOVE_DUPLICATES=true
echo -e "["$(date)"]\tFinish marking duplicates of $sample."
rm $bwa_dir/${sample}_mapQ30_sort.bam $bwa_dir/${sample}_mapQ30_sort.bam.bai

#==========================================================================
#      Step5: Add barcode tag (CB), split bam
#==========================================================================
# add cell barcode tag
#echo -e "["$(date)"]\tAdding cell barcode tag to bam.."
#$samtools view $bwa_dir/${sample}_mapQ30_sort_rmdup.bam -h | \
#    awk '{if ($1~/^@/) {print $0} else {split($1,a,"_"); bar=a[2]; print $0"\tCB:Z:"bar}}' | \
#    $samtools view -b > $bwa_dir/${sample}_mapQ30_sort_rmdup_tag.bam

#rm $bwa_dir/${sample}_mapQ30_sort_rmdup.bam $bwa_dir/${sample}_mapQ30_sort_rmdup.bai
#mv $bwa_dir/${sample}_mapQ30_sort_rmdup_tag.bam $bwa_dir/${sample}_mapQ30_sort_rmdup.bam

#$samtools index $bwa_dir/${sample}_mapQ30_sort_rmdup.bam

#echo -e "["$(date)"]\tFinish adding cell barcode tag."

# split bam by barcode
#echo -e "["$(date)"]\tSpliting bam by CB tag..."
#cd $bam_dir
#$sinto filterbarcodes \
#    -b $bwa_dir/${sample}_mapQ30_sort_rmdup.bam \
#    -c $bin_dir/ZY_barcode_withGAT_whitelist \
#    -p 20 \
#    --barcodetag CB

#rename sc ${sample}_ sc*
#cd ..

#==========================================================================
#      Step4: (Optional) bed from bam
#==========================================================================

/lustre1/tangfuchou_pkuhpc/Software/bedtools2/bin/bamToBed -i $bwa_dir/${sample}_mapQ30_sort_rmdup.bam |gzip > $cnv_dir/${sample}_mapQ30_sort_rmdup.bed.gz

# split bed
export R_HOME=/apps/bioinfo/R-3.2.0
export qATH=$R_HOME/bin:$PATH

cd $cnv_dir
$R_HOME/bin/Rscript $bin_dir/split_barcode_bed_CLS.R $sample $bin_dir

cat merge_sampleInfo |while read line;do seq=`echo $line |awk '{print $2}' ` && sam=`echo $line |awk '{print $3}' ` && echo "zcat ${sample}_mapQ30_sort_rmdup.bed.gz |grep $seq | gzip > $sam.bed.gz";done > $sample.grep_work.sh

perl $bin_dir/multi-process.pl -cpu 20 $sample.grep_work.sh
