for i in `cat sample.xls`
do
    echo 'sh /home/software/pipeline/Tanglab_Single_Cell_RNA_Seq/S01.scRNA_UMI_QC.sh' $i > $i/S01.work.sh && \
    echo 'sh /home/software/pipeline/Tanglab_Single_Cell_RNA_Seq/S02.scRNA_UMI_STAR.sh' $i > $i/S02.work.sh && \
    echo 'sh /home/software/pipeline/Tanglab_Single_Cell_RNA_Seq/S03.scRNA_UMI_UMI.sh' $i > $i/S03.work.sh
    echo 'qsub -cwd -l vf=1g,p=1 -V S01.work.sh o>S01.jid' > $i/scRNA_work.sh
    echo "id1=\`perl -n -e '/(\d+)/ && print \$1 ' S01.jid\`" >> $i/scRNA_work.sh
    echo 'qsub -cwd -l vf=28g,p=4 -V -hold_jid $id1 S02.work.sh o>S02.jid' >> $i/scRNA_work.sh
    echo "id2=\`perl -n -e '/(\d+)/ && print \$1 ' S02.jid\`" >> $i/scRNA_work.sh
    echo 'qsub -cwd -l vf=1g,p=1 -V -hold_jid $id2 S03.work.sh' >> $i/scRNA_work.sh
    cd $i && sh scRNA_work.sh && cd ..
done


