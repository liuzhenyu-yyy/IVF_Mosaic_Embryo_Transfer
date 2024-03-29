args=commandArgs(T)

RNA_PIPELINE <- Sys.getenv("RNA_PIPELINE")

print(RNA_PIPELINE)

if ( args[2] == "human" ) {
	gene_list <- "/GRCh38_gene_list.xls"
} else if ( args[2] == "mouse" ) {
	gene_list <- "/mm10_gene_list.xls"
}

print(gene_list)

gene <- read.table(paste(RNA_PIPELINE,gene_list,sep=""),header=T,row.names=1)
barcode <- read.table(paste(RNA_PIPELINE,"/barcode_96_8bp.txt",sep=""),row.names = 1)

count <- as.matrix(read.table(paste(args[1],"UMI_counts.tsv",sep="."),header=T,row.names = 1))

sample_name <- paste(args[1],as.character(barcode[as.character(colnames(count)),]),sep = ".")

colnames(count) <- sample_name

matrix_all_gene <- matrix(0,nrow = nrow(gene),ncol =length(sample_name))
colnames(matrix_all_gene) <- sample_name
rownames(matrix_all_gene) <- as.character(rownames(gene))

matrix_not <- matrix_all_gene[!rownames(matrix_all_gene) %in% rownames(count),]

UMI_count <- rbind(count,matrix_not)

UMI_count_final <- UMI_count[as.character(rownames(gene)),]

rownames(UMI_count_final) <- as.character(gene$Gene_id_unique)

write.table(cbind(Gene=as.character(rownames(UMI_count_final)),UMI_count_final),
            paste(args[1],"umi_counts.xls",sep="."),sep="\t",quote =F,row.names=F)

tran_num <- colSums(UMI_count_final)
gene_num <- colSums(UMI_count_final>0)

gene_tran_num <- cbind(Gene_Number=gene_num,
                       Transcript_Number=tran_num)
write.table(cbind(Sample=rownames(gene_tran_num),gene_tran_num),paste(args[1],"gene_tran_num.xls",sep="."),sep="\t",quote =F,row.names=F)


