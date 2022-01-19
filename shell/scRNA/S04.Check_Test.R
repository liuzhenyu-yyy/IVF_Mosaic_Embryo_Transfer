getwd()

args=commandArgs(T)

sample=args[1]

umi_c.n <- paste("ref_result/",sample,"/",sample,".umi_counts.xls",sep="")
umi_e.n <- paste(sample,"/",sample,".umi_counts.xls",sep="")

sum_c.n <- paste("ref_result/",sample,"/",sample,".gene_tran_num.xls",sep="")
sum_e.n <- paste(sample,"/",sample,".gene_tran_num.xls",sep = "")


umi_c <- read.csv(umi_c.n,row.names = 1,header = T,sep = "\t")
umi_e <- read.csv(umi_e.n,row.names = 1,header = T,sep = "\t")
sum_c <- read.csv(sum_c.n,sep = "\t")
sum_e <- read.csv(sum_e.n,sep = "\t")

#umi_c[1:5,1:5]
#umi_e[1:5,1:5]

#head(sum_c)
#head(sum_e)

#print("Which Umi_C != Umi_E")
#which(umi_c != umi_e)
#print("Which Sum_C != Sum_E")
#which(sum_c[,2:3] != sum_e[,2:3])

#print("Num Umi_C != Umi_E")
#sum(umi_c != umi_e)
#print("Num Sum_C != Sum_E")
#sum(sum_c[,2:3] != sum_e[,2:3])

result <- sum(umi_c != umi_e) + sum(sum_c[,2:3] != sum_e[,2:3])

if ( result == 0 ){
	out <- paste(sample,"OK",sep="\t")
}else {        
	out <- paste(sample,"ERROR",sep="\t")
}              

write.table(out,"check_result.txt",append=TRUE,quote=F,sep="\t",row.names=F,col.names=F)
