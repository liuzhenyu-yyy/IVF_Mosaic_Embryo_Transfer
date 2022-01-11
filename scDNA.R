setwd("E:/LabWork/Project/fetal_genome/Merge/01.scDNA/")
source("E:/LabWork/code/MyFunction.R")

library(ggplot2)
library(ggExtra)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(GenomicRanges)
library(ComplexHeatmap)
library(circlize)
library(ggsci)

# 1. load data ---------------------------------------------------------------

sample.ls <- dir("data/SegCopy/")
sample.ls <- gsub("_Seg.+?$","",sample.ls) %>% unique()

bin_name <- paste("bin_",1:2808,sep = '')

bin.info <- read.table("data/SegCopy/DNA_A1_lym_1_FKDL210056832-1a_SegCopy",
                       header = T,row.names = bin_name)[,1:3]

# read SegNorm (reads at each bin)
SegNorm <- bin.info
for (one in sample.ls) {
   temp <- read.table(paste("data/SegNorm/",one,"_SegNorm",sep = ''),
                      header = T,row.names = bin_name)
   SegNorm <- cbind(SegNorm, temp[,4:ncol(temp)])
}
SegNorm <- SegNorm[,4:ncol(SegNorm)]

# read SegCopy (reads at each bin)
SegCopy <- bin.info
for (one in sample.ls) {
   temp <- read.table(paste("data/SegCopy/",one,"_SegCopy",sep = ''),
                      header = T,row.names = bin_name)
   SegCopy <- cbind(SegCopy, temp[,4:ncol(temp)])
}
SegCopy <- SegCopy[,4:ncol(SegCopy)]

# read SegRef (reads at each bin)
SegRef <- bin.info
for (one in sample.ls) {
   temp <- read.table(paste("data/SegNorm_ref/",one,"_SegNorm_ref",sep = ''),
                      header = T,row.names = bin_name)
   SegRef <- cbind(SegRef, temp[,4:ncol(temp)])
}
SegRef <- SegRef[,4:ncol(SegRef)]

# read SegRef (reads at each bin)
sample.info <- data.frame()
for (one in sample.ls) {
   temp <- read.table(paste("data/SegStat/",one,"_SegStats.tsv",sep = ''),
                      header = T)
   sample.info <- rbind(sample.info, temp)
}

sample.info$Cell <- rownames(sample.info)

sample.info$Sample <- substr(sample.info$Cell,5,6)
table(sample.info$Sample)

sample.info$Sample_2 <-  sample.info$Sample
sample.info$Sample_2 <- gsub("A","P0",sample.info$Sample_2)
sample.info$Sample_2 <- gsub("B","N0",sample.info$Sample_2)
sample.info$Sample_2 <- gsub("P08","P07",sample.info$Sample_2)
table(sample.info$Sample_2)

sample.info$Sample <- factor(sample.info$Sample,
                             levels = c(paste0("A",c(1:6,8)),
                                        paste0("B",1:3)))
sample.info$Sample_2 <- factor(sample.info$Sample_2,
                             levels = c(paste0("P0",c(1:7)),
                                        paste0("N0",1:3)))

sample.rename <- levels(sample.info$Sample_2)
names(sample.rename) <- levels(sample.info$Sample)

saveRDS(sample.info,"data/sample.info.rds")
saveRDS(SegCopy,"data/SegCopy.rds")
saveRDS(SegNorm,"data/SegNorm.rds")
saveRDS(SegRef,"data/SegRef.rds")

bin.chr19 <- rownames(bin.info)[bin.info$CHR=="chr19"]

SegCopy.alt <- SegCopy
SegCopy.alt[bin.chr19,] <- 2

patient.color <- brewer.pal(10,"Set3")
names(patient.color) <- levels(sample.info$Sample_2)

# 2. QC -------------------------------------------------------------------
cutoff.reads <- 300000
cutoff.CV <- 1

pdf("QC.Stat.pdf",6,5)
p1 <- ggplot(sample.info)+
   geom_violin(aes(x=Sample_2,y=log10(Reads),fill=Sample_2),alpha=0.7,show.legend = F)+
   geom_boxplot(aes(x=Sample_2,y=log10(Reads),color=Sample_2),alpha=0,show.legend = F)+
   theme_bw()+theme(panel.grid = element_blank())+
   geom_hline(yintercept = log10(cutoff.reads), lty=2, color = "red")+
   scale_fill_manual(values = patient.color)+
   scale_color_manual(values = patient.color)
p2 <- ggplot(sample.info)+
   geom_violin(aes(x=Sample_2,y=log10(Disp+1),fill=Sample_2),alpha=0.7,show.legend = F)+
   geom_boxplot(aes(x=Sample_2,y=log10(Disp+1),color=Sample_2),alpha=0,show.legend = F)+
   theme_bw()+theme(panel.grid = element_blank())+
   geom_hline(yintercept = log10(cutoff.CV + 1), lty=2, color = "red")+
   scale_fill_manual(values = patient.color)+
   scale_color_manual(values = patient.color)
multiplot(p1, p2, cols = 1)
dev.off()

pdf("QC_scatter.pdf",4,4)
p <- ggplot(sample.info)+
   geom_point(aes(x=log10(Reads),y=log10(Disp)),
              color = ifelse(sample.info$Reads>cutoff.reads&sample.info$Disp<cutoff.CV,"red","gray"))+
   geom_vline(xintercept = log10(cutoff.reads), color="blue",lty=2)+
   geom_hline(yintercept = log10(cutoff.CV), color="blue",lty=2)+
   theme_bw()+theme(panel.grid = element_blank())
ggMarginal(p, type="density", fill="#f8766d", alpha=0.5)
dev.off()

table(ifelse(sample.info$Reads>cutoff.reads&sample.info$Disp<cutoff.CV,"red","gray")) # 479   1729  

cell.qc = sample.info[sample.info$Reads>cutoff.reads & 
                        sample.info$Disp< cutoff.CV,]$Cell

SegCopy <- SegCopy[,cell.qc]
SegCopy.alt <- SegCopy.alt[,cell.qc]
SegNorm <- SegNorm[,cell.qc]
SegRef <- SegRef[,cell.qc]

bin.auto <- rownames(bin.info[!bin.info$CHR%in%c("chrX","chrY"),])
hist(log10(colSums(SegCopy[bin.auto,]!=2)))

cell.qc = colnames(SegCopy)[colSums(SegCopy[bin.auto,]!=2) < 100]
QC.tab <- table(sample.info[cell.qc,]$Sample) 
QC.pct <- QC.tab / table(sample.info$Sample)

round(QC.pct,2)

sample.info <-sample.info[cell.qc,]
SegCopy <- SegCopy[,cell.qc]
SegCopy.alt <- SegCopy.alt[,cell.qc]
SegNorm <- SegNorm[,cell.qc]
SegRef<- SegRef[,cell.qc]

SegCopy.filter <- myFilterCNV(SegCopy.alt, bin.info, 10)
rm(SegCopy.alt)

# 3. CNV profile ----------------------------------------------------------
# plot CN profile for each cell
chr.color <- rep(c("gray85","gray75"),12)
names(chr.color) <- unique(bin.info$CHR)

Samples <- unique(sample.info$Sample_2)
sex <- c("XX","XY","XX","XY","XY","XY","XY","XX","XY","XX")
names(sex) <- Samples

for (one in Samples) {
   cells <- sample.info[sample.info$Sample_2 == one,]$Cell
   
   pdf(paste("CN_",one,"_all.pdf",sep = ""),10,2)
   for (cell in cells) {
      p <- plotCNprofile(segref = SegRef[,cell],segcopy=SegCopy.filter[,cell],
                    cell = cell, bininfo = bin.info, sex = sex[one],return = T)
      plot(p)
   }
   dev.off()
}

#CN cells

cell.cnv.all <- names(which(colSums(SegCopy.filter[bin.auto,]!=2)>0))

pdf("CN_CNV_all.Cells.pdf",10,2) 
for (cell in cell.cnv.all) {
   one <- substr(cell,5,6)
   p <- plotCNprofile(segref = SegRef[,cell],segcopy=SegCopy.filter[,cell],
                      cell = cell, bininfo = bin.info, sex = sex[sample.rename[one]],return = T)
   plot(p)
}
dev.off()

# 4. manully select --------------------------------------------------------
library(Vennerable)

cnv.manull <- read.csv("cell.cnv.merge.csv")
cnv.manull$nPerson <- rowSums(cnv.manull[,2:4])
   
cor(cnv.manull$LZY,cnv.manull$GY)
cor(cnv.manull$LZY,cnv.manull$QSY)
cor(cnv.manull$GY,cnv.manull$QSY)

pdf("Venn.CNV.manull.pdf",5,4)
plot(Venn(list("GY"=cnv.manull$Cell[cnv.manull$GY==1],
               "QSY"=cnv.manull$Cell[cnv.manull$QSY==1],
               "LZY"=cnv.manull$Cell[cnv.manull$LZY==1])))
dev.off()

cell.cnv.3A <- cnv.manull$Cell[cnv.manull$nPerson >= 3]

pdf("CN_CNV_3A_Cells.pdf",10,2) 
for (cell in cell.cnv.3A) {
   one <- substr(cell,5,6)
   p <- plotCNprofile(segref = SegRef[,cell],segcopy=SegCopy.filter[,cell],
                      cell = cell, bininfo = bin.info, sex = sex[sample.rename[one]],return = T)
   plot(p)
}
dev.off()

write.csv(cell.cnv.3A,"cell.cnv.3A.csv",row.names = F,col.names = F,quote = F)

pdf("Fraction.CNV.3A.pdf",3,2.5)
view.data <- as.data.frame(table(sample.info[cell.cnv.3A,]$Sample_2))
colnames(view.data) <- c("Patient", "nCell_CNV")
table(sample.info[cell.cnv.3A,]$Sample)

view.data$nCell_QC <- as.numeric(table(sample.info$Sample))
view.data$CNV_Precent <- round(view.data$nCell_CNV*100 / view.data$nCell_QC, 2)

ggplot(view.data)+
   geom_bar(aes(x=Patient, y=CNV_Precent, fill=Patient),
            stat = "identity")+
   geom_text(aes(x=Patient, y=CNV_Precent+0.1,
                 label=paste0(round(CNV_Precent,1),"%")),
             size = 2)+
   scale_fill_manual(values = patient.color)+
   theme_bw()+theme(panel.grid = element_blank())+
   theme(axis.text.x  = element_text(angle = 50,hjust = 1))
dev.off()

pdf("Fraction.CNV.All.pdf",4,5)
view.data <- as.data.frame(table(sample.info[cell.cnv.all,]$Sample_2))
colnames(view.data) <- c("Patient", "nCell_CNV")
table(sample.info[cell.cnv.all,]$Sample)

view.data$nCell_QC <- as.numeric(table(sample.info$Sample))
view.data$CNV_Precent <- round(view.data$nCell_CNV*100 / view.data$nCell_QC, 2)

multiplot(
   ggplot(view.data)+
      geom_bar(aes(x=Patient, y=nCell_QC, fill=Patient),
               stat = "identity")+
      geom_text(aes(x=Patient, y=nCell_QC+10, label=nCell_QC),size = 3)+
      scale_fill_manual(values = patient.color)+
      theme_bw()+theme(panel.grid = element_blank()),
   ggplot(view.data)+
      geom_bar(aes(x=Patient, y=CNV_Precent, fill=Patient),
               stat = "identity")+
      geom_text(aes(x=Patient, y=CNV_Precent+0.1, label=CNV_Precent),size = 3)+
      scale_fill_manual(values = patient.color)+
      theme_bw()+theme(panel.grid = element_blank()),
   cols = 1 
)
dev.off()

# 5. generating merged CNV Heatmap ----------------------------------------
# pre defined CNVs
CNV.info <- read.table("CNV.info.tsv")
colnames(CNV.info) <- c("Patient","chr","start","end","CN")
CNV.info$Patient <- sample.rename[CNV.info$Patient ]
CNV.info <- CNV.info[-8,]
CNV.info <- GRanges(CNV.info)

# results from CNV seq
SegCopy_bulk <- SegNorm[,1:10]
SegCopy_bulk[,] <- 2
colnames(SegCopy_bulk) <- c(paste0("A",c(1:6,8),"_CNV"),
                            paste0("B",c(1:3),"_CNV"))

bin.info.bed <- GRanges(bin.info)

SegCopy_bulk[countOverlaps(bin.info.bed,CNV.info[1])>0,]$A1_CNV <- 1.0
SegCopy_bulk[countOverlaps(bin.info.bed,CNV.info[2])>0,]$A2_CNV <- 1.0
SegCopy_bulk[countOverlaps(bin.info.bed,CNV.info[3])>0,]$A3_CNV <- 1.7
SegCopy_bulk[countOverlaps(bin.info.bed,CNV.info[4])>0,]$A4_CNV <- 1.7
SegCopy_bulk[countOverlaps(bin.info.bed,CNV.info[5])>0,]$A5_CNV <- 2.5
SegCopy_bulk[countOverlaps(bin.info.bed,CNV.info[6])>0,]$A5_CNV <- 2.5
SegCopy_bulk[countOverlaps(bin.info.bed,CNV.info[7])>0,]$A6_CNV <- 1.0
SegCopy_bulk[countOverlaps(bin.info.bed,CNV.info[8])>0,]$A8_CNV <- 2.4

bin.auto <- rownames(bin.info[!bin.info$CHR%in%c("chrX","chrY"),])
bin.X <- rownames(bin.info[bin.info$CHR%in%c("chrX"),])
bin.Y <- rownames(bin.info[bin.info$CHR%in%c("chrY"),])

SegCopy_bulk[bin.X,sex=="XX"] <- 2
SegCopy_bulk[bin.Y,sex=="XX"] <- 0
SegCopy_bulk[bin.X,sex=="XY"] <- 1
SegCopy_bulk[bin.Y,sex=="XY"] <- 1

# convert SegRef (from 50bwa to 150bwa)
SegRef_bulk.50bwa <- read.table("../../CNVSeq/data/bed_SegNorm_ref",
                                header = T,
                                row.names = paste("bin_",1:1048,sep = ''))

bin.info.50bwa <- SegRef_bulk.50bwa[,1:3]
bin.info.50bwa.bed <- GRanges(bin.info.50bwa)
SegRef_bulk.50bwa <- SegRef_bulk.50bwa[,-c(1:3)]
colnames(SegRef_bulk.50bwa) <- gsub("_DNA_","_CNV",colnames(SegRef_bulk.50bwa))
SegRef_bulk.50bwa <- SegRef_bulk.50bwa[,colnames(SegCopy_bulk)]

SegRef_bulk <- SegCopy_bulk
SegRef_bulk[,] <- 0

for (one in 1:nrow(SegRef_bulk)) {
   SegRef_bulk[one,] <- colMeans(SegRef_bulk.50bwa[countOverlaps(bin.info.50bwa.bed,
                                                                 bin.info.bed[one])>0,])
}

SegRef_bulk[bin.Y,sex=="XX"] <- 0

rm(bin.info.50bwa, bin.info.50bwa.bed, SegRef_bulk.50bwa)

# sort by location of first cnv 
sample.info.CNV <- sample.info[cell.cnv.3A,]
rownames(sample.info.CNV) %in% colnames(SegCopy.filter)
sample.info.CNV$Location <- as.numeric(lapply(apply(SegCopy.filter[,rownames(sample.info.CNV)]!=2,
                                                    2,
                                                    which),
                                   min))
cell.cnv.all <- rownames(sample.info.CNV)[order(sample.info.CNV$Location)]
cell.cnv.3A <- cell.cnv.all[cell.cnv.all%in%cell.cnv.3A]


# generate merged heatmap
col_fun = colorRamp2(c(1, 2, 3), c("blue", "white", "red"))
bin.info$CHR <- gsub("chr","",bin.info$CHR)
bin.info$CHR <- factor(bin.info$CHR, levels=unique(bin.info$CHR))

chr.color.2 <- rep(c("gray85","gray75"),12)
names(chr.color.2) <- unique(bin.info$CHR)


# sex chr relative value
SegCopy.filter.sex <- SegCopy.filter
cell.XX <- sample.info[sample.info$Sample_2 %in% names(which(sex=="XX")),]$Cell
cell.XY <- sample.info[sample.info$Sample_2 %in% names(which(sex=="XY")),]$Cell

SegCopy.filter.sex[bin.X,cell.XX] <- SegCopy.filter.sex[bin.X,cell.XX] / 2 * 2
SegCopy.filter.sex[bin.Y,cell.XX] <- 2
SegCopy.filter.sex[bin.X,cell.XY] <- SegCopy.filter.sex[bin.X,cell.XY] / 1 * 2
SegCopy.filter.sex[bin.Y,cell.XY] <- 2

cell.no.cnv <- names(which(colSums(SegCopy.filter.sex[bin.auto,]!=2)==0))
table(substr(cell.cnv.3A,5,6))

sample.info$is_CNV <- "normal"
sample.info[cell.cnv.3A,]$is_CNV <- "CNV"

ht_opt$TITLE_PADDING = unit(c(4, 4), "points")

for (one in levels(sample.info$Sample)) {
   cells <- grep(one,cell.cnv.all,value = T)
   
   if (length(cells) < 10) {
      cells.2 <- sample(grep(one,cell.no.cnv,value = T),10-length(cells))
      cells <- c(cells,cells.2)
      rm(cells.2)
   }
   
   point.cloud <- SegRef_bulk[,paste0(one,"_CNV")]
   point.cloud[point.cloud>4] <- 100000
   point.normal <- SegCopy_bulk[,paste0(one,"_CNV")]
   point.amp <- point.normal
   point.del <- point.normal
   point.normal[point.normal!=2] <- 100000
   point.amp[point.amp < 2.01] <- 100000
   point.del[point.del > 1.99] <- 100000

   ha.col <- HeatmapAnnotation(TE=anno_points(data.frame("Cloud"=point.cloud,
                                                         "normal"=point.normal,
                                                         "amp"=point.amp,
                                                         "del"=point.del), 
                                              size = unit(c(1.3,0.7,0.7,0.7), "mm"),
                                              gp = gpar(col = c("gray","black",
                                                                "red","blue")),
                                              ylim = c(0, 4)),
                               show_annotation_name = c(TE = TRUE),
                               annotation_name_gp = gpar(cex = 0.7),
                               annotation_name_side = "left",
                               annotation_name_rot = 0,
                               border = T)

   ht1 <- Heatmap(t(SegCopy.filter.sex[,cells]),
                  col = col_fun,
                  cluster_rows = F, cluster_columns = F,
                  show_row_names = F,show_column_names = F,
                  top_annotation = ha.col,
                  
                  row_split = rep("PBMC",10),
                  row_title_gp = gpar(fill = patient.color[sample.rename[one]],
                                      cex = 0.7),
                  row_gap = unit(0, "mm"), column_gap = unit(0, "mm"),
                  row_title_rot = 0,
                  
                  column_split = bin.info$CHR,
                  column_title_gp = gpar(fill = chr.color.2,
                                         cex = 0.7),
                  
                  
                  name="Copy Number",
                  border = T)
   
   if (one == "A1") {
      ht <- ht1
   }
   else {
      ht <- ht %v% ht1
   }

}
rm(ha.col,ht1)

lgd = Legend(labels = names(patient.color),
             title = "Sample",
             legend_gp = gpar(fill = patient.color))

pdf("Heatmap_CNV_mergeTE.All.2.pdf",10,10)
draw(ht, annotation_legend_list = list(lgd),
     annotation_legend_side = "right")
dev.off()


saveRDS(cell.cnv.all,"cell.cnv.all.rds")
saveRDS(cell.cnv.3A,"cell.cnv.3A.rds")
save.image("scDNA.RData")


# 6. generate supplementary table -----------------------------------------
saveRDS(sample.info,"sample.info.qc.rds")
sample.info.out <- sample.info.CNV

sample.info.out$Cell <- gsub("DNA_","",sample.info.out$Cell)
sample.info.out$Cell <- gsub("_FK.+?1a","",sample.info.out$Cell)
sample.info.out$Cell <- gsub(".","_",sample.info.out$Cell,fixed = T)

for (one in names(sample.rename)) {
   sample.info.out$Cell <- gsub(one,sample.rename[one],
                                sample.info.out$Cell,
                                fixed = T)
}

colnames(sample.info.out)
sample.info.out <- sample.info.out[,c(13,11)]
colnames(sample.info.out) <- c("Patient","Cell_ID")

sample.info.merge <- read.csv("../02.scRNA/Table S1. Summarize of scRNA-Seq.csv")
rownames(sample.info.merge) <- sample.info.merge$Cell_ID
sample.info.out$Cell_Type <- sample.info.merge[sample.info.out$Cell_ID,]$Cell_Type
table(sample.info.out$Cell_Type)

cells <- rownames(sample.info.out)

cnv.bin.list <- apply(SegCopy.filter[bin.auto, cells],2,function(x){return(which(x!=2))})
cnv.region.list <- lapply(cnv.bin.list,function(x){
   return(reduce(bin.info.bed[x]))
})

cnv.region.char <- lapply(cnv.region.list,function(x){
   return(paste(paste0(x@seqnames,":",start(x),"-",end(x)),collapse = "; "))
})
cnv.region.char <- unlist(cnv.region.char)

identical(names(cnv.region.char),rownames(sample.info.out))
sample.info.out$Region_of_SCNAs <- cnv.region.char

cnv.region.width <- lapply(cnv.region.list,function(x){
   return(sum(width(x)))
})
cnv.region.width <- unlist(cnv.region.width)
identical(names(cnv.region.width),rownames(sample.info.out))
sample.info.out$Size_of_SCNAs <- cnv.region.width

# rename cell
temp1 <- grep("P",sample.info.out$Cell_ID,value = T)
names(temp1) <- temp1
temp1 <- gsub("lym","plate",temp1)
temp2 <- grep("N",sample.info.out$Cell_ID,value = T)
names(temp2) <- temp2
temp2 <- gsub("lym","plate_1",temp2)

temp <- c(temp1,temp2)[sample.info.out$Cell_ID]
rm(temp1, temp2)

temp <- paste0(substr(temp,1,11),"_Cell_#",substr(temp,13,14))
sample.info.out$Cell_ID <- temp

write.csv(sample.info.out,"Table S2. Summarize of Cells with SCNAs.csv",row.names = F,col.names = T,quote = F)
