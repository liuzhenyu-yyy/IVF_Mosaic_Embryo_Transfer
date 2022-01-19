setwd("E:/LabWork/Project/fetal_genome/Merge/02.scRNA/")
source("E:/LabWork/code/MyFunction.R")

library(Seurat)
library(patchwork)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(dplyr)

# 1. load data ---------------------------------------------------------------
sample.ls <- gsub(".umi_counts.xls|.gene_tran_num.xls","",dir("data/")) %>% unique()

rc <- read.table(paste("data/RNA_A1_lym_1_FKDL210056829-1a.umi_counts.xls",sep = ""),
                 header = T,row.names =1)[,1:2]
colnames(rc) <- c("t1","t2")
sample.info <- data.frame()

for (one in sample.ls) {
  temp <- read.table(paste("data/",one,".umi_counts.xls",sep = ""),
                     header = T,row.names =1)
  rc <- cbind(rc, temp)
  temp <- read.table(paste("data/",one,".gene_tran_num.xls",sep = ""),
                     header = T,row.names =1)
  sample.info <- rbind(sample.info, temp)
}

rownames(sample.info) <- gsub("-",".",rownames(sample.info))
rc <- rc[ ,rownames(sample.info)]
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

dim(sample.info)
rc <- rc[,rownames(sample.info)]

write.table(sample.info,"sample.info.tsv",quote = F)
saveRDS(rc,"data/rc.merge.rds")

patient.color <- brewer.pal(10,"Set3")
names(patient.color) <- levels(sample.info$Sample_2)

pdf("01.gene_tran.pdf",6,3)
ggplot(data=sample.info)+
  geom_boxplot(aes(x=Sample_2,y=Gene_Number, color =Sample_2),
               outlier.alpha = 0)+
  geom_violin(aes(x=Sample_2,y=Gene_Number, fill =Sample_2),alpha=0.5)+
  scale_fill_manual(values = patient.color)+
  scale_color_manual(values = patient.color)+
  theme_bw()
ggplot(data=sample.info)+
  geom_boxplot(aes(x=Sample_2,y=log10(Transcript_Number), color =Sample_2),
               outlier.alpha = 0)+
  geom_violin(aes(x=Sample_2,y=log10(Transcript_Number), fill =Sample_2),alpha=0.5)+
  scale_fill_manual(values = patient.color)+
  scale_color_manual(values = patient.color)+
  theme_bw()
dev.off()

pdf("01.gene_tran.scatter.pdf",5,4)
ggplot(data=sample.info)+
  geom_point(aes(x=Transcript_Number,y=Gene_Number, color =Sample_2))+
  theme_bw()+
  scale_color_manual(values = patient.color)
dev.off()


aggregate(sample.info$Gene_Number,list(sample.info$Sample_2),mean) #577.7034 
aggregate(sample.info$Gene_Number,list(sample.info$Sample_2),median) # 300 

aggregate(sample.info$Transcript_Number,list(sample.info$Sample_2),mean) # 7192.025 
aggregate(sample.info$Transcript_Number,list(sample.info$Sample_2),median) # 599 

# 2. QC  ------------------------------------------------------------
dim(rc)

seu.obj <- CreateSeuratObject(rc, project = "fetal_genome",
                              assay = "RNA",
                              meta.data = sample.info)
seu.obj

seu.obj[["percent.mt"]] <- PercentageFeatureSet(seu.obj, pattern = "^MT-")



hist(seu.obj$nFeature_RNA,breaks = 60)
hist(seu.obj$percent.mt,breaks = 60)

pdf("02.QC.Scatter.pdf",10,4)
plot1 <- FeatureScatter(seu.obj, feature1 = "nFeature_RNA", feature2 = "percent.mt",
                        group.by = "Sample_2",cols = patient.color)
plot2 <- FeatureScatter(seu.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                        group.by = "Sample_2",cols = patient.color)
plot(plot1 + plot2)
plot1 = plot1 + geom_hline(yintercept = 30, color = "red")+
  geom_vline(xintercept = 300, color = "red")
plot2 = plot2 + geom_hline(yintercept = c(300), color = "red")
plot(plot1 + plot2)
dev.off()

pdf("02.QC.Scatter.2.pdf",4,3)
plot(plot1+theme_lzy()+coord_equal(40))
dev.off()


table(seu.obj$nFeature_RNA > 300 & seu.obj$percent.mt < 30)
seu.obj.qc <- subset(seu.obj, subset = nFeature_RNA > 300 & percent.mt < 30)

table(seu.obj.qc$Sample)
QC.tab <- round(table(seu.obj.qc$Sample) / table(sample.info$Sample),3)

seu.obj.qc <- NormalizeData(seu.obj.qc, normalization.method = "LogNormalize",
                            scale.factor = 10000)

seu.obj.qc <- FindVariableFeatures(seu.obj.qc, selection.method = "vst")

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seu.obj.qc), 10)

seu.obj.qc <- ScaleData(seu.obj.qc)

seu.obj.qc <- RunPCA(seu.obj.qc, features = VariableFeatures(object = seu.obj.qc))

DimPlot(seu.obj.qc, reduction = "pca")

ElbowPlot(seu.obj.qc,ndims = 30)
pc.use = 1:8

# t-SNE
seu.obj.qc <- RunUMAP(seu.obj.qc, dims = pc.use)

pdf("03.DR.raw.pdf",6,5)
DimPlot(seu.obj.qc, reduction = "umap",group.by = "Sample_2",
        cols = patient.color)+
  coord_equal(1.5)+theme_bw()+
  theme(panel.grid = element_blank())
dev.off()

pdf("03.DR.Quality.pdf",6,5)
FeaturePlot(seu.obj.qc, features = "percent.mt", 
            reduction = "umap",cols = c("gray","red"))+
  coord_equal(1.5)+theme_bw()+
  theme(panel.grid = element_blank())
FeaturePlot(seu.obj.qc, features = "nCount_RNA", 
            reduction = "umap",cols = c("gray","red"))+
  coord_equal(1.5)+theme_bw()+
  theme(panel.grid = element_blank())
FeaturePlot(seu.obj.qc, features = "nFeature_RNA", 
            reduction = "umap",cols = c("gray","red"))+
  coord_equal(1.5)+theme_bw()+
  theme(panel.grid = element_blank())
FeaturePlot(seu.obj.qc, features = "percent.mt", 
            reduction = "pca",cols = c("gray","red"))+
  coord_equal(1.2)+theme_bw()+
  theme(panel.grid = element_blank())
FeaturePlot(seu.obj.qc, features = "nCount_RNA", 
            reduction = "pca",cols = c("gray","red"))+
  coord_equal(1.2)+theme_bw()+
  theme(panel.grid = element_blank())
FeaturePlot(seu.obj.qc, features = "nFeature_RNA", 
            reduction = "pca",cols = c("gray","red"))+
  coord_equal(1.2)+theme_bw()+
  theme(panel.grid = element_blank())
dev.off()

# 3. identify cell type ---------------------------------------------------
# construct KNN graph
seu.obj.qc <- FindNeighbors(seu.obj.qc, dims = pc.use) 
seu.obj.qc <- FindClusters(seu.obj.qc, resolution = 1.0)

table(Idents(seu.obj.qc))

pdf("03.DR.cluster.pdf",6,5)
DimPlot(seu.obj.qc, reduction = "umap",label = T)+
  coord_equal(1.5)+theme_bw()+
  theme(panel.grid = element_blank())
DimPlot(seu.obj.qc, reduction = "pca",label = T)+
  coord_equal(1.2)+theme_bw()+
  theme(panel.grid = element_blank())
dev.off()

marker.list <- c("EPCAM","PECAM1","CDH1","PTPRC",
                 "CD3D","CD3G","CD4","CD8A", #T
                 "NCAM1","NKG7","KLRF1","XCL1","GZMA", # NK
                 "MS4A1", "CD79A","CD27",# B
                 "CD14","LYZ", #CD14+ Mono
                 "CST3","ITGAX", # DC
                 "FCGR3A", "MS4A7", # FCGR3A+ Mono
                 "CD68", # Monocytes
                 "MZB1", "IGHG1","IGLL1",# Plasma
                 "HBB","HBA1","HBA2",
                 "AIF1","CD14",# Macrophage
                 "PPBP" #megakaryocytes
                 ) 

pdf("03.Marker.pdf",6,5)
for (one in marker.list) {
  p <- FeaturePlot(seu.obj.qc, features = one, 
                   reduction = "umap",cols = c("gray","red"),order = T)
  p <- p+
    coord_equal(1.5)+theme_bw()+
    theme(panel.grid = element_blank())
  plot(p)
}
dev.off()

# run cell blast
write.table(rc[,colnames(seu.obj.qc)],"rc.merge.tsv",col.names = T,row.names = T,quote = F,sep = "\t")
cell.predict <- read.csv("pred_9e1107e7-314f-466e-9ef6-b7d9bb653e36.csv")

table(cell.predict$qid %in% colnames(seu.obj.qc))
rownames(cell.predict) <- cell.predict$qid
cell.predict <- cell.predict[colnames(seu.obj.qc),]
seu.obj.qc$cell_type_predict <- cell.predict$cell_ontology_class
seu.obj.qc$prdict_frac <- cell.predict$majority_frac.for.cell_ontology_class
seu.obj.qc$n_hits <- cell.predict$n_hits

table(seu.obj.qc$cell_type_predict)
table(is.na(seu.obj.qc$prdict_frac))

pdf("03.DR.predicted_cell_type.pdf",10,5)
DimPlot(seu.obj.qc, reduction = "umap",group.by = "cell_type_predict",
              )+
  coord_equal(1.5)+theme_bw()+
  theme(panel.grid = element_blank())
dev.off()
pdf("03.DR.predicted_frac.pdf",6,5)
FeaturePlot(seu.obj.qc, reduction = "umap",features = "prdict_frac",
            cols = c("gray","red"))+
  coord_equal(1.5)+theme_bw()+
  theme(panel.grid = element_blank())
dev.off()

DimPlot(seu.obj.qc, reduction = "umap",label = T)+
  coord_equal(1.5)+theme_bw()+
  theme(panel.grid = element_blank())

markers <- FindAllMarkers(seu.obj.qc,only.pos = T)

levels(seu.obj.qc)
new_id <- c("Monocytes","B","NK","T","T",
            "Monocytes","Monocytes","Unclassified","T","DC")
names(new_id) <- levels(seu.obj.qc)
seu.obj.qc <- RenameIdents(seu.obj.qc, new_id)

DimPlot(seu.obj.qc, reduction = "umap")
seu.obj.qc$cell_type <- Idents(seu.obj.qc)
seu.obj.qc$cell_type <- factor(seu.obj.qc$cell_type, levels = c("B","T","NK","Monocytes","DC","Unclassified"))

colors <- brewer.pal(5,"Set2")
colors <- c(colors,"gray")
names(colors) <- levels(seu.obj.qc$cell_type)

pdf("04.Cell_Type.pdf",6,5)
DimPlot(seu.obj.qc, reduction = "umap",
        group.by = "cell_type", cols = colors)+
  coord_equal(1.5)+theme_bw()+
  theme(panel.grid = element_blank())
dev.off()

view.data <- aggregate(seu.obj.qc$nFeature_RNA,list(seu.obj.qc$Sample),mean)
colnames(view.data) <- c("Patient","Average nGene")


# fraction of cell for each patient
plot.data <- as.data.frame(table(seu.obj.qc$Sample_2,
                                 seu.obj.qc$cell_type))
colnames(plot.data) <- c("Patient","Cell_type","Count")

plot.data$Fraction <- plot.data$Count / rep(table(seu.obj.qc$Sample))

plot.data$Cell_type


pdf("05.Fraction_Cell_type.pdf",5,3)
ggplot(plot.data,aes(x=Patient,y=Count,fill=Cell_type))+
  geom_bar(stat = 'identity')+
  # geom_text(aes(label =Count), 
  #           position = position_stack(), # 可以不设置该参数
  #           vjust = 0, hjust=0.5)+
  theme_bw()+theme(panel.grid = element_blank())+
  scale_fill_manual(values = colors)
ggplot(plot.data,aes(x=Patient,y=Fraction,fill=Cell_type))+
  geom_bar(stat = 'identity')+
  theme_bw()+theme(panel.grid = element_blank())+
  scale_fill_manual(values = colors)
dev.off()

sample.info.qc <- seu.obj.qc@meta.data

pdf("05.gene_tran.QC.pdf",5,2.5)
ggplot(data=sample.info.qc)+
  geom_boxplot(aes(x=Sample_2,y=Gene_Number, color =Sample_2),
               outlier.alpha = 0)+
  geom_violin(aes(x=Sample_2,y=Gene_Number, fill =Sample_2),alpha=0.5)+
  scale_fill_manual(values = patient.color)+
  scale_color_manual(values = patient.color)+
  theme_bw()
ggplot(data=sample.info.qc)+
  geom_boxplot(aes(x=Sample_2,y=log10(Transcript_Number), color =Sample_2),
               outlier.alpha = 0)+
  geom_violin(aes(x=Sample_2,y=log10(Transcript_Number), fill =Sample_2),alpha=0.5)+
  scale_fill_manual(values = patient.color)+
  scale_color_manual(values = patient.color)+
  theme_bw()
dev.off()

table(seu.obj.qc$Sample)

pdf("05.UMAP.one.patient.pdf",4,3)
for (one in levels(sample.info$Sample)) {
  
  seu.obj.qc.sub <- subset(seu.obj.qc, subset = Sample == one)
  if(ncol(seu.obj.qc.sub) < 15){
    next
  }
  seu.obj.qc.sub <- seu.obj.qc.sub %>% 
    NormalizeData()  %>% 
    FindVariableFeatures() %>% 
    ScaleData() %>%
    RunPCA(npcs = 10) %>% 
    RunUMAP(dims = 1:4)
  
  p <- DimPlot(seu.obj.qc.sub,cols = colors)+ggtitle(sample.rename[one])
  plot(p)
}
dev.off()
rm(seu.obj.qc.sub)

# 4. CNV profiles ---------------------------------------------------------
cell.cnv <- readRDS("../01.scDNA/cell.cnv.3A.rds")
cell.cnv <- gsub("DNA","RNA",cell.cnv)
cell.cnv <- gsub("_FKDL.+?1a","",cell.cnv)

cell.RNA <- rownames(sample.info)
names(cell.RNA) <- cell.RNA
cell.RNA <- gsub("_FKDL.+?1a","",cell.RNA)
cell.RNA <- gsub(".","_",cell.RNA,fixed = T)

table(cell.cnv %in% cell.RNA) # 20
length(cell.RNA[cell.RNA%in%cell.cnv]) # 20
cell.cnv <- cell.RNA[cell.RNA%in%cell.cnv]
cell.cnv <- names(cell.cnv)

table(cell.cnv %in% colnames(seu.obj.qc)) # 5 of 15

seu.obj.qc$CNV <- as.character(seu.obj.qc$cell_type)
seu.obj.qc@meta.data[!colnames(seu.obj.qc) %in% cell.cnv,]$CNV <- "normal"
table(seu.obj.qc$CNV)

seu.obj.qc$CNV <- factor(seu.obj.qc$CNV, levels = c("normal","B","T","NK","Monocytes","DC","Unclassified"))

pdf("05.CNV.QC.pdf",5,4)
DimPlot(seu.obj.qc,group.by = "CNV",
        cols = c(colors,
                 "normal"="gray90"),order = T)+
  coord_equal(1.5)+theme_bw()+
  theme(panel.grid = element_blank())
dev.off()


# 5. add CNV cells --------------------------------------------------------
cell.select <- unique(c(colnames(seu.obj.qc),cell.cnv))

seu.obj.cnv <- CreateSeuratObject(rc[,cell.select], project = "fetal_genome",
                                  assay = "RNA",
                                  meta.data = sample.info)

seu.obj.cnv <- subset(seu.obj.cnv,subset = nFeature_RNA>200)

seu.obj.cnv <- NormalizeData(seu.obj.cnv, normalization.method = "LogNormalize",
                            scale.factor = 10000)

seu.obj.cnv <- FindVariableFeatures(seu.obj.cnv, selection.method = "vst")

seu.obj.cnv <- ScaleData(seu.obj.cnv)

seu.obj.cnv <- RunPCA(seu.obj.cnv, features = VariableFeatures(object = seu.obj.cnv))

DimPlot(seu.obj.cnv, reduction = "pca")

ElbowPlot(seu.obj.cnv,ndims = 30)
pc.use = 1:7

# UMAP
seu.obj.cnv <- RunUMAP(seu.obj.cnv, dims = pc.use)

pdf("06.DR.raw.pdf",4,3)
DimPlot(seu.obj.cnv, reduction = "umap",group.by = "Sample_2",
        cols = patient.color,pt.size = 0.5)+
  coord_equal(1)+theme_bw()+
  theme(panel.grid = element_blank())
dev.off()

pdf("06.DR.Quality.pdf",6,5)
FeaturePlot(seu.obj.cnv, features = "nCount_RNA", 
            reduction = "umap",cols = c("gray","red"))+
  coord_equal(1)+theme_bw()+
  theme(panel.grid = element_blank())
FeaturePlot(seu.obj.cnv, features = "nFeature_RNA", 
            reduction = "umap",cols = c("gray","red"))+
  coord_equal(1)+theme_bw()+
  theme(panel.grid = element_blank())
dev.off()

seu.obj.cnv <- FindNeighbors(seu.obj.cnv)
seu.obj.cnv <- FindClusters(seu.obj.cnv,resolution = 1)

table(Idents(seu.obj.cnv))

pdf("06.DR.cluster.pdf",6,5)
DimPlot(seu.obj.cnv, reduction = "umap",label = T)+
  coord_equal(1)+theme_bw()+
  theme(panel.grid = element_blank())
dev.off()

levels(seu.obj.cnv)
new_id <- c("Monocytes","Unclassified","B","T","T",
            "NK","Monocytes","T","DC")
names(new_id) <- levels(seu.obj.cnv)
seu.obj.cnv <- RenameIdents(seu.obj.cnv, new_id)
seu.obj.cnv$cell_type <- Idents(seu.obj.cnv)

seu.obj.cnv$cell_type <- factor(seu.obj.cnv$cell_type,
                                levels = c("B","T","NK","Monocytes","DC","Unclassified"))

p <- DimPlot(seu.obj.cnv, reduction = "umap",
             group.by = "cell_type", cols = colors,order = T)+
  coord_equal(1.2)+theme_bw()+
  theme(panel.grid = element_blank())

cell.b <- CellSelector(p)
cell.amb <- CellSelector(p)
cell.dc <- CellSelector(p)

seu.obj.cnv@meta.data[cell.b,]$cell_type <- "B"
seu.obj.cnv@meta.data[cell.amb,]$cell_type <- "Unclassified"
seu.obj.cnv@meta.data[cell.dc,]$cell_type <- "DC"

pdf("06.Cell_Type_CNV.pdf",4,3)
DimPlot(seu.obj.cnv, reduction = "umap",
        group.by = "cell_type", cols = colors,order = T,pt.size = 0.5)+
  coord_equal(1)+theme_bw()+
  theme(panel.grid = element_blank())
dev.off()

# fraction of cell for each patient
plot.data <- as.data.frame(table(seu.obj.cnv$Sample_2,
                                 seu.obj.cnv$cell_type))
colnames(plot.data) <- c("Patient","Cell_type","Count")

plot.data$Fraction <- plot.data$Count / rep(table(seu.obj.qc$Sample))

plot.data$Cell_type

pdf("06.Fraction_Cell_type.pdf",5,3)
ggplot(plot.data,aes(x=Patient,y=Count,fill=Cell_type))+
  geom_bar(stat = 'identity')+
  # geom_text(aes(label =Count), 
  #           position = position_stack(), # 可以不设置该参数
  #           vjust = 0, hjust=0.5)+
  theme_bw()+theme(panel.grid = element_blank())+
  scale_fill_manual(values = colors)
ggplot(plot.data,aes(x=Patient,y=Fraction,fill=Cell_type))+
  geom_bar(stat = 'identity')+
  theme_bw()+theme(panel.grid = element_blank())+
  scale_fill_manual(values = colors)
dev.off()

pdf("06.UMAP.one.patient.pdf",4,3)
for (one in levels(sample.info$Sample)) {
  
  seu.obj.cnv.sub <- subset(seu.obj.cnv, subset = Sample == one)
  if(ncol(seu.obj.cnv.sub) < 15){
    next
  }
  seu.obj.cnv.sub <- seu.obj.cnv.sub %>% 
    NormalizeData()  %>% 
    FindVariableFeatures() %>% 
    ScaleData() %>%
    RunPCA(npcs = 10) %>% 
    RunUMAP(dims = 1:4)
  
  p <- DimPlot(seu.obj.cnv.sub,cols = colors)+ggtitle(sample.rename[one])
  plot(p)
}
dev.off()
rm(seu.obj.cnv.sub)

# CNV 
seu.obj.cnv$CNV <- as.character(seu.obj.cnv$cell_type)
seu.obj.cnv@meta.data[!colnames(seu.obj.cnv) %in% cell.cnv,]$CNV <- "normal"
table(seu.obj.cnv$CNV)

seu.obj.cnv$CNV <- factor(seu.obj.cnv$CNV, levels = c("normal","B","T","NK","Monocytes","DC","Unclassified"))

pdf("06.CNV.QC.pdf",4,3)
DimPlot(seu.obj.cnv,group.by = "CNV",
        cols = c(colors,
                 "normal"="gray90"),order = T)+
  coord_equal(1)+theme_bw()+
  theme(panel.grid = element_blank())
dev.off()

sample.info.qc <- seu.obj.qc@meta.data

pdf("06.gene_tran.QC.pdf",5,2.5)
ggplot(data=seu.obj.cnv@meta.data)+
  geom_boxplot(aes(x=cell_type,y=Gene_Number, color =cell_type),
               outlier.alpha = 0)+
  geom_violin(aes(x=cell_type,y=Gene_Number, fill =cell_type),alpha=0.5)+
  scale_fill_manual(values = colors)+
  scale_color_manual(values = colors)+
  coord_equal(0.0005)+
  ylab("nGene")+
  theme_bw()
ggplot(data=seu.obj.cnv@meta.data)+
  geom_boxplot(aes(x=cell_type,y=log10(Transcript_Number), color =cell_type),
               outlier.alpha = 0)+
  geom_violin(aes(x=cell_type,y=log10(Transcript_Number), fill =cell_type),alpha=0.5)+
  scale_fill_manual(values = colors)+
  scale_color_manual(values = colors)+
  coord_equal(0.6)+
  ylab("log10(nUMI)")+
  theme_bw()
dev.off()

pdf("06.Marker.pdf",4,3)
for (one in marker.list) {
  p <- FeaturePlot(seu.obj.cnv, features = one, 
                   reduction = "umap",cols = c("gray","red"),order = F)
  p <- p+
    coord_equal(1)+theme_bw()+
    theme(panel.grid = element_blank())+
    theme(plot.title = element_text(hjust = 0.5,face = "italic",family = ))
  plot(p)
}
dev.off()

marker.list.select <- c("PTPRC","CD3D","CD3G","NKG7","KLRF1",
                        "MS4A1", "CD79A","CD14","LYZ", #CD14+ Mono
                        "FCGR3A", "MS4A7", # FCGR3A+ Mono
                        "CD68", # Monocytes
                        "CST3","ITGAX" # DC
) 


plist <- VlnPlot(seu.obj.cnv,features = marker.list.select,
        cols = colors,group.by = "cell_type",ncol = 3,combine = F)
plist <- lapply(plist, function(x){
  return(x +theme(axis.text.x = element_blank()) + 
           theme(legend.position = "none") + 
           theme(axis.ticks.x = element_blank())+
           theme(plot.title = element_text(hjust = 0.5,face = "italic"))+
           theme(axis.title.x = element_blank())+
           theme(axis.title.y = element_blank()))
})

pdf("06.Violin.marker.pdf",12,5)
multiplot(plotlist = plist,cols = 5,layout = "row")
dev.off()

rm(p,plot.data, plot1,plot2,p1,temp)
save.image("scRNA.RData")


# 6. output supplement table ----------------------------------------------

sample.info.RNA <- seu.obj@meta.data
colnames(sample.info.RNA)
sample.info.RNA[,c("orig.ident","nCount_RNA","nFeature_RNA","Sample")] <- NULL
sample.info.RNA$percent.mt <- round(sample.info.RNA$percent.mt,2)

sample.info.RNA$Cell <- gsub("RNA_","",sample.info.RNA$Cell)
sample.info.RNA$Cell <- gsub("_FK.+?1a","",sample.info.RNA$Cell)
sample.info.RNA$Cell <- gsub(".","_",sample.info.RNA$Cell,fixed = T)

for (one in names(sample.rename)) {
  sample.info.RNA$Cell <- gsub(one,sample.rename[one],
                               sample.info.RNA$Cell,
                              fixed = T)
}

for (one in 1:9) {
  sample.info.RNA$Cell <- gsub(paste0("_",one,"$"),paste0("_0",one),
                               sample.info.RNA$Cell)
}

sample.info.RNA$Pass_RNA_QC <- "No"
sample.info.RNA[colnames(seu.obj.cnv),]$Pass_RNA_QC="Yes"
table(sample.info.RNA$Pass_RNA_QC )

sample.info.RNA$Cell_Type <- ""
sample.info.RNA[colnames(seu.obj.cnv),]$Cell_Type=as.character(seu.obj.cnv$cell_type)
table(sample.info.RNA$Cell_Type )

# write renamed rc
temp <- as.data.frame(seu.obj@assays$RNA@counts)
identical(colnames(temp),rownames(sample.info.RNA))
colnames(temp) <- sample.info.RNA$Cell
temp$Gene <- rownames(temp)
temp <- temp[,c(2206,1:2205)]
write.csv(temp,"read.count.csv",row.names = F,quote = F)
rm(temp)

rownames(sample.info.RNA) <- sample.info.RNA$Cell


# read DNA info
sample.info.DNA <- readRDS("E:/LabWork/Project/fetal_genome/Merge/01.scDNA/sample.info.qc.rds")
cell.qc.dna <- sample.info.DNA$Cell
sample.info.DNA <- readRDS("E:/LabWork/Project/fetal_genome/Merge/01.scDNA/data/sample.info.rds")

sample.info.DNA$Pass_DNA_QC <- "no"
sample.info.DNA[cell.qc.dna,]$Pass_DNA_QC <- "Yes"
table(sample.info.DNA$Pass_DNA_QC)

sample.info.DNA$Cell <- gsub("DNA_","",sample.info.DNA$Cell)
sample.info.DNA$Cell <- gsub("_FK.+?1a","",sample.info.DNA$Cell)
sample.info.DNA$Cell <- gsub(".","_",sample.info.DNA$Cell,fixed = T)

for (one in names(sample.rename)) {
  sample.info.DNA$Cell <- gsub(one,sample.rename[one],
                               sample.info.DNA$Cell,
                              fixed = T)
}
rownames(sample.info.DNA) <- sample.info.DNA$Cell
table(sample.info.RNA$Cell %in% sample.info.DNA$Cell)
table(sample.info.DNA$Cell %in% sample.info.RNA$Cell)

# add missing values
cell.empty <- sample.info.DNA$Cell[!sample.info.DNA$Cell %in% sample.info.RNA$Cell]
temp <- sample.info.RNA[1:3,]
rownames(temp) <- cell.empty
temp$Gene_Number <- 0
temp$Transcript_Number <- 0
temp$Cell <- cell.empty
temp$Sample_2 <- "P03"
temp$percent.mt <- 0
temp$Pass_RNA_QC <- "No" 
temp$Cell_Type <- "" 
sample.info.RNA <- rbind(sample.info.RNA, temp)

# merge sample information
sample.info.RNA <- sample.info.RNA[rownames(sample.info.DNA),]

colnames(sample.info.RNA)
sample.info.RNA[,c("Cell","Sample_2")] <- NULL
sample.info.RNA <- sample.info.RNA[,c(2,1,3,4,5)]
colnames(sample.info.DNA)
sample.info.DNA[,c("Min","X25th","Median","X75th","Max","Sample")] <- NULL
sample.info.DNA <- sample.info.DNA[,c(7,6,1,2,3,4,5,8)]

colnames(sample.info.RNA) <- c("RNA_UMI_Number",
                               "RNA_Gene_Number",
                               "RNA_Mitochondrial_Ratio",
                               "Pass_RNA_QC","Cell_Type")
colnames(sample.info.DNA) <- c("Patient",
                               "Cell_ID",
                               "DNA_Unique_Reads_Number",
                               "DNA_Bin_Number",
                               "Average_Reads_per_Bin",
                               "SD_of_Reads_per_Bin",
                               "CV",
                               "Pass_DNA_QC")

sample.info.merge <- cbind(sample.info.DNA, sample.info.RNA)

temp1 <- grep("P",sample.info.merge$Cell_ID,value = T)
names(temp1) <- temp1
temp1 <- gsub("lym","plate",temp1)
temp2 <- grep("N",sample.info.merge$Cell_ID,value = T)
names(temp2) <- temp2
temp2 <- gsub("lym","plate_1",temp2)

temp <- c(temp1,temp2)[sample.info.merge$Cell_ID]
rm(temp1, temp2)

temp <- paste0(substr(temp,1,11),"_Cell_#",substr(temp,13,14))
sample.info.merge$Cell_ID <- temp
saveRDS(temp, "cell.rename.rds")
write.csv(sample.info.merge,"Table S1. Summarize of Sequenced Cells.csv",row.names = F,col.names = T,quote = F)