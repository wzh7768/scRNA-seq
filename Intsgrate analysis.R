library(Seurat)
library(cowplot)
library(dplyr)
library(pheatmap)
library(monocle)
library(data.table)
library(ROTS)
library(clusterProfiler)
library(org.Mm.eg.db)
library(forcats)
library(ggplot2)
library(GOplot)
library(stringr)
library(GSEABase)
library(enrichplot)
aud1 <- read.csv('/home/wzh/Desktop/aud/GSE140883_P10_1.tsv',header=T, row.names=1, sep="\t") 
aud_1_1 <- read.csv('/home/wzh/Desktop/aud/GSE140883_P10_2.tsv',header=T, row.names=1, sep="\t")
aud1 <- t(aud1)
aud_1_1 <- t(aud_1_1)
aud_1 <- cbind(aud1,aud_1_1)
aud2 <- read.csv('/home/wzh/Desktop/aud/GSE140883_P15NR_1.tsv',header=T, row.names=1, sep="\t") 
aud_2_1 <- read.csv('/home/wzh/Desktop/aud/GSE140883_P15NR_2.tsv',header=T, row.names=1, sep="\t")
aud2 <- t(aud2)
aud_2_1 <- t(aud_2_1)
aud_2 <- cbind(aud2,aud_2_1)
aud3 <- read.csv('/home/wzh/Desktop/aud/GSE140883_P20NR_1.tsv',header=T, row.names=1, sep="\t") 
aud_3_1 <- read.csv('/home/wzh/Desktop/aud/GSE140883_P20NR_2.tsv',header=T, row.names=1, sep="\t")
aud3 <- t(aud3)
aud_3_1 <- t(aud_3_1)
aud_3 <- cbind(aud3,aud_3_1)
#Find mt's, compute the percent mt, and drop them from the raw data
drop_mt <- function(data) {
  mts <- grep(pattern = "^mt.", x = rownames(x = data), value = TRUE)
  percent.mt <- Matrix::colSums(data[mts, ])/Matrix::colSums(data)
  mt.index <- grep(pattern = "^mt.", x = rownames(x = data), value = FALSE)
  data_1 <- data[-mt.index,]
  return(data_1)
}
aud_1_dropmt <- drop_mt(aud_1)
dim(aud_1_dropmt)
aud_2_dropmt <- drop_mt(aud_2)
dim(aud_2_dropmt)
aud_3_dropmt <- drop_mt(aud_3)
dim(aud_3_dropmt)
# Create the Seurat object with all the data (unfiltered)
aud_P10_NR <- CreateSeuratObject(counts = aud_1_dropmt,project = "P10_NR")
aud_P15_NR <- CreateSeuratObject(counts = aud_2_dropmt,project = "P15_NR")
aud_P20_NR <- CreateSeuratObject(counts = aud_3_dropmt,project = "P20_NR")
#Filter cells 
audP10NR_filtered <- subset(x=aud_P10_NR, subset = nCount_RNA < 15000 & nFeature_RNA > 500)
audP15NR_filtered <- subset(x=aud_P15_NR, subset = nCount_RNA < 15000 & nFeature_RNA > 500)
audP20NR_filtered <- subset(x=aud_P20_NR, subset = nCount_RNA < 15000 & nFeature_RNA > 500)
#查看过滤后数据质量
plot <- FeatureScatter(audP20NR_filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot(plot)
VlnPlot(audP20NR_filtered, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, pt.size = 0)
#Seurat包整合分析(用于同一样本不同批次)
object.list <- c(audP10NR_filtered,audP15NR_filtered,audP20NR_filtered)
object.list <- lapply(X=object.list,FUN = function(x) {
  x <- NormalizeData(x,normalization.method = "LogNormalize", 
                     scale.factor = 10000)
  x <- FindVariableFeatures(object = x, selection.method = "vst", nfeatures = 2000)
})
sce.anchors <- FindIntegrationAnchors(object.list = object.list,dims=1:20)
sce.Int <- IntegrateData(anchorset = sce.anchors,dims=1:20)
#简单数据叠加(用于不同样本整合)
sce.comb <- merge(audP10NR_filtered,c(audP15NR_filtered,audP20NR_filtered),add.cell.ids = c("P10_NR","P15_NR","P20_NR"))
#归一化降维聚类
sce=sce.comb
sce <- NormalizeData(sce, normalization.method = "LogNormalize", 
                     scale.factor = 10000)
GetAssay(sce,assay = "RNA")
sce <- FindVariableFeatures(sce, 
                            selection.method = "vst", nfeatures = 2000) 
sce <- ScaleData(sce) 
sce <- RunPCA(object = sce, pc.genes = VariableFeatures(sce)) 
DimHeatmap(sce, dims = 1:12, cells = 100, balanced = TRUE)
ElbowPlot(sce)
sce <- FindNeighbors(sce, dims = 1:15)
sce <- FindClusters(sce, resolution = 0.2)
table(sce@meta.data$RNA_snn_res.0.2) 
sce <- FindClusters(sce, resolution = 0.8)
table(sce@meta.data$RNA_snn_res.0.8)
#看不同resolution分群细胞数的差异
library(gplots)
tab.1=table(sce@meta.data$RNA_snn_res.0.2,sce@meta.data$RNA_snn_res.0.8) 
balloonplot(tab.1)
#不同resolution条件下的分群情况
res.used <- seq(0.1,1,by=0.2)
res.used
# Loop over and perform clustering of different resolutions 
for(i in res.used){
  sce <- FindClusters(object = sce, verbose = T, resolution = res.used)
}
# Make plot 
library(clustree)
clus.tree.out <- clustree(sce) +
  theme(legend.position = "bottom") + 
  scale_color_brewer(palette = "Set1") +
  scale_edge_color_continuous(low = "grey80", high = "red")

clus.tree.out
#tSNE聚类
set.seed(1234)
sce <- RunTSNE(object = sce, dims = 1:15, do.fast = TRUE)
DimPlot(sce,reduction = "tsne",label=T)
phe=data.frame(cell=rownames(sce@meta.data),
               cluster =sce@meta.data$seurat_clusters)
head(phe)
table(phe$cluster)
tsne_pos=Embeddings(sce,'tsne') 
DimPlot(sce,reduction = "tsne",group.by ='orig.ident')
DimPlot(sce,reduction = "tsne",label=T,split.by ='orig.ident')
head(phe)
table(phe$cluster)
head(tsne_pos) 
dat=cbind(tsne_pos,phe)
pro='first'
save(dat,file=paste0(pro,'_for_tSNE.pos.Rdata')) 
load(file=paste0(pro,'_for_tSNE.pos.Rdata')) 
library(ggplot2)
p=ggplot(dat,aes(x=tSNE_1,y=tSNE_2,color=cluster))+geom_point(size=0.95)
p=p+stat_ellipse(data=dat,aes(x=tSNE_1,y=tSNE_2,fill=cluster,color=cluster),
                 geom = "polygon",alpha=0.2,level=0.9,type="t",linetype = 2,show.legend = F)+coord_fixed()
print(p) 
theme= theme(panel.grid =element_blank()) + ## 删去网格
  theme(panel.border = element_blank(),panel.background = element_blank()) + ## 删去外层边框
  theme(axis.line = element_line(size=1, colour = "black")) 
p=p+theme+guides(colour = guide_legend(override.aes = list(size=5)))
print(p)
ggplot2::ggsave(filename = paste0(pro,'_tsne_res0.8.pdf'))
#UMAP聚类
set.seed(123)
sce <- RunUMAP(object = sce, dims = 1:15, do.fast = TRUE)
DimPlot(sce,reduction = "umap",label=T)
phe=data.frame(cell=rownames(sce@meta.data),
               cluster =sce@meta.data$seurat_clusters)
head(phe)
table(phe$cluster)
tsne_pos=Embeddings(sce,'umap') 
DimPlot(sce,reduction = "umap",group.by ='orig.ident')
DimPlot(sce,reduction = "umap",label=T,split.by ='orig.ident')
head(phe)
table(phe$cluster)
head(tsne_pos) 
dat=cbind(tsne_pos,phe)
pro='first'
save(dat,file=paste0(pro,'_for_tSNE.pos.Rdata')) 
load(file=paste0(pro,'_for_tSNE.pos.Rdata')) 
library(ggplot2)
p=ggplot(dat,aes(x=UMAP_1,y=UMAP_2,color=cluster))+geom_point(size=0.95)
p=p+stat_ellipse(data=dat,aes(x=UMAP_1,y=UMAP_2,fill=cluster,color=cluster),
                 geom = "polygon",alpha=0.2,level=0.9,type="t",linetype = 2,show.legend = F)+coord_fixed()
print(p) 
theme= theme(panel.grid =element_blank()) + ## 删去网格
  theme(panel.border = element_blank(),panel.background = element_blank()) + ## 删去外层边框
  theme(axis.line = element_line(size=1, colour = "black")) 
p=p+theme+guides(colour = guide_legend(override.aes = list(size=5)))
print(p)
ggplot2::ggsave(filename = paste0(pro,'_umap_res0.8.pdf'))
#查看不同群的基因信息
#plot1 <- FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#CombinePlots(plots = list(plot1, plot2))
plot(plot2)
ggplot2::ggsave(filename = paste0(pro,'_CombinePlots.pdf'))
#查看不同群的基因表达信息
#VlnPlot(sce, features = c("percent.ribo", "percent.mt"), ncol = 2)
#ggplot2::ggsave(filename = paste0(pro,'_mt-and-ribo.pdf'))
VlnPlot(sce, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2,pt.size = 0)
ggplot2::ggsave(filename = paste0(pro,'_counts-and-feature.pdf'))
#VlnPlot(sce, features = c("percent.ribo", "nCount_RNA"), ncol = 2)
table(sce@meta.data$seurat_clusters)
#找到每一群的Markergene并绘制feature plot和vlnplot
for( i in unique(sce@meta.data$seurat_clusters) ){
  markers_df <- FindMarkers(object = sce, ident.1 = i, min.pct = 0.25)
  print(x = head(markers_df))
  markers_genes = rownames(head(x = markers_df, n = 5))
  VlnPlot(object = sce, features =markers_genes,log =T,pt.size = 0 )
  ggsave(filename=paste0(pro,'_VlnPlot_subcluster_',i,'_sce.markers_heatmap.pdf'))
  FeaturePlot(object = sce, features=markers_genes )
  ggsave(filename=paste0(pro,'_FeaturePlot_subcluster_',i,'_sce.markers_heatmap.pdf'))
}
#亚群依据生物学背景命名
FeaturePlot(object = sce, 
            features ='Cx3cr1', 
            cols= c("grey", "blue"), 
            reduction = "tsne")
head(as.numeric(as.character(sce@active.ident)))
tmp=sce@meta.data
a=read.table('/home/wzh/Desktop/SCN.txt')
labers=a[match(as.numeric(as.character(sce@active.ident)),a[,1]),2]
sce <- AddMetaData(object = sce, 
                    metadata = labers, 
                    col.name = 'labers')
colP<-c('green4', 
        'pink', 
        '#FF7F00', 
        'orchid', 
        '#99c9fb', 
        'dodgerblue2', 
        'yellow'
)
table(labers)
head(labers)
labers=as.factor(labers)
head(labers)
sce@meta.data$labers=labers 
DimPlot(sce, group.by = 'labers',
        cols =  colP,
        label = T,label.size=3,pt.size=0.6,reduction = "tsne")
#提取感兴趣的亚群进行分析(Excitatory Neuron)
sce_exci = subset(sce,idents =c(0,3,6,14,12,10,1,4,2))
sce_exci <- FindVariableFeatures(sce_exci, 
                            selection.method = "vst", nfeatures = 2000) 
sce_exci <- ScaleData(sce_exci) 
sce_exci <- RunPCA(object = sce_exci, pc.genes = VariableFeatures(sce_exci))
sce_exci <- FindNeighbors(sce_exci, dims = 1:15)
sce_exci <- FindClusters(sce_exci, resolution = 0.8)
sce_exci <- RunTSNE(object = sce_exci, dims = 1:15, do.fast = TRUE)
DimPlot(sce_exci,reduction = "tsne",label=T)
DimPlot(sce_exci,reduction = "tsne",group.by ='orig.ident')
FeaturePlot(object = sce_exci, 
            features ='Abcc5', 
            cols= c("grey", "red"), 
            reduction = "tsne")
VlnPlot(object = sce_exci, features ="Gad1",group.by = 'orig.ident',log =T,pt.size = 0 )
#提取感兴趣的亚群进行分析(Inhibitory Neuron)
sce_inhi = subset(sce,idents =c(5,7))
sce_inhi <- FindVariableFeatures(sce_inhi, 
                                 selection.method = "vst", nfeatures = 2000) 
sce_inhi <- ScaleData(sce_inhi) 
sce_inhi <- RunPCA(object = sce_inhi, pc.genes = VariableFeatures(sce_inhi))
sce_inhi <- FindNeighbors(sce_inhi, dims = 1:15)
sce_inhi <- FindClusters(sce_inhi, resolution = 0.5)
sce_inhi <- RunTSNE(object = sce_inhi, dims = 1:15, do.fast = TRUE)
DimPlot(sce_inhi,reduction = "tsne",label=T)
DimPlot(sce_inhi,reduction = "tsne",group.by ='orig.ident')
FeaturePlot(object = sce_inhi, 
            features ='Mapt', 
            cols= c("grey", "red"), 
            reduction = "tsne")
#绘制基因表达量的箱型图
data <- sce_inhi@assays$RNA@data[c("Gad1","Gad2","Pvalb"),]
data <- as.data.frame(data)
P10 <- data[,1:622]
P10 <- as.data.frame(t(P10))
P10['type'] <- 'P10'
P15 <- data[,623:1163]
P15 <- as.data.frame(t(P15))
P15['type'] <- 'P15'
P20 <- data[,1164:1604]
P20 <- as.data.frame(t(P20))
P20['type'] <- 'P20'
Pcom <- rbind(P10,P15)
Pcombine <- rbind(Pcom,P20)
P_Gad1 <- Pcombine[,c(1,4)]
P_Gad1 <- P_Gad1[P_Gad1$Gad1!=0,]
ggboxplot(P_Gad1,"type","Gad1",bxp.errorbar=T,fill = "type",palette = c("#00AFBB", "#E7B800", "#FC4E07"))
P_Gad2 <- Pcombine[,c(2,4)]
P_Gad2 <- P_Gad2[P_Gad2$Gad2!=0,]
ggboxplot(P_Gad2,"type","Gad2",bxp.errorbar=T,fill = "type",palette = c("#00AFBB", "#E7B800", "#FC4E07"))
P_Pvalb <- Pcombine[,c(3,4)]
P_Pvalb <- P_Pvalb[P_Pvalb$Pvalb!=0,]
ggboxplot(P_Pvalb,"type","Pvalb",bxp.errorbar=T,fill = "type",palette = c("#00AFBB", "#E7B800", "#FC4E07"))
#除去异常值函数
outlierKD <- function(dt,var){
  var_name <- eval(substitute(var),eval(dt))
  tot <- sum(!is.na(var_name))
  na1 <- sum(is.na(var_name))
  m1 <- mean(var_name,na.rm = T)
  par(mfrow=c(2,2),oma=c(0,0,3,0))
  boxplot(var_name,main="with outliers")
  hist(var_name,main="with outliers",xlab=NA,ylab=NA)
  outlier <- var_name[var_name>5]
  mo <- mean(outlier)
  var_name <- ifelse(var_name %in% outlier,NA,var_name)
  boxplot(var_name,main="without outliers")
  hist(var_name,main="without outliers",xlab=NA,ylab=NA)
  title("outlier check",outer=T)
  na2 <- sum(is.na(var_name))
  cat("outliers identified:",na2 - na1,"\n")
  cat("proportion (%) of outliers:",round((na2 - na1)/tot*100,1),"\n")
  cat("mean of the outliers:",round(mo,2),"\n")
  m2 <- mean(var_name,na.rm = T)
  cat("mean without removing outliers:",round(m1,2),"\n")
  cat("mean if we remove outliers:",round(m1,2),"\n")
  response <- readline(prompt = "Do you want to remove outliers and to replace with NA?[yes/no]:")
  if(response=="y" | response=="yes"){
    dt[as.character(substitute(var))] <- invisible(var_name)
    assign(as.character(as.list(match.call())$dt),dt,envir = .GlobalEnv)
    cat("outliers successfully removed","\n")
    return(invisible(dt))
  } else{
    cat("nothing changed","\n")
    return(invisible(var_name))
  }
}
P_Gad1 <- outlierKD(P_Gad1,Gad1)
#绘制markergene热图
sce.markers <- FindAllMarkers(object = sce, only.pos = TRUE, min.pct = 0.25, 
                              thresh.use = 0.25)
DT::datatable(sce.markers)
write.csv(sce.markers,file=paste0(pro,'_sce.markers.csv'))
library(dplyr) 
top10 <- sce.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(sce,top10$gene,size=3)
ggsave(filename=paste0(pro,'_sce.markers_heatmap.pdf'))
save(sce,sce.markers,file = 'first_sce.Rdata')
#查看关键基因在每群细胞中的表达情况（此处为神经元表达基因）
genes_to_check = c("Slc17a7","Satb2","Gad1","Gad2","Slc32a1","Slc6a1")
# All on Dotplot 
p <- DotPlot(sce, features = genes_to_check) + coord_flip()
p
# Annotate Inhineuron vs Non-Inhineuron clusters
# At this point we dont care for a more detailed annotation as we will annotate immune and non-immune separately later
dat=p$data 
Inhi=dat[dat$features.plot=='Gad1',]
fivenum(Inhi$avg.exp.scaled)
inhibitory=Inhi[Inhi$avg.exp.scaled > -0.5,]$id
inhibitory
sce@meta.data$inhi_annotation <-ifelse(sce@meta.data$RNA_snn_res.0.8 %in% inhibitory ,'Inhineuron','non-Inhineuron')
table(sce@meta.data$inhi_annotation)
p <- UMAPPlot(object = sce, group.by = 'inhi_annotation')
p
#将Inhibitoryneuron继续分群
table(sce@meta.data$inhi_annotation)
phe=sce@meta.data
table(phe$inhi_annotation)
cells.use <- row.names(sce@meta.data)[which(phe$inhi_annotation=='Inhineuron')]
length(cells.use)
sce_1 <-subset(sce, cells=cells.use) 
sce_1 <- NormalizeData(sce_1, normalization.method = "LogNormalize", 
                       scale.factor = 10000)
GetAssay(sce_1,assay = "RNA")
sce_1 <- FindVariableFeatures(sce_1, 
                              selection.method = "vst", nfeatures = 2000) 
sce_1 <- ScaleData(sce_1) 
sce_1 <- RunPCA(object = sce_1, pc.genes = VariableFeatures(sce_1)) 
sce_1 <- FindClusters(object = sce_1, verbose = T, resolution = 0.3)
set.seed(123)
sce_1 <- RunUMAP(object = sce_1, dims = 1:15, do.fast = TRUE)
DimPlot(sce_1,reduction = "umap",label=T)
#DimPlot(sce,reduction = "tsne",label=T, group.by = "patient_id")
table(sce_1@meta.data$seurat_clusters)
#singleR自动注释
sce_for_SingleR <- GetAssayData(sce_1, slot="data")
sce_for_SingleR
library(SingleR)
mpca.se <- MouseRNAseqData()
mpca.se
clusters=sce_1@meta.data$seurat_clusters
pred.mesc <- SingleR(test = sce_for_SingleR, ref = mpca.se, labels = mpca.se$label.main,
                     method = "cluster", clusters = clusters, 
                     assay.type.test = "logcounts", assay.type.ref = "logcounts")
table(pred.mesc$labels)
celltype = data.frame(ClusterID=rownames(pred.mesc), celltype=pred.mesc$labels, stringsAsFactors = F) 
sce_1@meta.data$singleR=celltype[match(clusters,celltype$ClusterID),'celltype']
DimPlot(sce_1, reduction = "umap", group.by = "singleR")
phe=sce@meta.data
table(phe$singleR)
save(phe,file = 'phe-of-subtypes-Immune-by-singleR.Rdata')
FeaturePlot(object = sce, 
            features ='Slc17a7', 
            cols= c("grey", "blue"), 
            reduction = "tsne")
#GSEA分析
cluster1.markers <- FindMarkers(sce_inhi, ident.1 = 1, min.pct = 0.01,logfc.threshold = 0.01)
geneList <- cluster1.markers$avg_logFC
names(geneList)=toupper(rownames(cluster1.markers))
geneList = sort(geneList,decreasing = T)
gmtfile='/home/wzh/Desktop/c5.go.v7.3.symbols.gmt'
geneset <- read.gmt( gmtfile ) 
egmt <- GSEA(geneList,TERM2GENE = geneset,minGSSize = 1,pvalueCutoff = 0.99,verbose = F)
gsea_results_df <- egmt@result
View(gsea_results_df)
gseaplot2(egmt,geneSetID ="GOBP_AXO_DENDRITIC_TRANSPORT", pvalue_table=T)
coregene <- gsea_results_df[gsea_results_df$ID=="GOBP_AXO_DENDRITIC_TRANSPORT",]$core_enrichment 
FeaturePlot(object = sce_inhi, 
            features ='Cnih2', 
            cols= c("grey", "red"), 
            reduction = "tsne")
# ＧＯ分析
# 寻找差异表达基因(所有分群)
# 该方法的优点是可以知道差异表达基因来自于哪个分群
dat=as.data.frame(sce_inhi@assays$RNA@data)
count=as.data.frame(sce_inhi@assays$RNA@counts)
clus=sce_inhi@active.ident
#也可以对特定两群进行比较
#sce1_for_DEG = subset(sce1,idents=c(1,2))
#dat=as.data.frame(sce1_for_DEG@assays$RNA@data)
#count=as.data.frame(sce1_for_DEG@assays$RNA@counts)
#clus=sce1_for_DEG@active.ident
prepare_for_DE <- function(count_matrix=count_matrix, clustering=clustering){
  expr_matrix <- as.matrix(count_matrix)
  sample_sheet <- data.frame(cells=colnames(count_matrix),  
                             cellType=clustering)
  rownames(sample_sheet)<- colnames(count_matrix)
  gene_annotation <- as.data.frame(rownames(count_matrix))
  rownames(gene_annotation)<- rownames(count_matrix)
  colnames(gene_annotation)<- "genes"
  pd <- new("AnnotatedDataFrame", data = sample_sheet)
  fd <- new("AnnotatedDataFrame", data = gene_annotation)
  HSMM <- newCellDataSet(
    as(expr_matrix,"sparseMatrix"),
    phenoData = pd,
    featureData = fd,
    lowerDetectionLimit=0.5,
    expressionFamily=negbinomial.size()
  )
  HSMM <- detectGenes(HSMM, min_expr = 5)
  HSMM <- HSMM[fData(HSMM)$num_cells_expressed > 10, ]
  HSMM <- estimateSizeFactors(HSMM)
  HSMM <- estimateDispersions(HSMM)
  return(HSMM)
}

findDEgenes <- function(HSMM=HSMM, qvalue=qvalue){
  diff_test_res <- differentialGeneTest(
    HSMM,
    fullModelFormulaStr="~cellType"
  )
  
  sig_genes_0.05 <- subset(diff_test_res, qval < 0.05)
  sig_genes_0.01 <- subset(diff_test_res, qval < 0.01)
  
  print(paste(nrow(sig_genes_0.05), " significantly DE genes (FDR<0.05).", sep=""))
  print(paste(nrow(sig_genes_0.01), " significantly DE genes (FDR<0.01).", sep=""))
  
  diff_test_res <- subset(diff_test_res, qval< qvalue)
  
  return(diff_test_res)
}

get_up_reg_clusters <- function(count, clustering, DE_genes){
  cluster_nb <- unique(clustering)
  mean_per_cluster <- vector()
  DE_genes <- DE_genes[order(rownames(DE_genes)),]
  count <- count[order(rownames(count)),]
  count_de_genes <- count[rownames(count) %in% DE_genes$genes,]
  print(dim(count_de_genes))
  for (clusters in cluster_nb) {
    # print(head(count_de_genes[,
    # 		colnames(count_de_genes) %in% names(clustering[clustering==clusters])
    # 	]))
    mean <- rowMeans(
      as.matrix(count_de_genes[,
                               colnames(count_de_genes) %in% names(clustering[clustering==clusters])
                               ])
    )
    names(mean) <- clusters
    mean_per_cluster <- cbind(
      mean_per_cluster,
      mean
    )
  }
  colnames(mean_per_cluster) <- cluster_nb
  up_reg_cluster <- colnames(mean_per_cluster)[apply(mean_per_cluster,1,which.max)]
  de_genes_table <- data.frame(
    DE_genes,
    mean_per_cluster,
    cluster=up_reg_cluster
  )
  
  return(de_genes_table)
}

DE_scn <- prepare_for_DE (
  count, 
  clus
)

scn_DE_genes <- findDEgenes(
  DE_scn, 
  qvalue=0.05
)

de_clusters <- get_up_reg_clusters(
  dat, 
  clus, 
  scn_DE_genes
)
#常规气泡图
de_genes <- de_clusters
gene_names <- subset(de_genes, qval<0.05)
gene_names <- gene_names$genes
entrez_genes <- bitr(gene_names, fromType="SYMBOL", 
                     toType="ENTREZID", 
                     OrgDb="org.Mm.eg.db")
entrez_genes <- entrez_genes[!entrez_genes$ENTREZID %in% "101055843",]

de_gene_clusters <- de_genes[de_genes$genes %in% entrez_genes$SYMBOL,
                             c("genes", "cluster")]
de_gene_clusters <- data.frame(
  ENTREZID=entrez_genes$ENTREZID[entrez_genes$SYMBOL %in% de_gene_clusters$genes],
  cluster=de_gene_clusters$cluster
)
table(de_gene_clusters$cluster)
list_de_gene_clusters <- split(de_gene_clusters$ENTREZID, 
                               de_gene_clusters$cluster)
formula_res <- compareCluster(
  ENTREZID~cluster, 
  data=de_gene_clusters, 
  fun="enrichGO", 
  OrgDb="org.Mm.eg.db",
  ont		   = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.01,
  qvalueCutoff  = 0.05
)
dotplot(formula_res, showCategory=5)
#查看第n群中的差异基因并进行GO分析
DE_1 <- de_clusters
#DE_1 <- de_clusters[de_clusters$cluster==c(1,2,3,4,5),]
DE_1_gene <- DE_1$genes
DE_1_gene <- as.character(DE_1_gene)
entrez_genes <- bitr(DE_1_gene, fromType="SYMBOL", 
                     toType="ENTREZID", 
                     OrgDb="org.Mm.eg.db")
GO_1 <- enrichGO(gene=entrez_genes$ENTREZID,
                 OrgDb = 'org.Mm.eg.db',
                 ont='BP',
                 pAdjustMethod = 'BH',
                 pvalueCutoff = 0.01,
                 qvalueCutoff = 0.05)
BH<-function(ratio){
  sapply(ratio,function(x) as.numeric(strsplit(x,"/")[[1]][1])/as.numeric(strsplit(x,"/")[[1]][2]))
}
GO_1@result$FOLD=BH(GO_1@result$GeneRatio)/BH(GO_1@result$BgRatio)
dim(GO_1)
GO_1_2 <- simplify(GO_1,cutoff=0.7,by="p.adjust",select_fun=min)
dim(GO_1_2)
bar <- GO_1[,c(2,10)]
bar <- bar[order(-bar[,2]),][c(1,9,13,18,19),]
p <- ggplot(data=bar, aes(x=fct_reorder(Description,FOLD), y=FOLD))+
  geom_bar(stat='identity')+
  labs(x=' ',y='Fold change')+
  coord_flip()
### 寻找差异表达基因(所有分群)方法２
#sce1_for_DEG = subset(sce1,idents=c(1,2))
count_matrix = sce1@assays$RNA@counts
cluster=sce1@active.ident
table(cluster)
expr_matrix <- as.matrix(count_matrix)
sample_sheet <- data.frame(cells=colnames(count_matrix),  
                           cellType=cluster)
rownames(sample_sheet)<- colnames(count_matrix)
gene_annotation <- as.data.frame(rownames(count_matrix))
rownames(gene_annotation)<- rownames(count_matrix)
colnames(gene_annotation)<- "genes"
pd <- new("AnnotatedDataFrame", data = sample_sheet)
fd <- new("AnnotatedDataFrame", data = gene_annotation)
HSMM <- newCellDataSet(
  as(expr_matrix,"sparseMatrix"),
  phenoData = pd,
  featureData = fd,
  lowerDetectionLimit=0.5,
  expressionFamily=negbinomial.size())
HSMM <- detectGenes(HSMM, min_expr = 5)
HSMM <- HSMM[fData(HSMM)$num_cells_expressed > 10, ]
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)

cds=HSMM
# 单细胞转录组最重要的就是把细胞分群啦，这里可供选择的算法非常多，我们首先演示PCA结果。
# 并不是所有的基因都有作用，所以先进行挑选，合适的基因用来进行聚类。
#disp_table <- dispersionTable(cds)
#unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
#cds <- setOrderingFilter(cds, unsup_clustering_genes$gene_id)
#plot_ordering_genes(cds) 
# plot_pc_variance_explained(cds, return_all = F) # norm_method='log'
# 其中 num_dim 参数选择基于上面的PCA图
#cds <- reduceDimension(cds, max_components = 2, num_dim = 6,
#reduction_method = 'tSNE', verbose = T)
#cds <- clusterCells(cds, num_clusters = 5) 
#plot_cell_clusters(cds, color = "cellType")
#table(pData(cds)$Cluster,cluster)
#plot_cell_clusters(cds, 1, 2 )
#计算差异表达基因
diff_test_res <- differentialGeneTest(cds,
                                      fullModelFormulaStr = "~cellType")
# Select genes that are significant at an FDR < 10%
sig_genes <- subset(diff_test_res, qval < 0.1)
dim(sig_genes)
# 热图可视化
count_matrix[1:4,1:4]
dim(count_matrix)
table(cluster)
htmapGenes='Gpr17,Meg3,Olig1,Pdgfra,Ptprz1,Snhg11'
htmapGenes=strsplit(htmapGenes,',')[[1]]
htmapGenes %in% rownames(sig_genes)
library(pheatmap)
dat=count_matrix[htmapGenes,]
pheatmap(dat)
n=t(scale(t(dat)))
n[n>2]=2 #限定上限，使表达量大于2的等于2
n[n< -2]= -2 #限定下限，使表达量小于-2的等于-2
n[1:4,1:4]
pheatmap(n,show_colnames =F,show_rownames = F)
ac=data.frame(group=cluster)
#ac = ac[order(ac$group),,drop=F]
rownames(ac)=colnames(n)
pheatmap(n,annotation_col = ac,
         show_colnames =F,show_rownames = T)
n[n< -1]= -1 #限定下限，使表达量小于-2的等于-2
n[1:4,1:4] 
pheatmap(n,annotation_col = ac,
         show_colnames =F,show_rownames = T,cluster_cols = F,cluster_rows = F)
DE_1_gene <- sig_genes$genes
DE_1_gene <- as.character(DE_1_gene)
entrez_genes <- bitr(DE_1_gene, fromType="SYMBOL", 
                     toType="ENTREZID", 
                     OrgDb="org.Mm.eg.db")
GO_1 <- enrichGO(gene=entrez_genes$ENTREZID,
                 OrgDb = 'org.Mm.eg.db',
                 ont='BP',
                 pAdjustMethod = 'BH',
                 pvalueCutoff = 0.01,
                 qvalueCutoff = 0.05)



