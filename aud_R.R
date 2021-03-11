library(Seurat)
library(pheatmap)
library(dplyr)
library(monocle)
library(data.table)
library(ROTS)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
library(GOplot)

###P10降维聚类
aud1 <- read.csv('/home/wzh/Desktop/aud/GSE140883_P10_1.tsv',header=T, row.names=1, sep="\t") 
aud_1_1 <- read.csv('/home/wzh/Desktop/aud/GSE140883_P10_2.tsv',header=T, row.names=1, sep="\t")
aud1 <- t(aud1)
aud_1_1 <- t(aud_1_1)
aud_1 <- cbind(aud1,aud_1_1)
#按照列，对每一个细胞进行内部归一化，主要是统一文库大小。
#dataaud_1 <- log2(1 + sweep(aud_1, 2, median(colSums(aud_1))/colSums(aud_1), '*'))
#dim(dataaud_1)
#统计基因在多少个细胞里面表达
fivenum(apply(aud_1,1,function(x) sum(x>0)))
#统计细胞表达多少基因
fivenum(apply(aud_1,2,function(x) sum(x>0)))
aud_P10 <- CreateSeuratObject(counts = aud_1,min.cell=1,min.features=200)  
sce1 <- NormalizeData(aud_P10)
sce1 <- ScaleData(sce1, display.progress  = F)
sce1 <- FindVariableFeatures(object = sce1, mean.function = ExpMean, dispersion.function = LogVMR, mean.cutoff = c(0.0125,3), dispersion.cutoff = c(0.5,Inf)) 
sce1 <- RunPCA(object = sce1, pc.genes = VariableFeatures(sce1))
sce1 <- FindNeighbors(sce1, dims = 1:10)
sce1 <- FindClusters(sce1, resolution = 0.5)
sce1 <- RunUMAP(object = sce1, dims=1:10)
#存储Ｒ数据
save(sce1,file='aud_p10.Rdata')
DimPlot(sce1,reduction = "umap",label = T,label.size=5,pt.size=0.6)

### 寻找差异表达基因(所有分群)
#该方法的优点是可以知道差异表达基因来自于哪个分群
dat=as.data.frame(sce1@assays$RNA@data)
count=as.data.frame(sce1@assays$RNA@counts)
clus=sce1@active.ident
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
#查看第１群中的差异基因并进行富集分析
DE_1 <- de_clusters[de_clusters$cluster==1,]
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

###marker基因表达热图
sce.markers <- FindAllMarkers(sce1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top5 <- sce.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(sce1, features = top5$gene) + NoLegend()

###Feature plot
cluster0.markers <- FindMarkers(sce1, ident.1 = 0, min.pct = 0.25)
markers_genes =  rownames(head(x = cluster0.markers, n = 5))
FeaturePlot(object = sce1, 
            features ='Meg3', 
            cols= c("grey", "blue"), 
            reduction = "umap")

