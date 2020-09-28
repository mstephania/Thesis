#####################################################################################
## Script used for PhD thesis titled:
## "Cancer Immunology of Cutaneous Melanoma A Systems Biology Approach"
## Author:
## Mindy Muñoz
## Bioinformatic Program, University of São Paulo
#####################################################################################
#####################################################################################
## Feel free to copy any section of the script, just remember to reference the work
#####################################################################################
library(Seurat)
library(ggplot2)
library(cowplot)
library(patchwork)
library(DESeq2)
melanoma_seurat <- read.table("melanoma_seurat.txt",sep=",",header=TRUE,row.names=1)
dim(melanoma_seurat)
head (melanoma_seurat)

metadata <- read.csv(file = "metadatamelanoma", header = TRUE, sep = " ", quote = "\"")
melanoma <- CreateSeuratObject(melanoma_seurat,assay = "RNA", meta.data = metadata)

melanoma.norm<- NormalizeData(melanoma, verbose = FALSE)
melanoma_var <- FindVariableFeatures(melanoma.norm, selection.method = "vst", nfeatures = 2000)
melanoma_scale <- ScaleData(melanoma_var, verbose = FALSE)

melanoma_pca <- RunPCA(melanoma_scale, npcs = 30, verbose = FALSE)
tsne.integrated <- RunTSNE(object = melanoma_pca, approx=FALSE, dims.use = 1:10)
melanoma_umap <- RunUMAP(melanoma_pca, reduction = "pca", dims = 1:30)
melanoma_vc <- FindNeighbors(melanoma_pca, dims = 1:20, verbose = FALSE)
melanoma_cluster <- FindClusters(melanoma_vc, verbose = FALSE)

DimPlot(melanoma_cluster, label = TRUE) + NoLegend()
FeaturePlot(melanoma_cluster, features = c("CD14", "CD4", "CD8A","DNMT1", "S100A9", "CD79A","LYZ","CD3E","GNLY"))
DimPlot(melanoma_cluster, group.by = c("tumor.Stage", "cell.subtype"))
DimPlot(tsne.integrated, group.by = c("tumor.Stage", "cell.subtype"))

FeaturePlot(tsne.integrated, features = c("PTGS2", "THBS1", "TNFAIP3","CCL3", "FKBP5", "ANPEP","PTX3","CCL4","C3AR1"))

melanoma_var <- FindVariableFeatures(melanoma.norm, selection.method = "mvp")
melanoma_scale <- ScaleData(melanoma_var, verbose = FALSE)
melanoma_pca <- RunPCA(melanoma_scale, npcs = 30, verbose = FALSE)
tsne.integrated <- RunTSNE(object = melanoma_pca, approx=FALSE, dims.use = 1:10)
melanoma_vc <- FindNeighbors(tsne.integrated, dims = 1:20, verbose = FALSE)
melanoma_cluster <- FindClusters(melanoma_vc, verbose = FALSE)

markers.to.plot <- c("CD3D", "HSPH1", "SELL", "CACYBP", "GNLY", "NKG7", "CCL5",
                     "CD8A", "MS4A1", "CD79A", "NME1", "FCGR3A", "S100A9", "HLA-DQA1",
                     "GPR183", "PPBP", "GNG11", "HBA2", "HBB", "TSPAN13", "IL3RA")

FeaturePlot(melanoma_cluster, features = c("S100A8", "PI3","HAL","PDK4"))

VlnPlot(melanoma_cluster, c("S100A8", "PI3","HAL","PDK4","S100A9","S100A4","S100A10"), group.by = "tumor.Stage")

FeaturePlot(tsne.integrated, features = c("CD8A", "CD4","CD14", "GNLY", "CD3E","LYZ", "CD40","MITF","FCGR3A"))
              
new.cluster.ids <- c("CD8+ T", "CD14+ Monocyte","CD4+ T")
names(new.cluster.ids) <- levels(melanoma_cluster)
melanoma_new_names <- RenameIdents(melanoma_cluster, new.cluster.ids)

DotPlot(melanoma_new_names, features = rev(markers.to.plot), cols = c("blue", "red"), dot.scale = 8)+ RotatedAxis()
markers.to.plot <- c("CXCL1","CD14","CSF1","CD4","CCL5","CD8A","CXCL5","CXCL8","RAB32")
markersCD8.to.plot <- c("RPS18P12","RPS28P7","TRBV13","JUND","JUNB","CX3CR1","ITGB1","AC022149.1","TRBV2","CCR4")
markersCD14.to.plot <- c("RNA5-8SN2.4", "FP671120.4","FP236383.2","RN7SK","MALAT1","PTGS2","THBS1","TNFAIP3","CCL3","FKBP5") 
markersCD4.to.plot <- c("HBB","LYZ","PPBP","HBA2","MT-TA", "AC015911.1","AC090602.1","AL590135.1","AC016700.2","AL034379.1")
all.markers.to.plot <- c("RNA5-8SN2.4", "FP671120.4","FP236383.2","RN7SK","MALAT1","PTGS2","THBS1","TNFAIP3","CCL3","FKBP5","HBB","LYZ","PPBP","HBA2","MT-TA", "AC015911.1","AC090602.1","AL590135.1","AC016700.2","AL034379.1","RPS18P12","RPS28P7","TRBV13","JUND","JUNB","CX3CR1","ITGB1","AC022149.1","TRBV2","CCR4")
DotPlot(melanoma_new_names, features = rev(markersCD4.to.plot), cols = c("blue", "red"), dot.scale = 8)+ RotatedAxis()
DotPlot(melanoma_new_names, features = rev(markersCD8.to.plot), cols = c("blue", "red"), dot.scale = 8)+ RotatedAxis()
DotPlot(melanoma_new_names, features = rev(all.markers.to.plot), cols = c("blue", "red"), dot.scale = 8)+ RotatedAxis()


RColorBrewer::brewer.pal.info
levels(melanoma_cluster)
melanoma_cluster
melanoma_cluster$cell.subtype
levels(melanoma_cluster$cell.subtype)
melanoma_cluster$tumor.Stage
  
melanoma_ft <- FindVariableFeatures(melanoma_cluster, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(melanoma_ft), 10)
top10
plot1 <- VariableFeaturePlot(melanoma_ft)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = FALSE)
CombinePlots(plots = list(plot1, plot2))

m.markers <- FindMarkers(melanoma_cluster, ident.1 = 0)
head(m.markers)
write.csv(m.markers, "melanoma_markers_CD8.csv", row.names=TRUE)
m1.markers <- FindMarkers(melanoma_cluster, ident.1 = 1)
head(m1.markers)
write.csv(m1.markers, "melanoma_markers_CD14.csv", row.names=TRUE)
m2.markers <- FindMarkers(melanoma_cluster, ident.1 = 2)
head(m2.markers)
write.csv(m2.markers, "melanoma_markers_CD4.csv", row.names=TRUE)

CD14_stage.markers <- FindMarkers(melanoma_cluster, ident.1 = "Stage-IV-melanoma", group.by = "tumor.Stage", subset.ident = "1")

CD14_subset <- subset(melanoma_cluster, idents = '1')
write.csv(CD14_subset, "CD14_subset.csv", row.names=TRUE)
CD8_subset <- subset(melanoma_cluster, idents = '0')
write.csv(CD8_subset, "CD8_subset.csv", row.names=TRUE)
CD4_subset <- subset(melanoma_cluster, idents = '2')
write.csv(CD4_subset, "CD4_subset.csv", row.names=TRUE)
VlnPlot(CD14_subset, c("RNA5-8SN2.1","RNA5-8SN2.2","RNA5-8SN2.3","RNA5-8SN2.4","FP671120.4","FP236383.2","FP236383.3","FP671120.3", "HBB","RN7SK","MIR6087", "JUN"), group.by = "tumor.Stage")
VlnPlot(CD4_subset, c("PTEN","MITF","NRAS","BRAF", "NF1","CDK4"),group.by = "tumor.Stage")
VlnPlot(CD14_subset, c("PTEN"),group.by = "tumor.Stage")
VlnPlot(CD8_subset, c("PTEN"),group.by = "tumor.Stage")
VlnPlot(CD8_subset, c("TNFSF14"),group.by = "tumor.Stage")
VlnPlot(CD14_subset, c("CD74"),group.by = "tumor.Stage")

m_stage.markers <- ?FindMarkers(melanoma_cluster, ident.1 = "Stage-IV-melanoma", group.by = "tumor.Stage", subset.ident = "0")
write.csv(m_stage.markers, "melanoma_markers_CD8_stage.csv", row.names=TRUE)
subset <- m_stage.markers[order(m_stage.markers$avg_logFC),]
subset2 <- subset[ subset$p_val< 0.05 ,]
subset3 <- subset2[order(subset2$p_val),] 
subset4 <- subset3[order(subset3$avg_logFC),] 
write.csv(subset4, "melanoma_markers_CD8_stage_filtered.csv", row.names=TRUE)
m1_stage.markers <- FindMarkers(melanoma_cluster, ident.1 = "Stage-IV-melanoma", group.by = "tumor.Stage", subset.ident = "1")
write.csv(m1_stage.markers, "melanoma_markers_CD14_stage.csv", row.names=TRUE)
subset <- m1_stage.markers[order(m1_stage.markers$avg_logFC),] 
subset2 <- subset[ subset$p_val< 0.05 ,]
subset3 <- subset2[order(subset2$p_val),] 
subset4 <- subset3[order(subset3$avg_logFC),] 
write.csv(subset4, "melanoma_markers_CD14_stage_filtered.csv", row.names=TRUE)
m2_stage.markers <- FindMarkers(melanoma_cluster, ident.1 = "Stage-IV-melanoma",ident.2 = "Health", group.by = "tumor.Stage", subset.ident = "2")
write.csv(m2_stage.markers, "melanoma_markers_CD4_stage.csv", row.names=TRUE)
subset <- m2_stage.markers[order(m2_stage.markers$avg_logFC),] 
subset2 <- subset[ subset$p_val < 0.05 ,]
subset3 <- subset2[order(subset2$p_val),] 
subset4 <- subset3[order(subset3$avg_logFC),] 
write.csv(subset4, "melanoma_markers_CD4_stage_filtered.csv", row.names=TRUE)

DoHeatmap(object = melanoma_cluster)
DoHeatmap(object = melanoma_cluster, group.by = c("cell.subtype"), group.bar = TRUE)
DoHeatmap(object = melanoma_cluster, group.by = c("tumor.Stage"), group.bar = TRUE)
DoHeatmap(object = melanoma_cluster, group.by = c("tumor.Stage","cell.subtype"), label = FALSE)

DoHeatmap(object = melanoma_cluster, group.by = c("tumor.Stage"), label=FALSE, group.bar = TRUE)
DoHeatmap(object = melanoma_cluster, group.by = c("tumor.Stage","cell.subtype"), label=FALSE, group.colors = c("red", "blue"))


VlnPlot(CD8_subset, c("ABCD1","CCR4", "TNFSF14","LAG3", "ABCD3"),group.by = "tumor.Stage")
VlnPlot(CD8_subset, c("CXCL12","HSDF1","TPAR1","SDF1A","SDF1B","SDF1","TLSF","HIRH","CCR4"),group.by = "tumor.Stage")

head(m_stage.markers)
head(m1_stage.markers)
head(m2_stage.markers)
FeaturePlot(melanoma_cluster, features = c("RPL35P2", "RPS28P7", "ARL6IP5","THBS1", "FKBP5", "TNFAIP3","AC010240.2","RPL10P18","AC091839.1"))
melanoma.markers <- FindAllMarkers(melanoma_cluster, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

mh_stage.markers <- FindMarkers(melanoma_cluster, ident.1 = "Health", group.by = "tumor.Stage", subset.ident = "0")
write.csv(mh_stage.markers, "markers_sobreh_CD8_stage.csv", row.names=TRUE)
m1h_stage.markers <- FindMarkers(melanoma_cluster, ident.1 = "Health", group.by = "tumor.Stage", subset.ident = "1")
write.csv(m1h_stage.markers, "markers_sobreh_CD14_stage.csv", row.names=TRUE)
m2h_stage.markers <- FindMarkers(melanoma_cluster, ident.1 = "Health", group.by = "tumor.Stage", subset.ident = "2")
write.csv(m2h_stage.markers, "markers_sobreh_CD4_stage.csv", row.names=TRUE)

head(mh_stage.markers)
head(m1h_stage.markers)
head(m2h_stage.markers)

head(mh_stage.markers[order(mh_stage.markers$avg_logFC),])
head(m_stage.markers[order(m_stage.markers$avg_logFC),])
head(m1h_stage.markers[order(m1h_stage.markers$avg_logFC),])
head(m1_stage.markers[order(m1_stage.markers$avg_logFC),])
head(m2h_stage.markers[order(m2h_stage.markers$avg_logFC),])

subset <- m2_stage.markers[order(m2_stage.markers$avg_logFC),] 
subset2 <- subset[ subset$p_val_adj < 0.01 ,]
subset3 <- subset2[order(subset2$p_val),] 
subset4 <- subset3[order(subset3$avg_logFC),] 

head(subset4)
FeaturePlot(melanoma_cluster, features = c("AL034379.1", "AC090602.1", "RPL38P4","THBS1", "FKBP5", "CXCL8","RPS28P7","AC012618.1","AL450998.1"))
subset <- m2h_stage.markers[order(m2h_stage.markers$avg_logFC),] 
subset2 <- subset[ subset$p_val_adj< 0.01 ,]
subset3 <- subset2[order(subset2$p_val),] 
subset4 <- subset3[order(subset3$avg_logFC),] 
head(subset4)

FeaturePlot(melanoma_cluster, features = c("CCR4", "ANXA2P2", "EIF4EBP2","RNA5-8SN2.4", "FKBP5", "THBS1","HBB","CCR4","CXCR4"))
FeaturePlot(melanoma_cluster, features = c("MALAT1", "PTGS2", "RN7SK","RNA5-8SN2.4", "MIR6087", "THBS1","HBB","CCR4","CXCR4"))


FeaturePlot(melanoma_cluster, features = c("CD14", "CD44"))
VlnPlot(melanoma_cluster, c("CD8A", "CD14", "CCL5", "S100A4", "ANXA1", "CCR7", "CD74", "CD3D","ISG15"), group.by = "cell.subtype")

melanoma.markers
write.csv(melanoma.markers, "melanoma_markers_seurat.csv", row.names=TRUE)

subset <- m2h_stage.markers[order(m2h_stage.markers$avg_logFC),] 
subset2 <- subset[ subset$p_val_adj< 0.05 ,]
subset3 <- subset2[order(subset2$p_val),] 
subset4 <- subset3[order(subset3$avg_logFC),] 
head(subset4)
write.csv(subset4, "markers_sobreh_filtropvaj_menorq1_CD4_stage.csv", row.names=TRUE)
subset <- m1h_stage.markers[order(m1h_stage.markers$avg_logFC),] 
subset2 <- subset[ subset$p_val< 1 ,]
subset3 <- subset[order(subset$p_val),] 
subset4 <- subset3[order(subset3$avg_logFC),] 
head(subset4)
write.csv(subset4, "markers_sobreh_filtropv_menorq1_CD14_stage.csv", row.names=TRUE)
subset <- mh_stage.markers[order(mh_stage.markers$avg_logFC),] 
subset2 <- subset[ subset$p_val_adj< 0.05 ,]
subset3 <- subset2[order(subset2$p_val),] 
subset4 <- subset3[order(subset3$avg_logFC),] 
head(subset4)
write.csv(subset4, "markers_sobreh_filtropvaj_menorq1_CD8_stage.csv", row.names=TRUE)

FeaturePlot(melanoma_cluster, features = c("CCR4","RPL38P4", "SESN3", "B2M", "RECQL", "UHMK1","C3AR1","CCL4","ANPEP"))
FeaturePlot(melanoma_cluster, features = c("EZR", "ITGA2", "ITGA3", "ITGB4", "SMAD3","TGFBR2", "RECQL"))

metadata <- read.csv(file = "metadatamelanoma", header = TRUE, sep = " ", quote = "\"")
melanoma2 <- CreateSeuratObject(melanoma_seurat, meta.data = metadata)
melanoma.list <- SplitObject(melanoma2, split.by = "tumor.Stage")
head(melanoma.list)
for (i in 1:length(melanoma.list)) {
  melanoma.list[[i]] <- NormalizeData(melanoma.list[[i]], verbose = FALSE)
  melanoma.list[[i]] <- FindVariableFeatures(melanoma.list[[i]], selection.method = "vst",
                                             nfeatures = 2000, verbose = FALSE)
}

