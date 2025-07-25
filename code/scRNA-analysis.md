# Loading R packages
```
#library(packages)
library(Seurat)
library(harmony)
library(ggsci)
library(BayesSpace)
library(stringr)
library(plyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(cowplot)
library(scales)
library(dittoSeq)
library(pheatmap)
library(ggpubr)
library(dplyr)
library(RColorBrewer)
library(scales)
library(viridis)
library(SingleCellExperiment)
library(reshape2)
library(EnhancedVolcano)
library(msigdbr)
library(GSVA)
library(rstatix)
library(reshape2)
library(SpaGene)
library(ggpointdensity)
library(cowplot)
library(tidyr)
library(rjags)

#Set the colors to be used
cols <- colorRampPalette(brewer.pal(10, "RdBu"))(50) %>% rev()
col1 <- c(pal_npg()(10),pal_igv()(9),pal_uchicago("light")(9),pal_futurama()(12), pal_lancet()(9))[-8]#UMAP
col2 <- viridis(50)#marker
col3 <- inferno(50)#heatmap
col4 <- rocket(50)#dotplot
```
# Integrated the scRNA dataset
### integrate single cell RNA data into a object
```
setwd("Path to your data/JackyWu_PDAC/analysis")
load("CRA001160_seurat.Rdata")
load("GSE154778_seurat.Rdata")
load("GSE155698_seurat.Rdata")
load("GSE197177_seurat.Rdata")
load("GSE212966_seurat.Rdata")
CRA001160$orig.ident <- "CRA001160"
CRA001160@active.ident <- as.factor(CRA001160$orig.ident)
PDAC_scRNA<-merge(CRA001160,c(GSE154778,GSE155698,GSE197177,GSE212966))
save(PDAC_scRNA,file = "PDAC_scRNA.Rdata")
#
remove(CRA001160)
remove(GSE154778)
remove(GSE155698)
remove(GSE197177)
remove(GSE212966)
gc()

```
# scRNA-seq filtering
### scRNA-data filtering
```
load('PDAC_scRNA.Rdata')
#calculate mitochondrial, hemoglobin and ribosomal gene counts
PDAC_scRNA <- PercentageFeatureSet(PDAC_scRNA, pattern = "^MT-", col.name = "pMT")
p<-VlnPlot(PDAC_scRNA, features = c("nFeature_RNA", "nCount_RNA", "pMT"), cols = col1,pt.size = 0, group.by = "orig.ident", ncol = 1, log = T)

#filtered scRNA object
nFeature_lower <- 200
nFeature_upper <- 10000
nCount_lower <- 1000
nCount_upper <- 100000
pMT_lower <- 0
pMT_upper <- 10
PDAC_scRNA_sub <- subset(PDAC_scRNA, subset = nFeature_RNA > nFeature_lower & nFeature_RNA < nFeature_upper & nCount_RNA > nCount_lower & nCount_RNA < nCount_upper & pMT < pMT_upper )
```
### 2.Major cell type identification
```
marker_features<-c('PTPRC',#immune type
                   "CD3E",'CD3D',#T cell
                   'CD4','IL7R',#CD4 T cells
                   'CD8A','CD8B',#CD8T
                   'NKG7','GNLY',#NK
                   'CD79A','MS4A1',#B cells
                   'IGLC2','IGHA1',#B plasma cells
                   'CD1C',#Dendritic cell
                   'CD68','MARCO',#Macrophage
                   'CD14','S100A8',#monocyte
                   'CPA3','KIT',#Mast
                   'COL1A1','DCN',#Fibroblast
                   'PECAM1','VWF',#Endothelial
                   'ACTA2','CSPG4',#Pericyte
                   'EPCAM','KRT8'#Epithelial                  
)


DotPlot(PDAC_scRNA,features = marker_features,group.by = 'seurat_clusters',assay = "SCT")+scale_color_gradient2(low = 'white',
                                                                                                                mid = col4[50],
                                                                                                                high = col4[20])+coord_flip()
```


### scRNA CAF
```
load("analysis/PDAC_TME_sc.Rdata")
table(PDAC_TME_sc$celltype)
Fibro <- subset(PDAC_TME_sc, celltype == "Fibroblast")
Fibro <- RunPCA(Fibro)
Fibro <- RunHarmony(Fibro,group.by.vars = "sample")
#PDAC_scRNA <- JackStraw(PDAC_scRNA, num.replicate = 10,reduction = "pca",assay = "RNA")
ElbowPlot(Fibro,reduction = "harmony")

Fibro <- RunUMAP(Fibro,reduction = "harmony",dims = 1:30)
Fibro <- FindNeighbors(Fibro, reduction = "harmony", dims = 1:30)
Fibro <- FindClusters(Fibro,
                      resolution =0.3,
                      method ='igraph',
                      verbose = T)# 6
#defined celltype subtype
table(Fibro$SCT_snn_res.0.3)
cluster<- c("Fibro_c1",
            "Fibro_c2",
            "Fibro_c3",
            "Fibro_c4",
            "Fibro_c5",
            "Fibro_c6"
            )
new.cluster.ids <- cluster
names(new.cluster.ids)<- levels(Fibro)
Fibro <- RenameIdents(Fibro, new.cluster.ids)
Fibro@meta.data$celltype_subtype<-Fibro@active.ident

table(Fibro@active.ident)
Fibro$site <- factor(Fibro$site,levels = c("Normal pancreas",'Adjacent',"Liver metastases","PDAC primary"))
#UMAP plot
#Fig 3A
DimPlot(Fibro,group.by = "celltype_subtype",cols = col1[1:6],raster = T,split.by = "site")

#Dotplot
markers<-c("IL6",'PDGFRA','CXCL12','CFD','CXCL1',#iCAF
           "ACTA2",'MMP11','MYL9','HOPX','TPM1',#myCAF
           "HLA-B",'SAA3','CD74'#apCAF antigen presenting CAFs
           )
DotPlot(Fibro,features = markers,group.by = "celltype_subtype")+scale_color_gradient2(low = 'white',mid = col4[50],high = col4[20])+
  theme_bw()+theme(panel.grid.major = element_blank() ,axis.text.x = element_text(size = 15,color = "black",angle = 45,hjust = 1,vjust = 1),
                   axis.text.y = element_text(size = 15,color = "black"),
                   axis.title = element_text(size = 15,color = 'black'),
                   legend.text = element_text(size = 12,color = "black"))+labs(x="",y="")

```

### scRNA Macro/Mono
```
Macro <- subset(PDAC_TME_sc, celltype == "Macro/Mono")
Macro <- RunPCA(Macro)
Macro <- RunHarmony(Macro,group.by.vars = "sample")
#PDAC_scRNA <- JackStraw(PDAC_scRNA, num.replicate = 10,reduction = "pca",assay = "RNA")
ElbowPlot(Macro,reduction = "harmony")

Macro <- RunUMAP(Macro,reduction = "harmony",dims = 1:30)
Macro <- FindNeighbors(Macro, reduction = "harmony", dims = 1:30)
Macro <- FindClusters(Macro,
                      resolution =0.2,
                      method ='igraph',
                      verbose = T)# 6
#defined celltype subtype
table(Macro$SCT_snn_res.0.2)
cluster<- c("Macro_c1",
            "Macro_c2",
            "Macro_c3",
            "Macro_c4",
            "Macro_c5",
            "Macro_c6",
            "Macro_c7",
            "Macro_c8"
)
new.cluster.ids <- cluster
names(new.cluster.ids)<- levels(Macro)
Macro <- RenameIdents(Macro, new.cluster.ids)
Macro@meta.data$celltype_subtype<-Macro@active.ident

table(Macro@active.ident)
Macro$site <- factor(Macro$site,levels = c("Normal pancreas",'Adjacent',"Liver metastases","PDAC primary"))
#UMAP plot
DimPlot(Macro,group.by = "celltype_subtype",cols = col1[8:16],raster = T,split.by = "site")


#Dotplot
markers <- c("CD68","CD163",'CD14','FCN1')
 DotPlot(Macro,features = markers,group.by = "celltype_subtype")+scale_color_gradient2(low = 'white',mid = col4[50],high = col4[20])+
  theme_bw()+theme(panel.grid.major = element_blank() ,axis.text.x = element_text(size = 15,color = "black",angle = 45,hjust = 1,vjust = 1),
                   axis.text.y = element_text(size = 15,color = "black"),
                   axis.title = element_text(size = 15,color = 'black'),
                   legend.text = element_text(size = 12,color = "black"))+labs(x="",y="")

```

### scRNA analysis of T/NK subset
```
load('analysis/PDAC_scRNA_Runharmony_for_using.Rdata')

Tcell <- subset(PDAC_scRNA_sub,idents="T/NK")
Tcell<-subset(Tcell, subset = PTPRC > 0,slot="counts")

load("analysis/TME/Tcell.Rdata")
DotPlot(Tcell,features = c('EPCAM','KRT19','KRT18','KRT8','KRT7','CD79A','CD68','CD14','PECAM1','COL1A1','FN1'),group.by = "sample")

Tcell<-subset(Tcell, subset = EPCAM == 0,slot="counts")
Tcell<-subset(Tcell, subset = KRT81 == 0,slot="counts")
Tcell<-subset(Tcell, subset = KRT86 == 0,slot="counts")
Tcell<-subset(Tcell, subset = CD79A == 0,slot="counts")


Tcell <- SCTransform(Tcell,assay = "SCT")
Tcell <- FindVariableFeatures(Tcell,assay = "SCT")

varfeature <- Tcell@assays[["SCT"]]@var.features
varfeature <- varfeature[-grep('^TRA',varfeature)]
varfeature <- varfeature[-grep('^TRB',varfeature)]

Tcell@assays$SCT@var.features<-varfeature

VariableFeatures(Tcell) <- setdiff(Tcell@assays$SCT@var.features, row.names(Tcell) %>% grep(pattern = "^MT-|^RPL|^RPS", v = T))

Tcell <- RunPCA(Tcell, features = VariableFeatures(Tcell), npcs = 50, verbose = TRUE)
Tcell <- RunHarmony(Tcell,group.by.vars = "sample")
#Tcell <- JackStraw(Tcell, num.replicate = 10,reduction = "pca",assay = "RNA")
ElbowPlot(Tcell,reduction = "harmony")

Tcell <- RunUMAP(Tcell,reduction = "harmony",dims = 1:30,min.dist = 0.5)
Tcell <- FindNeighbors(Tcell, reduction = "harmony", dims = 1:30)
Tcell <- FindClusters(Tcell, 
                      resolution =0.3,
                      method ='igraph',
                      verbose = T)
```

### scRNA analysis of B/Plasma
```
Bcell <- subset(PDAC_scRNA_sub,idents=c('B/Plasma'))
remove(PDAC_scRNA_sub)

# Bcell<-subset(Bcell, subset = PTPRC > 0,slot="counts")
#load("analysis/TME/Bcell.Rdata")
DotPlot(Bcell,features = c('EPCAM','KRT19','KRT18','KRT8','KRT7','CD79A','CD68','CD14','PECAM1','COL1A1','FN1'),group.by = "sample")

Bcell<-subset(Bcell, subset = EPCAM == 0,slot="counts")
Bcell<-subset(Bcell, subset = CD79A > 0,slot="counts")

Bcell <- SCTransform(Bcell,assay = "SCT")
Bcell <- FindVariableFeatures(Bcell,assay = "SCT")

varfeature <- Bcell@assays[["SCT"]]@var.features
varfeature <- varfeature[-grep('^TRA',varfeature)]
varfeature <- varfeature[-grep('^TRB',varfeature)]

Bcell@assays$SCT@var.features<-varfeature

VariableFeatures(Bcell) <- setdiff(Bcell@assays$SCT@var.features, row.names(Bcell) %>% grep(pattern = "^MT-|^RPL|^RPS", v = T))

Bcell <- RunPCA(Bcell, features = VariableFeatures(Bcell), npcs = 50, verbose = TRUE)
Bcell <- RunHarmony(Bcell,group.by.vars = "sample")
#Bcell <- JackStraw(Bcell, num.replicate = 10,reduction = "pca",assay = "RNA")
ElbowPlot(Bcell,reduction = "harmony")

Bcell <- RunUMAP(Bcell,reduction = "harmony",dims = 1:30,min.dist = 0.5)
Bcell <- FindNeighbors(Bcell, reduction = "harmony", dims = 1:30)
Bcell <- FindClusters(Bcell, 
                      resolution =0.15,
                      method ='igraph',
                      verbose = T)
```
