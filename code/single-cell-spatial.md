# Xenium
### POSTN Fibro Macro
```
library(Seurat)#V5
#sample PanIN
path <- "path to/Xenium/PanIN/PanIN1132"
# Load the Xenium data
xenium.obj <- LoadXenium(path, fov = "fov",assay = "Xenium")
# remove cells with 0 counts
xenium.obj <- subset(xenium.obj, subset = nCount_Xenium > 0)

VlnPlot(xenium.obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)
xenium.obj <- SCTransform(xenium.obj, assay = "Xenium")
xenium.obj <- RunPCA(xenium.obj, npcs = 30, features = rownames(xenium.obj))
xenium.obj <- RunUMAP(xenium.obj, dims = 1:30)
xenium.obj <- FindNeighbors(xenium.obj, reduction = "pca", dims = 1:30)
xenium.obj <- FindClusters(xenium.obj, resolution = 0.1)

ImageFeaturePlot(xenium.obj, features = c("POSTN","COL1A1"), blend=T,border.color = "black",axes = TRUE, max.cutoff = "q90")[[3]]
ImageFeaturePlot(xenium.obj, features = c('TFF1', 'MUC5AC'), blend=T,border.color = "black",axes = TRUE, max.cutoff = "q90")[[3]]

crop <- Crop(xenium.obj[["fov"]], x = c(2500, 4500), y = c(3000, 5000))
xenium.obj[["crop"]] <- crop
ImageDimPlot(xenium.obj,fov = "crop", cols = "polychrome", size = 1,border.color = "black",axes=T,crop = F)
#zoom 1
crop3 <- Crop(xenium.obj[["fov"]], x = c(3250,3750), y = c(3500,4000))
xenium.obj[["crop3"]] <- crop3
p<-ImageDimPlot(xenium.obj, boundaries='segmentation'  ,cols = c("yellow",col1[2],"orange","purple","green",col1[6],col1[1]), fov = "crop3",crop = F, axes=T,alpha = 0.3,border.size=0.1, size=1,mols.size=1.5,mols.cols=c("green",col1[11]),molecules = c("COL1A1","POSTN"))
p
p<-ImageDimPlot(xenium.obj, boundaries='segmentation' ,cols = c("yellow",col1[2],"orange","purple","green",col1[6],col1[1]), fov = "crop3",crop = F, axes=T,alpha = 0.3,border.size=0.1, size=1,mols.size=1.5,mols.cols=c("purple"),molecules = c("CD163"))
p
p<-ImageDimPlot(xenium.obj, boundaries='segmentation' ,cols = c("yellow",col1[2],"orange","purple","green",col1[6],col1[1]), fov = "crop3",crop = F, axes=T,alpha = 0.3,border.size=0.1, size=1,mols.size=1.5,mols.cols=c('yellow'),molecules = c("TFF1"))
p
p<-ImageDimPlot(xenium.obj, boundaries='segmentation' ,cols = c("yellow",col1[2],"orange","purple","green",col1[6],col1[1]), fov = "crop3",crop = F, axes=T,alpha = 0.3,border.size=0.1, size=1,mols.size=1.5,mols.cols=c("green",col1[11],"purple","yellow"),molecules = c("COL1A1","POSTN","CD163","TFF1"))
p


# zoom 2
crop4 <- Crop(xenium.obj[["fov"]], x = c(3000,3500), y = c(4000,4500))
xenium.obj[["crop4"]] <- crop4
p<-ImageDimPlot(xenium.obj, boundaries='segmentation'  ,cols = c("yellow",col1[2],"orange","purple","green",col1[6],col1[1]), fov = "crop4",crop = F, axes=T,alpha = 0.3,border.size=0.1, size=1,mols.size=1.5,mols.cols=c("green",col1[11]),molecules = c("COL1A1","POSTN"))
p
p<-ImageDimPlot(xenium.obj, boundaries='segmentation' ,cols = c("yellow",col1[2],"orange","purple","green",col1[6],col1[1]), fov = "crop4",crop = F, axes=T,alpha = 0.3,border.size=0.1, size=1,mols.size=1.5,mols.cols=c("purple"),molecules = c("CD163"))
p
p<-ImageDimPlot(xenium.obj, boundaries='segmentation' ,cols = c("yellow",col1[2],"orange","purple","green",col1[6],col1[1]), fov = "crop4",crop = F, axes=T,alpha = 0.3,border.size=0.1, size=1,mols.size=1.5,mols.cols=c('yellow'),molecules = c("TFF1"))
p
p<-ImageDimPlot(xenium.obj, boundaries='segmentation' ,cols = c("yellow",col1[2],"orange","purple","green",col1[6],col1[1]), fov = "crop4",crop = F, axes=T,alpha = 0.3,border.size=0.1, size=1,mols.size=1.5,mols.cols=c("green",col1[11],"purple","yellow"),molecules = c("COL1A1","POSTN","CD163","TFF1"))
p


```


### CCL4+ CD8+ T and IGHG1+ B Plasma
```
library(Seurat)#V5
options(future.globals.maxSize = 10 * 1024 * 1024 * 1024)
#sample1
path <- "path to/Xenium/Xenium_V1_Human_Ductal_Adenocarcinoma_FFPE_outs"
# Load the Xenium data
xenium.obj <- LoadXenium(path, fov = "fov",assay = "Xenium")
# remove cells with 0 counts
xenium.obj <- subset(xenium.obj, subset = nCount_Xenium > 0)

VlnPlot(xenium.obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)
ImageDimPlot(xenium.obj, fov = "fov", molecules = c("CD8A", "CD28"),crop = T)
xenium.obj <- SCTransform(xenium.obj, assay = "Xenium")
xenium.obj <- RunPCA(xenium.obj, npcs = 30, features = rownames(xenium.obj))
xenium.obj <- RunUMAP(xenium.obj, dims = 1:30)
xenium.obj <- FindNeighbors(xenium.obj, reduction = "pca", dims = 1:30)
xenium.obj <- FindClusters(xenium.obj, resolution = 0.1)

ImageFeaturePlot(xenium.obj, features = c("CD8A","CCL4"), blend=T,border.color = "black",axes = TRUE, max.cutoff = "q90")[[3]]
crop <- Crop(xenium.obj[["fov"]], x = c(1000, 3000), y = c(6000, 8000))
xenium.obj[["crop"]] <- crop
ImageDimPlot(xenium.obj,fov = "crop",group.by = "seurat_clusters", cols = "polychrome", size = 1,border.color = "black",axes=T,crop = F)

#Zoom 
#CCL4 CD8T and IGHG1 Plasma
#zoom 1
crop3 <- Crop(xenium.obj[["fov"]], x = c(1750,2000), y = c(6750, 7000))
xenium.obj[["crop3"]] <- crop3
p<-ImageDimPlot(xenium.obj, boundaries='segmentation' ,cols = c("orange","purple","green",col1[10]), fov = "crop3",crop = F, axes=T,alpha = 0.3,border.size=0.1, size=1,mols.size=1.5,mols.cols=c(col1[2],col1[2],col1[13],col1[11]),molecules = c("CD3D","CD3E","CD8A","CCL4"))
p
p2<-ImageDimPlot(xenium.obj, boundaries='segmentation' ,cols = c("orange","purple","green",col1[10]), fov = "crop3",crop = F, axes=T,alpha = 0.3,border.size=0.1, size=1,mols.size=1.5,mols.cols=c(col1[10],col1[18]),molecules = c("CD38","IGHG1"))
p2
#zoom 2
crop4 <- Crop(xenium.obj[["fov"]], x = c(1750,2000), y = c(6500, 6750))
xenium.obj[["crop4"]] <- crop4
p<-ImageDimPlot(xenium.obj, boundaries='segmentation' ,cols = c("orange","purple","green",col1[10]), fov = "crop4",crop = F, axes=T,alpha = 0.3,border.size=0.1, size=1,mols.size=1.5,mols.cols=c(col1[2],col1[2],col1[13],col1[11]),molecules = c("CD3D","CD3E","CD8A","CCL4"))
p
p2<-ImageDimPlot(xenium.obj, boundaries='segmentation' ,cols = c("orange","purple","green",col1[10]), fov = "crop4",crop = F, axes=T,alpha = 0.3,border.size=0.1, size=1,mols.size=1.5,mols.cols=c(col1[10],col1[18]),molecules = c("CD38","IGHG1"))
p2
#zoom 3
crop5 <- Crop(xenium.obj[["fov"]], x = c(2000,2250), y = c( 6750,7000))
xenium.obj[["crop5"]] <- crop5
p<-ImageDimPlot(xenium.obj, boundaries='segmentation' ,cols = c("orange","purple","green",col1[10]), fov = "crop5",crop = F, axes=T,alpha = 0.3,border.size=0.1, size=1,mols.size=1.5,mols.cols=c(col1[2],col1[2],col1[13],col1[11]),molecules = c("CD3D","CD3E","CD8A","CCL4"))
p
p2<-ImageDimPlot(xenium.obj, boundaries='segmentation' ,cols = c("orange","purple","green",col1[10]), fov = "crop5",crop = F, axes=T,alpha = 0.3,border.size=0.1, size=1,mols.size=1.5,mols.cols=c(col1[10],col1[18]),molecules = c("CD38","IGHG1"))
p2



```
