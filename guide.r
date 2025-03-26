YOUR SEURAT<-Read10X()
YOUR SEURAT <- CreateSeuratObject(counts = YOUR SEURAT, project = "YOUR SEURAT", min.cells = 3, min.features = 200)

data<-list( "YOUR SEURAT")
Seurat <- merge(x=data[[1]], y = data[-1], project = "harmony")
Seurat <- NormalizeData(Seurat, normalization.method = "LogNormalize", scale.factor = 10000)
Seurat <- FindVariableFeatures(Seurat, selection.method = "vst", nfeatures = 2000)
#QC
table(Seurat@meta.data$orig.ident)
Seurat[["percent.mt"]] <- PercentageFeatureSet(Seurat, pattern = "^MT-") 
metadata <-Seurat@meta.data
Idents(Seurat)<-"Sample_id"
names(Seurat@active.ident) <- colnames(Seurat)
print(dim(Seurat))
metadata <-Seurat@meta.data
#小提琴图质控
Idents(Seurat)<-"orig.ident"
p1<-VlnPlot(Seurat, features = c("nFeature_RNA"),y.max =6000,pt.size = 0)
p2<-VlnPlot(Seurat, features = c("nCount_RNA"),y.max = 12000,pt.size = 0)
p3<-VlnPlot(Seurat, features = c( "percent.mt"),y.max = 20,pt.size = 0)


Seurat <- subset(Seurat,subset = nFeature_RNA >700 & nFeature_RNA <15000  & percent.mt <5) 
dim(Seurat)

Seurat<- NormalizeData(Seurat) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE)
#下游分析
Seurat<- RunHarmony(Seurat, group.by.vars = "orig.ident")
Seurat<- RunUMAP(Seurat, reduction = "harmony", dims = 1:17, min.dist=0.2, return.model = TRUE)
Seurat <- FindNeighbors(Seurat, reduction = "harmony", dims = 1:18)
Seurat <- FindClusters(Seurat, resolution =0.75)
p4 <- DimPlot(Seurat, reduction = "umap", group.by = "seurat_clusters",label = TRUE)
p6 <- DimPlot(Seurat, reduction = "umap", label = TRUE, repel = TRUE)
p6

new.cluster.ids <- c("")#Celltype

names(new.cluster.ids) <- levels(Seurat)
Seurat <- RenameIdents(Seurat, new.cluster.ids)

#在metadata中，添加Celltype信息
Seurat$celltype <- Idents(Seurat)

head(Seurat@meta.data)
mycolor3<-c("#DCDEC8","#e5a7b4","#2681b6","#f2bf9e","#dbedd3","#8cba54","#97cebf", "#D8E9F7","#542788","#982c2c","#e1703c","#2d66a5","#d6b54c","#9d93e5","#aeceb6")
p2 <- DimPlot(Seurat, reduction = "umap", label = TRUE,cols=mycolor3,raster=F)
p2