library(pagoda2)
library(Matrix)
library(conos)


##  read raw  count matrix 

rawlist=readRDS('all.raw.list')

## run pagoda2 
p2list=lapply(raw,function(x) basicP2proc(x,n.cores = 10,min.transcripts.per.cell =0))



## using conos to integrate mutiple samples 

runConos=function(datlp2,appname,n.cores=8){
  print(names(datlp2))
  con <- Conos$new(datlp2,n.cores=n.cores)
  con$buildGraph()
  con$findCommunities()
  
  
  f1=paste(appname,'_conos_clustering.pdf',sep='')    # largvis embedding 
  p1=con$plotGraph()
  ggsave(f1,p1)
  
  f1=paste(appname,'_conos.rds',sep='')
  saveRDS(con,f1)
  
  
  con$embedGraph(method="UMAP", n.cores=30,min.dist=1e-10,n.neighbors=50)  # UMAP embedding 
  
  p2a <- con$plotGraph(alpha=0.1,size=0.1,plot.na=T)
  
  saveRDS(con$embedding,paste(appname,'.umap.rds',sep=''))
  
  ggsave(paste(appname,'.umap.png',sep=''),p2a,height=7,width=7)
  
  return(con)
  
}



con=runConos(p2list)

saveRDS(con,'con.rds')






# Seurat basic process  

library(Seurat)


pancreas.list=readRDS('../raw.list.all.rds')

print(names(pancreas.list))


pancreas.list=lapply(pancreas.list,function(x) CreateSeuratObject(counts = x, project = "BMET", min.cells = 3, min.features = 200))




for (i in 1:length(pancreas.list)) {
  pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = FALSE)
  pancreas.list[[i]] <- FindVariableFeatures(pancreas.list[[i]], selection.method = "vst",
                                             nfeatures = 2000, verbose = FALSE)
}

names(pancreas.list)

reference.list=pancreas.list
pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)



library(ggplot2)
library(cowplot)
library(patchwork)
# switch to integrated assay. The variable features of this assay are automatically
# set during IntegrateData
DefaultAssay(pancreas.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:30)

pancreas.integrated <- FindNeighbors(pancreas.integrated, reduction = "pca", dims = 1:30)
pancreas.integrated <- FindClusters(pancreas.integrated, resolution = 0.5)

saveRDS(pancreas.integrated,'Res.Seurat.rds')



