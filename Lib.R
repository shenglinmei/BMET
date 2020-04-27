library(reshape2)
library(xtable)
library(pagoda2)
library(nbHelpers)
library(reshape2)
library(Matrix)
#library(clusterMatch)
library(igraph)
library(ggplot2)
#library(nbHelpers)                  
library(parallel)
#library(fastSave)
library(pagoda2)
#library(nbHelpers)
library(reshape2)
library(pheatmap)
library(conos)
#library('betareg')
library(DESeq2)
library(Matrix)
library(conos)
library(GOstats)
library(RSQLite)
library(biomaRt)
suppressPackageStartupMessages(library(org.Hs.eg.db))


#  run pagoda apps 
runPagoda=function(cd,appname='none', n.cores = 1, batch = NULL, n.odgenes = 3000, nPcs = 100, 
                          k = 30, perplexity = 50, log.scale = TRUE, trim = 10, keep.genes = NULL, 
                          min.cells.per.gene = 30, get.largevis = TRUE, get.tsne = TRUE, 
                          make.geneknn = TRUE) {

    rownames(cd) <- make.unique(rownames(cd))
    p2 <- Pagoda2$new(cd, n.cores = n.cores, batch = batch, keep.genes = keep.genes, 
                      trim = trim, log.scale = log.scale, min.cells.per.gene = min.cells.per.gene)
    
    
    pv=paste(appname,'.adjustVariance.pdf',sep='')
    pdf(pv)
    p2$adjustVariance(plot = T, gam.k = 10)
    dev.off()
    
    p2$calculatePcaReduction(nPcs = nPcs, n.odgenes = n.odgenes, 
                             maxit = 1000)
    p2$makeKnnGraph(k = k, type = "PCA", center = TRUE, weight.type = "none", 
                    n.cores = n.cores, distance = "cosine")
    p2$getKnnClusters(method = igraph::infomap.community, type = "PCA", 
                      name = "infomap")
    p2$getKnnClusters(method = igraph::multilevel.community, 
                      type = "PCA", name = "multilevel")
    
    p2$getKnnClusters(method = igraph::walktrap.community, type = 'PCA', name='walktrap')
    
    p2$getDifferentialGenes(type='PCA',verbose=T,clusterType='multilevel')
    
    
    p2$getEmbedding(type = 'PCA', embeddingType = 'tSNE', perplexity = 50)
#    M <- 30
#    p2$getEmbedding(type = 'PCA', embeddingType = 'largeVis', M = M, perplexity = perplexity, gamma = 1 / M, alpha = 1 )
    
    
    pdf1=paste(appname,'.tsn.pdf',sep='') 
    pdf(pdf1)
    p2$plotEmbedding(type='PCA',embeddingType='tSNE',show.legend=F,mark.clusters=T,min.group.size=1,shuffle.colors=F,mark.cluster.cex=1,alpha=0.1,main='clusters (tSNE)')
    dev.off()
    
    
    
    return(p2)
}



##    merge data from  list oblects
mergeDat=function(raw.mats){
  
  genelists <- lapply(raw.mats, function(x) rownames(x))
  str(genelists)
  commongenes <- Reduce(intersect,genelists)
  
  
  matrices2 <- mapply(function(m, name) {
    colnames(m) <- paste(name, colnames(m), sep='_');
    m[commongenes,]
  },
  raw.mats,
  names(raw.mats))
  cellNum=unlist(lapply(matrices2, function(x) ncol(x)))
  print('#commongenes')
  print(length(commongenes))
  
  print('#cell numbers')
  print(cellNum)
  
  return(matrices2)
}

#matrices=mergeDat(raw.mats)

creatLableList=function(matrices,res,key){
#  key='allName'
  celllists <- lapply(matrices, function(x) colnames(x))
  celss=unlist(celllists)
  names(celss)=celss
  
  for ( i in seq(length(matrices))){
    sampID=as.character(res[i,1])
    value=as.character(res[i,key])
    celss[celllists[[sampID]]]=value
  }
  print(table(celss))
  return(celss)
}

#ceLable=creatLableList(matrices,res,'allName')




makeWebLable2=function(appName,p2,jfac2,cell_ano2,n.cores=5,save=NULL){
  
  jfac2=as.factor(jfac2)
  cells.in.app <- rownames(p2$counts,jfac2)
  ## make hierarchical aspects
  hdea <- p2$getHierarchicalDiffExpressionAspects(type='PCA',clusterName='multilevel',z.threshold=3, n.cores = 5)
  ## make metadata   
  metadata.forweb <- list();
  metadata.forweb$multilevel <- p2.metadata.from.factor(p2$clusters$PCA$multilevel,displayname='Multilevel')
  ## Set colors manually
  #  myinfo=p2$clusters$PCA$infomap
  #  p1 <- rainbow(length(levels(myinfo)), s=1, v=1)
  #  names(p1) <- levels(myinfo)
  #  metadata.forweb$infomap <- p2.metadata.from.factor(myinfo[cells.in.app], pal=p1,displayname='Infomap')
  metadata.forweb$infomap <- p2.metadata.from.factor(p2$clusters$PCA$infomap, displayname = 'Infomap', start=0, end=0.5, s = 1, v=0.7)
  
  
  p1 <- rainbow(length(levels(jfac2)), s=1, v=1)
  names(p1) <- levels(jfac2)
  metadata.forweb$Label <- p2.metadata.from.factor(jfac2[cells.in.app], pal=p1,displayname='Label')
  ## get de sets
  
  if(!is.null(cell_ano2)){
    cell_ano2=as.factor(cell_ano2)
    p1 <- rainbow(length(levels(cell_ano2)), s=1, v=1)
    names(p1) <- levels(cell_ano2)
    metadata.forweb$cell_ano2 <- p2.metadata.from.factor(cell_ano2[cells.in.app], pal=p1,displayname='cell_ano2')
  } 
  
  deSets <- get.de.geneset(p2, groups=p2$clusters$PCA$multilevel, prefix='de_')
  ## Collect the genesets
  genesets = c(deSets, hierDiffToGenesets(hdea))
  
  #  genesets = hierDiffToGenesets(hdea)
  appmetadata = list(apptitle=appName)
  p2$makeGeneKnnGraph();
  ## make app
  wp <- make.p2.app(p2, additionalMetadata = metadata.forweb, geneSets = genesets,
                    dendrogramCellGroups=p2$clusters$PCA[[1]],show.clusters=F,
                    appmetadata=appmetadata)
  wp$serializeToStaticFast(paste0(appName,'.bin'))
  if(!is.null(save)){
    fout=paste(appName,'.rds',sep='')
    saveRDS(p2,fout)
  }
}



makeWebLable=function(appName,p2,jfac2,n.cores=5){

  jfac2=as.factor(jfac2)
  cells.in.app <- rownames(p2$counts,jfac2)
  ## make hierarchical aspects
  hdea <- p2$getHierarchicalDiffExpressionAspects(type='PCA',clusterName='multilevel',z.threshold=3, n.cores = 5)
  ## make metadata   
  metadata.forweb <- list();
  metadata.forweb$multilevel <- p2.metadata.from.factor(p2$clusters$PCA$multilevel,displayname='Multilevel')
  ## Set colors manually
  #  myinfo=p2$clusters$PCA$infomap
  #  p1 <- rainbow(length(levels(myinfo)), s=1, v=1)
  #  names(p1) <- levels(myinfo)
  #  metadata.forweb$infomap <- p2.metadata.from.factor(myinfo[cells.in.app], pal=p1,displayname='Infomap')
  metadata.forweb$infomap <- p2.metadata.from.factor(p2$clusters$PCA$infomap, displayname = 'Infomap', start=0, end=0.5, s = 1, v=0.7)
  
  
  p1 <- rainbow(length(levels(jfac2)), s=1, v=1)
  names(p1) <- levels(jfac2)
  metadata.forweb$Label <- p2.metadata.from.factor(jfac2[cells.in.app], pal=p1,displayname='Label')
  ## get de sets
  
  #    p1 <- rainbow(length(levels(S_cells_ano)), s=1, v=1)
  #    names(p1) <- levels(S_cells_ano)
  
  #    metadata.forweb$jointsamp <- p2.metadata.from.factor(S_cells_ano[cells.in.app], pal=p1,displayname='jointsamp')
  
  deSets <- get.de.geneset(p2, groups=p2$clusters$PCA$multilevel, prefix='de_')
  ## Collect the genesets
  genesets = c(deSets, hierDiffToGenesets(hdea))
  
  #  genesets = hierDiffToGenesets(hdea)
  appmetadata = list(apptitle=appName)
  p2$makeGeneKnnGraph();
  ## make app
  wp <- make.p2.app(p2, additionalMetadata = metadata.forweb, geneSets = genesets,
                    dendrogramCellGroups=p2$clusters$PCA[[1]],show.clusters=F,
                    appmetadata=appmetadata)
  wp$serializeToStaticFast(paste0(appName,'.bin'))
  
  fout=paste(appName,'.rds',sep='')
  saveRDS(p2,fout)
}




makeWeb=function(p2,appName,n.cores=5){
  
  hdea <- p2$getHierarchicalDiffExpressionAspects(type = "PCA", 
                                                  clusterName = "multilevel", z.threshold = 3, n.cores = 5)
  metadata.forweb <- list()
  metadata.forweb$multilevel <- p2.metadata.from.factor(p2$clusters$PCA$multilevel, 
                                                        displayname = "Multilevel")
  
  metadata.forweb$infomap <- p2.metadata.from.factor(p2$clusters$PCA$infomap, displayname = 'Infomap', start=0, end=0.5, s = 1, v=0.7)
  #    metadata.forweb$walktrap <- p2.metadata.from.factor(p2$clusters$PCA$walktrap, displayname = 'Walktrap', s = 0.5)
  extraWebMetadata = NULL
  metadata.forweb <- c(metadata.forweb, extraWebMetadata)
  
  deSets <- get.de.geneset(p2, groups = p2$clusters$PCA[[1]], prefix = 'de_')
  
  
  genesets <- hierDiffToGenesets(hdea)
  
  genesets <- deSets
  
  
  appmetadata = list(apptitle = appName)
  cat("Making KNN graph...\\n")
  p2$makeGeneKnnGraph(n.cores = n.cores)
  myPagoda2WebObject=make.p2.app(p2, additionalMetadata = metadata.forweb,
                                 geneSets = genesets, 
                                 dendrogramCellGroups = p2$clusters$PCA$multilevel,
                                 show.clusters = F, 
                                 appmetadata = appmetadata)
  
  
  # Save serialised web object, RDS app and session image
  #myPagoda2WebObject$serialiseToStatic(text.file.directory = './tmp', binary.filename = paste0(appName,'.bin'))
  myPagoda2WebObject$serializeToStaticFast(binary.filename = paste0(appName,'.bin'))
  
  #    invisible(p2)
  
  
}



GOanalysis=function(markers,n){

  ENTREZID=unlist(mget(markers, org.Hs.egSYMBOL2EG, ifnotfound=NA))
  ENTREZID=ENTREZID[!is.na(ENTREZID)]
  
  
  for(function_type in c("BP", "CC", "MF")){
    
    param <- new("GOHyperGParams", geneIds=ENTREZID,
                 #universe=universe,
                 annotation="org.Hs.eg.db", ontology=function_type,pvalueCutoff=0.1,
                 conditional=FALSE, testDirection="over")
    hyp <- hyperGTest(param)
    sumTable <- summary(hyp)
    
    
    
    david=sumTable[1:20,]
    david$Pvalue=-log(david[,2])
    termNumber=nrow(david)
    
    
    library(ggplot2)
    p1 <- ggplot(data=david, aes(x=Pvalue, y=Term, size=Count, colour = factor(david$Count)))
    p1 <- p1 + geom_point()
    p1 <- p1 + guides(color = FALSE)
    p1 <- p1 + theme(panel.grid.major = element_line(colour='blue'),
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank())
    p1 <- p1 + xlab(paste("-log10(Pvalue)", sep="")) + ylab("")
    p1 <- p1 + labs(title=paste("DAVID:", function_type, sep=""))
    p1 <- p1 + theme(axis.text.x=element_text(size=10, face="plain", colour ='black'))
    p1 <- p1 + theme(axis.text.y=element_text(size=6, face="plain", colour ='black'))
    #p1=p1+theme(axis.title.y=element_text(size=10,face="plain",colour ='black'))
    #p1=p1+theme(legend.position="bottom")
    p1 <- p1 + xlim(min(david$Pvalue), max(david$Pvalue))
    #p1=p1+ylim(-10,+15)
    print(p1)
    ggsave(file=paste(n,'_' ,function_type, ".png", sep=""), scale=0.8, dpi=600, width = 7, height=1+0.25*termNumber) 
  } 
  
}







#  draw figure for single genes
#plotWithGroupsGene(p2ens$p2objs,'PLAUR',jcl3.coarse, file='Myoloid_emb.jcl3.PLAUR.png',verbose=T,panel.size=300)
plotWithGroupsGene=function(p2objs,gene,groups = NULL, filename = NULL, panel.size = 300, mark.cluster.cex = 0.8, 
                            mar = c(0.5, 0.5, 0.5, 0.5), mgp = c(2, 0.65, 0), cex = 0.85, 
                            type = "PCA", embeddingType = "tSNE", mark.clusters = TRUE, 
                            verbose = FALSE) 
{
  require(Cairo)
  panel.dims <- getParMfrow(length(p2objs))
  if (verbose) 
    cat("Panel dimensions are ", panel.dims[1], " by ", panel.dims[2], 
        "\\n")
  if (!is.null(filename)) 
    CairoPNG(file = filename, height = panel.dims[1] * panel.size, 
             width = panel.dims[2] * panel.size)
  par(mfrow = c(panel.dims[1], panel.dims[2]), mar = mar, mgp = mgp, 
      cex = cex)
  lapply(names(p2objs), function(dn) {
    d <- p2objs[[dn]]
    g1 <- as.factor(groups)
    colors <- NULL
    if (!any(names(g1) %in% rownames(d$counts))) {
      g1 <- NULL
      cell.names <- rownames(d$counts)
      colors <- rep("grey70", length(cell.names))
      names(colors) <- cell.names
    }
    d$plotEmbedding(type = type, embeddingType = embeddingType, 
                    alpha = 0.2, min.group.size = 0, mark.clusters = mark.clusters, 
                    mark.cluster.cex = mark.cluster.cex, do.par = F, 
                    colors = d$counts[,gene])
    
    
    
    legend(x = "topleft", bty = "n", legend = dn)
  })
  if (!is.null(filename)) 
    dev.off()
  invisible(NULL)
}





runP2Genes=function(p2,gs,nPcs=12,k=15,n.cores=10 ,alpha = 0.3 ){
  p2$adjustVariance(plot = F, gam.k = 10,alpha =alpha)
  g2=intersect(p2$misc$odgenes,gs)
  print(length(g2))
  p2$calculatePcaReduction(nPcs = 12, odgenes = g2, 
                           maxit = 1000)
  p2$makeKnnGraph(k = k, type = "PCA", center = TRUE, weight.type = "none", 
                  n.cores = n.cores, distance = "cosine")                             
  p2$getKnnClusters(method = igraph::multilevel.community, 
                    type = "PCA", name = "multilevel")
  p2$getEmbedding(type = 'PCA', embeddingType = 'tSNE', perplexity = 50)
  return(p2)
}


draw_conos_sgdold=function(con,appname,sgd_batches,id,cell_ano=NULL){
  con$embedGraph(sgd_batches = sgd_batches)
  f1=paste(appname,'_',id,'_conos_clustering.png',sep='')
  a=con$plotGraph(groups=cell_ano_sampleType)
  ggsave(f1,a,width = 7,height=7)
  
  if(is.null(cell_ano)){
    f1=paste(appname,'_',id,'_conos_clustering2.png',sep='')
    a2= con$plotGraph(groups=cell_ano)
    ggsave(f1,a2,width = 7,height=7)
    f1=paste(appname,'_',id,'_conos_clustering2.rds',sep='')
    saveRDS(con,f1)
  } 
  
}

#draw_fig(con,appname2,1e+09,109)
#draw_fig(con,appname2,3e+08,308)
#draw_fig(con,appname2,6e+08,608)
#draw_fig(con,appname2,2e+09,209)
#draw_fig(con,appname2,4e+09,409)
#draw_fig(con,appname2,9e+09,909)








# drawfigureP2(p2,appname,jcl3.coarse=NULL,cell_ano_sampleType=NULL,cell_ano_sample=NULL,saveRDS=NULL)

drawfigureP2=function(p2,appname,jcl3.coarse=NULL,cell_ano_sampleType=NULL,cell_ano_sample=NULL,saveRDS=NULL,alpha=0.1){
  pdf1=paste(appname,'.clustering.tsn.png',sep='') 
  png(pdf1,height=600,width=600)
  p2$plotEmbedding(groups=p2$clusters$PCA$multilevel,type='PCA',embeddingType='tSNE',show.legend=F,mark.clusters=T,min.group.size=1,shuffle.colors=F,mark.cluster.cex=1,alpha=alpha,main='clusters (tSNE)')
  dev.off()
  
  
  if (!is.null(jcl3.coarse)){
    pdf1=paste(appname,'.cells.tsn.png',sep='') 
    png(pdf1,height=600,width=600)
    p2$plotEmbedding(type='PCA',groups =jcl3.coarse,embeddingType='tSNE',show.legend=F,mark.clusters=T,min.group.size=1,shuffle.colors=F,mark.cluster.cex=1,alpha=alpha,main='clusters (tSNE)')
    dev.off()
  }
  
  
  if (!is.null(cell_ano_sampleType)){ 
    pdf1=paste(appname,'.cells.tsn.sampleType.png',sep='') 
    png(pdf1,height=600,width=600)
    p2$plotEmbedding(type='PCA',groups =cell_ano_sampleType,embeddingType='tSNE',show.legend=F,mark.clusters=T,min.group.size=1,shuffle.colors=F,mark.cluster.cex=1,alpha=alpha,main='clusters (tSNE)')
    dev.off()
  }
  
  
  
  if (!is.null(cell_ano_sample)){
    pdf1=paste(appname,'.cells.tsn.sample.png',sep='') 
    png(pdf1,height=600,width=600)
    p2$plotEmbedding(type='PCA',groups =cell_ano_sample,embeddingType='tSNE',show.legend=F,mark.clusters=T,min.group.size=1,shuffle.colors=F,mark.cluster.cex=1,alpha=alpha,main='clusters (tSNE)')
    dev.off()
  }
  
  
  if (!is.null(saveRDS)){
    f1=paste(appname,'_p2combined.rds',sep='')
    saveRDS(p2,f1)
  } 
  
}


#drawfigureP2(p2,appname,jcl3.coarse=cell_ano_cell,cell_ano_sampleType=cell_ano_sampleType,cell_ano_sample=cell_ano_sample,saveRDS=NULL)
  



geneExpP2=function(p2,gs,lcol,lrow,appanme){
#  gs3=c('CD8A','CD8B','CD4','FCGR3A','FGFBP2','CD3E','CD3D','SH2D1A','TNFRSF4','CTLA4','TIGIT','GZMA','GZMB','GZMM','GZMK','GZMH','LAG3')
  filename=paste(appname,'exp.png',sep='')
  mar = c(0.5, 0.5, 0.5, 0.5)
  mgp = c(2, 0.65, 0)
  panel.size=300
  cex=1.1
  mark.cluster.cex = 0.8
  mark.clusters = TRUE
  
  png(file = filename, height = lrow * panel.size, 
      width = lcol * panel.size)
  par(mfrow = c(lrow, lcol), mar = mar, mgp = mgp, 
      cex = cex)
  
  for (gene in gs){
    p2$plotEmbedding(type='PCA',embeddingType='tSNE',colors=p2$counts[,gene],shuffle.colors=F,mark.cluster.cex=1,alpha=0.2,main=gene)
  }
  dev.off()
  
}

#gs=c('CD38','IGLL1','SOX4','CD9','CD52')
#geneExpP2(p2,gs,3,2,appanme)


geneExpConos=function(con,p2,gs,lrow,lcol,appanme,alpha =0.3){
  lis=list()
  for (gene in gs){
    print(gene)
    t=p2$counts[,gene]
    a=con$plotGraph(colors =t,title=gene,alpha =0.3)
    lis[[gene]]=a
  }
  b=  cowplot::plot_grid(plotlist=lis, ncol=lrow, nrow=lcol)
  fout=paste(appname,'expConos.png',sep='')
  ggsave(fout,b,width = 2.3*lcol,height=2.3*lrow)
}

#gs=c('CD38','CD8A','CD8B')
#geneExpConos(con,p2,gs,3,2,appanme)





#drawfigureConos=function(p2,appname,jcl3.coarse=NULL,cell_ano_sampleType=NULL,cell_ano_sample=NULL,saveRDS=NULL)
drawfigureConos=function(p2,appname,jcl3.coarse=NULL,cell_ano_sampleType=NULL,cell_ano_sample=NULL,saveRDS=NULL){
  pdf1=paste(appname,'.clustering.tsn.conos.png',sep='') 
  a1=con$plotGraph()
  ggsave(pdf1,a1,width = 7,height=7)
  
  
  if (!is.null(jcl3.coarse)){
    pdf1=paste(appname,'.cells.conos.png',sep='') 
    a1=con$plotGraph(groups =jcl3.coarse,show.legend=TRUE,title=appname)
    ggsave(pdf1,a1,width = 7,height=7)
  }
  
  
  if (!is.null(cell_ano_sampleType)){ 
    pdf1=paste(appname,'.sampleType.conos.png',sep='') 
    a1=con$plotGraph(groups =cell_ano_sampleType,show.legend=TRUE,title=appname)
    ggsave(pdf1,a1,width = 7,height=7)
  }
  
  
  
  if (!is.null(cell_ano_sample)){
    pdf1=paste(appname,'.sample.conos.png',sep='') 
    a1=con$plotGraph(groups =cell_ano_sample,show.legend=TRUE,title=appname)
    ggsave(pdf1,a1,width = 9,height=7)
  }
  
  
  if (!is.null(saveRDS)){
    f1=paste(appname,'_conos.rds',sep='')
    saveRDS(p2,f1)
  } 
  
}

#drawfigureConos(con,appname,jcl3.coarse=cell_ano_cell,cell_ano_sampleType=cell_ano_sampleType,cell_ano_sample=cell_ano_sample,saveRDS=NULL)






# DEcaculate(p2,appname,conosCluster,removeGene=TRUE)
# caculate differnential expressed gene based on cluster group
DEcaculate=function(p2,appname,conosCluster,removeGene=NULL){
  
  de1 <- p2$getDifferentialGenes(groups=conosCluster)
  
  
  f1=paste(appname,'_diffGene.rds',sep='')
  saveRDS(de1,f1)
  for(n in names(de1)){
    x=de1[[n]]
    z <- x[order(-x$Z),]
    if(!is.null(removeGene)){
      index=grepl('^RP[LKS]',rownames(z))
      z=z[!index,]
    }
    markers=rownames(z)[1:100]
    x <- as.matrix(p2$counts[names(conosCluster),markers])
    ## trim outliers
    x <- apply(x, 2, function(xp) {
      qs <- quantile(xp,c(0.01,0.98))
      xp[xp<qs[1]] <- qs[1]
      xp[xp>qs[2]] <- qs[2]
      xp
    })
    x <- x[,(apply(x,2,sd) != 0)]
    x <- t(scale(x))
    ## sample 2000 cells for plotting
    #x <- x[,sample(ncol(x),2000)]
    o <- order(as.numeric(as.character(conosCluster[colnames(x)])))
    annot <- data.frame(cluster=conosCluster[colnames(x)],row.names = colnames(x))
    
    annot$dtype='other'
    annot[as.character(annot[,1])==n,'dtype']=n
    annot$dtype=as.factor(annot$dtype)
    
    
    pal <- colorRampPalette(c('navy','white','firebrick3'))(50)
    ## draw heatmap
    
    fout=paste(appname,'_',n,'_marker.heatmap.new.png',sep='')
    rgb.palette <- colorRampPalette(c("blue","white","red"), space = "rgb" )
    pheatmap(x[,o],labels_col=FALSE,cluster_cols=FALSE,annotation_col = annot,show_colnames = F,
             cluster_rows = FALSE,color=rgb.palette(100),filename=fout,fontsize_row =4,width=5,height=6,   #3.5*0.02*length(markers),
             breaks = c(seq(min(x),-0.01,length.out = 50),seq(0.01,max(x),length.out = 50)))
    
    
    iid= paste(appname,'_cluster_',n,sep='')
    GOanalysis(markers,iid) 
  }

}

#conosCluster=con$clusters$multi$groups
#DEcaculate(p2,appname,conosCluster,removeGene=TRUE)






# MakeWebConos(p2,appanme,conosTSN,conosCluster,jcl3.coarse2)
#conosTSN=con$embedding
#conosCluster=con$clusters$multi$groups


MakeWebConos<-function(p2,appname,conosTSN,conosCluster,jcl3.coarse2,cell_ano_sample=NULL,cell_ano_sampleType=NULL,combine=NULL,grade=NULL){
  n=appname
  p2$embeddings$PCA$conosEb=t(conosTSN)
  jfac2=as.factor(conosCluster)
  
  S_cells_ano=as.factor(jcl3.coarse2)  
  
  cells.in.app <- rownames(p2$counts,jfac2)
  ## make hierarchical aspects
  hdea <- p2$getHierarchicalDiffExpressionAspects(type='PCA',clusterName='multilevel',z.threshold=3)
  ## make metadata   
  metadata.forweb <- list();
  metadata.forweb$multilevel <- p2.metadata.from.factor(p2$clusters$PCA$multilevel,displayname='Multilevel')
  ## Set colors manually
  ##p1 <- colorRamps::primary.colors(n = nlevels(jfac2))
  p1 <- rainbow(length(levels(jfac2)), s=1, v=1)
  names(p1) <- levels(jfac2)
  metadata.forweb$conos <- p2.metadata.from.factor(jfac2[cells.in.app], pal=p1,displayname='conos')
  ## get de sets
  
  
  p1 <- rainbow(length(levels(S_cells_ano)), s=1, v=1)
  names(p1) <- levels(S_cells_ano)
  metadata.forweb$jointsamp <- p2.metadata.from.factor(S_cells_ano[cells.in.app], pal=p1,displayname='jointsamp')
  
  if(!is.null(cell_ano_sample)){
    cell_ano_sample=as.factor(cell_ano_sample)
    p1 <- rainbow(length(levels(cell_ano_sample)), s=1, v=1)
    names(p1) <- levels(cell_ano_sample)
    metadata.forweb$sample <- p2.metadata.from.factor(cell_ano_sample[cells.in.app], pal=p1,displayname='sample')
  } 
  
  if(!is.null(cell_ano_sample)){
    cell_ano_sampleType=as.factor(cell_ano_sampleType)
    p1 <- rainbow(length(levels(cell_ano_sampleType)), s=1, v=1)
    names(p1) <- levels(cell_ano_sampleType)
    metadata.forweb$sampleType <- p2.metadata.from.factor(cell_ano_sampleType[cells.in.app], pal=p1,displayname='sampleType')
  } 
  
  
  if(!is.null(combine)){
    cell=jcl3.coarse2
    sampleType=cell_ano_sampleType[names(jcl3.coarse2)]
    
    ano_combin=paste(cell,as.character(sampleType),sep='_')
    names(ano_combin)=names(jcl3.coarse2)
    
    ano_combin=as.factor(ano_combin)
    
    print(table(ano_combin))
    p1 <- rainbow(length(levels(ano_combin)), s=1, v=1)
    names(p1) <- levels(ano_combin)
    metadata.forweb$CellType_DiseaseStatus <- p2.metadata.from.factor(ano_combin[cells.in.app], pal=p1,displayname='Cell_DiseaseStatus')
  } 
  
   if(!is.null(grade)){
    grade=as.factor(grade)
    p1 <- rainbow(length(levels(grade)), s=1, v=1)
    names(p1) <- levels(grade)
    metadata.forweb$status <- p2.metadata.from.factor(grade[cells.in.app], pal=p1,displayname='status')
  } 
  
   
  deSets <- get.de.geneset(p2, groups=jfac2, prefix='conosDE_')
  ## Collect the genesets
  genesets = c(deSets, hierDiffToGenesets(hdea))
  appmetadata = list(apptitle=n)
  p2$makeGeneKnnGraph();
  ## make app
  wp <- make.p2.app(p2, additionalMetadata = metadata.forweb, geneSets = genesets,
                    dendrogramCellGroups=p2$clusters$PCA[[1]],show.clusters=F,
                    appmetadata=appmetadata)
  wp$serializeToStaticFast(paste0(n,'.bin'))
  
}  


#conosTSN=con$embedding
#conosCluster=con$clusters$multi$groups
#MakeWebConos(p2,appname,conosTSN,conosCluster,cell_ano_cell,cell_ano_sample=cell_ano_sample,cell_ano_sampleType=cell_ano_sampleType)
  



creatLableList2=function(jcl3.coarse){
  #  key='allName'
  cell_ano_cell=as.character(jcl3.coarse)
  names(cell_ano_cell)=names(jcl3.coarse)
  cell_ano_sample=names(cell_ano_cell)
  names(cell_ano_sample)=names(cell_ano_cell)
  cell_ano_sample=apply(data.frame(cell_ano_sample),1,function(x) strsplit(x[1],'_')[[1]][1])
  cell_ano_sampleType=apply(data.frame(cell_ano_sample),1,function(x) strsplit(x[1],'-')[[1]][2])
  cellano=list('cell'=cell_ano_cell,'sample'=cell_ano_sample,'sampltType'=cell_ano_sampleType)
  print(table(cell_ano_cell))
  print(table(cell_ano_sampleType))
  print(table(cell_ano_sample))
  
  return(cellano)
}


mergeDat2=function(raw.mats){
  
  genelists <- lapply(raw.mats, function(x) rownames(x))
  str(genelists)
  commongenes <- Reduce(intersect,genelists)
  
  
  matrices2 <- mapply(function(m, name) {
 #   colnames(m) <- paste(name, colnames(m), sep='_');
    m[commongenes,]
  },
  raw.mats,
  names(raw.mats))
  cellNum=unlist(lapply(matrices2, function(x) ncol(x)))
  print('#commongenes')
  print(length(commongenes))
  
  print('#cell numbers')
  print(cellNum)
  bigM2 <- Reduce(cbind, matrices2)
  
  return(bigM2)
}





draw_conos_sgd=function(con,appname,sgd_batches,id,cell_ano=NULL){
  con$embedGraph(sgd_batches = sgd_batches)
  f1=paste(appname,'_',id,'_conos_clustering.png',sep='')
  a=con$plotGraph()
  ggsave(f1,a,width = 7,height=7)
  
  if(is.null(cell_ano)){
    f1=paste(appname,'_',id,'_conos_clustering2.png',sep='')
    a2= con$plotGraph(groups=cell_ano)
    ggsave(f1,a2,width = 7,height=7)
    #   f1=paste(appname,'_',id,'_conos_clustering2.rds',sep='')
    #  saveRDS(con,f1)
  } 
  
}


# Run DEseq 
# @input count table and categrey tab1 
# count1=do.call(cbind,lapply(raw.mats,function(x) rowSums(x)))
# tab1=apply(data.frame(colnames(count1)),1,function(x) strsplit(x,'-')[[1]][2])

runDeseq=function(count1,tab1) {
  
  ntype=table(tab1)
  lres=list()
  for (cell in names(ntype[ntype>1])){
    
    cm=count1
    ntab=tab1
    ntab[ntab!=cell]=0
    ntab[ntab==cell]=1
    meta <- data.frame(
      sample.id= colnames(cm),
      group= ntab
    )
    print(meta)
    dds1 <- DESeqDataSetFromMatrix(cm,meta,design=~group)
    dds1 <- DESeq(dds1)
    res1 <- results(dds1,independentFiltering = TRUE,pAdjustMethod = "BH")
    res1=res1[!is.na(res1$pvalue),]
    res1=res1[order(res1$pvalue),]
    lres[[cell]]=res1
    
  }
  return(lres)  
}  



# Run differential expressed gene with pagoda2
 # Run differential expressed gene with pagoda2
DEcaculate2=function(p2,appname,conosCluster,removeGene=NULL,cutoff=3,num=100,GO=NULL){
  allg=NULL  
  de1 <- p2$getDifferentialGenes(groups=conosCluster,z.threshold = cutoff)
  
  
  f1=paste(appname,'_diffGene.rds',sep='')
  saveRDS(de1,f1)
  for(n in names(de1)){
    x=de1[[n]]
    z <- x[order(-x$Z),]
    if(!is.null(removeGene)){
      index=grepl('^RP[LKS]',rownames(z))
      z=z[!index,]
    }
    markers=rownames(z)[1:num]
    markers=markers[!is.na(markers)]
    allg=cbind(allg,markers)
    x <- as.matrix(p2$counts[names(conosCluster),markers])
    ## trim outliers
    x <- apply(x, 2, function(xp) {
      qs <- quantile(xp,c(0.01,0.98))
      xp[xp<qs[1]] <- qs[1]
      xp[xp>qs[2]] <- qs[2]
      xp
    })
    x <- x[,(apply(x,2,sd) != 0)]
    x <- t(scale(x))
    ## sample 2000 cells for plotting
    #x <- x[,sample(ncol(x),2000)]
    o <- order(conosCluster[colnames(x)])
    annot <- data.frame(cluster=conosCluster[colnames(x)],row.names = colnames(x))
    
    annot$dtype='other'
    annot[as.character(annot[,1])==n,'dtype']=n
    annot$dtype=as.factor(annot$dtype)
    
    
    pal <- colorRampPalette(c('navy','white','firebrick3'))(50)
    ## draw heatmap
    
    fout=paste(appname,'_',n,'_marker.heatmap.new.png',sep='')
    rgb.palette <- colorRampPalette(c("blue","white","red"), space = "rgb" )
    pheatmap(x[,o],labels_col=FALSE,cluster_cols=FALSE,annotation_col = annot,show_colnames = F,
             cluster_rows = FALSE,color=rgb.palette(100),filename=fout,fontsize_row =4,width=5,height=8,   #3.5*0.02*length(markers),
             breaks = c(seq(min(x),-0.01,length.out = 50),seq(0.01,max(x),length.out = 50)))
    
    
    iid= paste(appname,'_cluster_',n,sep='')
    if(!is.null(GO)){
      GOanalysis(markers,iid) }
  }
  colnames(allg)=names(de1)
  
  
  write.table(allg,paste(appname,'_Diff_gene.xls',sep=''),sep='\t',col.names=T,row.names=F,quote=F)
  
}  

 
 




# draw mutiple figure in the one page conos 

drawfigureConos2=function(p2,appname,jcl3.coarse=NULL,cell_ano_sampleType=NULL,cell_ano_sample=NULL,saveRDS=NULL){

  lis=list()
  if (!is.null(jcl3.coarse)){
    pdf1=paste(appname,'.cells.conos.png',sep='') 
    a1=con$plotGraph(groups =jcl3.coarse,show.legend=TRUE,title=appname)
#    ggsave(pdf1,a1,width = 7,height=7)
    lis[[1]]=a1
  }
  
  
  if (!is.null(cell_ano_sampleType)){ 
    pdf1=paste(appname,'.sampleType.conos.png',sep='') 
    a2=con$plotGraph(groups =cell_ano_sampleType,show.legend=TRUE,title=appname)
#    ggsave(pdf1,a1,width = 7,height=7)
    lis[[2]]=a2
  }
  
  
  
  if (!is.null(cell_ano_sample)){
    pdf1=paste(appname,'.sample.conos.png',sep='') 
    a3=con$plotGraph(groups =cell_ano_sample,show.legend=TRUE,title=appname)
#    ggsave(pdf1,a1,width = 9,height=7)
    lis[[3]]=a3
  }
  
  
  if (!is.null(saveRDS)){
    f1=paste(appname,'_conos.rds',sep='')
    saveRDS(p2,f1)
  } 
  
    b=  cowplot::plot_grid(plotlist=lis, ncol=length(lis), nrow=1)

	return(b)
  
}

#b=drawfigureConos2(con,appname,jcl3.coarse=anoCell,cell_ano_sampleType=anoSampleType,cell_ano_sample=anoCell2,saveRDS=NULL)

#  Run robust ranking
library(RobustRankAggreg)
robust_ranking=function(res,decreasing=NULL) {
  
  DataOrder = list()
  DataOrderName = list()
  nset=ncol(res)
  for (i in seq(nset)){
    order_tmp =order(as.numeric(as.character(res[,i])),decreasing=F)
    if (!is.null(decreasing)){ order_tmp =order(as.numeric(as.character(res[,i])),decreasing=TRUE) }
    tmp_name = rownames(res[order_tmp,])
    print(i)
    print(tmp_name[1:4])
    print(order_tmp[1:4])
    DataOrderName[[i]]=tmp_name
    DataOrder[[i]] = rank(as.numeric(as.character(res[,i])))}
  
  glist = 	DataOrderName
  r = RobustRankAggreg::rankMatrix(glist) 
  roubustRanking = aggregateRanks(rmat = r)
  return(roubustRanking)
}



listToOverlapMatrix=function(res,counts=NULL){
  
  
  lc=length(res)
  stat=matrix(rep(0,lc*lc),lc,lc)
  for( i in seq(lc)){
    for (j in seq(lc)){
      stat[i,j]=length(intersect(res[[i]],res[[j]]))/length(union(res[[i]],res[[j]]))
      if (!is.null(counts)){stat[i,j]=length(intersect(res[[i]],res[[j]])) }
    }
  }
  colnames(stat)=names(res)
  rownames(stat)=names(res)
  return(stat)
  
}



listTomatrix=function(res){
  # list to a matrix,  row is union genes
  unionPos=NULL
  for ( i in names(res)){
    unionPos=union(unionPos,res[[i]])
  }
  #  mutation matix
  lc=length(res)
  lr=length(unionPos)
  stat=matrix(rep(0,lc*lr),lr,lc)
  for( i in seq(lc)){
    #	index=match(unionPos,tmp$pos)
    #	stat[!is.na(index),i]=1
    index=match(res[[i]],unionPos)
    stat[index,i]=1
    #	stat[index,i]=1
  }
  colnames(stat)=names(res)
  rownames(stat)=unionPos
  return(stat)
}



#  lcl <- con$clusters$multi$groups;
#  x <- jdf(as.factor(lcl),as.factor(anoCell))

  #comparison plot
  # calculate Jaccard comparison data frame with $overlap and $jacard columns
  jdf <- function(f1,f2,n.cores=30) {
    do.call(rbind,mclapply(levels(f1),function(l1) {
      do.call(rbind,lapply(levels(f2),function(l2) {
        n1 <- names(f1)[f1==l1]; n2 <- names(f2)[f2==l2];
        ov <- length(intersect(n1,n2)); tv <- length(union(n1,n2))
        data.frame(l1=l1,l2=l2,overlap=ov,jaccard=ov/tv)
      }))
    },mc.cores=n.cores,mc.preschedule=TRUE))
  }
  





rawToPagoda=function(datraw,appname,tSNE=TRUE){

	genelists <- lapply(datraw, function(x) rownames(x))
	str(genelists)
	commongenes <- Reduce(intersect,genelists)
	length(commongenes)


	matrices_raw <- mapply(function(mm, name) {
	 mm[commongenes,]
	},
	datraw,
	names(datraw))


bigM2 <- Reduce(cbind, matrices_raw)


p2=runPagoda(bigM2,appname,n.cores = 12,get.tsne=tSNE)

  f1=paste(appname,'_p2combined.rds',sep='')
  saveRDS(p2,f1)


return(p2)
}



#p2_normal=rawToPagoda(raw.mats[Normal],'all_normal')




runConos=function(datlp2,appname,n.cores=8){
  print(names(datlp2))
  con <- Conos$new(datlp2,n.cores=n.cores)
  con$buildGraph()
  con$findCommunities()


  f1=paste(appname,'_conos_clustering.pdf',sep='')
  p1=con$plotGraph()
  ggsave(f1,p1)

  f1=paste(appname,'_conos.rds',sep='')
  saveRDS(con,f1)

con$embedGraph(method = 'UMAP')

p=con$plotGraph()
ggsave(paste(appname,'.umap.png',sep=''),p,height=7,width=7)


saveRDS(con$embedding,paste(appname,'.umap.rds',sep=''))


con$embedGraph(method="UMAP", n.cores=30,min.dist=1e-10,n.neighbors=50)

p2a <- con$plotGraph(alpha=0.1,size=0.1,plot.na=T)

saveRDS(con$embedding,paste(appname,'.umap2.rds',sep=''))


ggsave(paste(appname,'.umap409.png',sep=''),p2a,height=7,width=7)


  return(con)


}
     






#

#f1='/home/meisl/Workplace/BMME/a.data/conos/Whole_all/Whole_all_conos.rds'
#tracking_conos(f1,'Whole_all',conosCluster,i,appname,jcl3.coarse)

tracking_conos=function(fin,fname,conosCluster,i,appname,jcl3.coarse){
  
  #  tcon=readRDS('/home/meisl/Workplace/BMME/a.data/conos/Whole_all/Whole_all_conos.rds')
  #  i=3
  print(fin)
  tcon=readRDS(fin)
  sle=conosCluster[conosCluster %in% as.character(i)]
  jcl3.coarse_cell=as.character(jcl3.coarse)
  names(jcl3.coarse_cell)=names(jcl3.coarse)
  jcl3.coarse_cell[!is.na(jcl3.coarse_cell)]='other'
  jcl3.coarse_cell[names(sle)]=paste('group',i,sep='_')
  
  a1=tcon$plotGraph(groups =jcl3.coarse_cell,show.legend=TRUE,title=appname)
  a2=tcon$plotGraph(groups =jcl3.coarse,show.legend=TRUE,title=appname)
  
  b=  cowplot::plot_grid(plotlist=list(a1,a2), ncol=2, nrow=1)
  
  pdf1=paste(appname,'_',i,'_',fname,'.conos.png',sep='') 
  ggsave(pdf1,b,width = 13,height=6)
  
}



   
 drawfigureP2_muti=function(p2,appname,jcl3.coarse=NULL,cell_ano_sampleType=NULL,cell_ano_sample=NULL,alpha=0.1){
   
   if (!is.null(cell_ano_sampleType)){
     lrow=2
   }else{
     lrow=1
   }   
   
   lcol=2
   filename=paste(appname,'.p2.tSNE.png',sep='')
   mar = c(0.5, 0.5, 0.5, 0.5)
   mgp = c(2, 0.65, 0)
   panel.size=500
   cex=1.2
   mark.cluster.cex = 1
   mark.clusters = TRUE
  
   png(file = filename, height = lrow * panel.size, 
       width = lcol * panel.size)
   par(mfrow = c(lrow, lcol), mar = mar, mgp = mgp, 
       cex = cex)
   
   p2$plotEmbedding(groups=p2$clusters$PCA$multilevel,type='PCA',embeddingType='tSNE',show.legend=F,mark.clusters=T,min.group.size=1,shuffle.colors=F,mark.cluster.cex=1,alpha=alpha,main='mutil level')

   
   if (!is.null(jcl3.coarse)){
     p2$plotEmbedding(type='PCA',groups =cell_ano_sampleType,embeddingType='tSNE',show.legend=F,mark.clusters=T,min.group.size=1,shuffle.colors=F,mark.cluster.cex=1,alpha=alpha,main='Cell annotation')
   }
   
   
   if (!is.null(cell_ano_sampleType)){ 
     p2$plotEmbedding(type='PCA',groups =cell_ano_sampleType,embeddingType='tSNE',show.legend=F,mark.clusters=T,min.group.size=1,shuffle.colors=F,mark.cluster.cex=1,alpha=alpha,main='Sample Type')
   }
   
   
   if (!is.null(cell_ano_sample)){
     p2$plotEmbedding(type='PCA',groups =cell_ano_sample,embeddingType='tSNE',show.legend=F,mark.clusters=T,min.group.size=1,shuffle.colors=F,mark.cluster.cex=1,alpha=alpha,main='Sample')
   }
   
   dev.off()
 }
 
 
# drawfigureP2_muti(p2,appname,jcl3.coarse=anoCell,cell_ano_sampleType=anoCell,cell_ano_sample=NULL,alpha=0.1)
   




#  run DAVID  source actviate py34 
runDAVID<-function(genelist,appname){
  genelist=genelist[!is.na(genelist)]
  ENTREZID=unlist(mget(genelist, org.Hs.egSYMBOL2EG, ifnotfound=NA))
  ENTREZID=ENTREZID[!is.na(ENTREZID)]
  
  tmpF=paste(appname,'.ENTREZID',sep='')
  fout=paste(appname,'.DAVID.txt',sep='')
  write.table(ENTREZID,tmpF,sep='\t',row.names=F,col.names=F,quote=F)
  cmd=paste("python /d0/home/meisl/Workplace/Archive/DAVID/DAVID_annotation_chart.py --genelist=", tmpF ,' --email="ma.tongji@gmail.com" --out=',fout,sep='')  
  system(cmd)
  
  cmd2=paste('Rscript /d0/home/meisl/Workplace/Archive/DAVID/DAVID_to_figure_vison.r ',fout,' --termNumber=18  -o ',appname,sep='')
  system(cmd2)
  
  print(cmd2)
  ENTREZID2=names(ENTREZID)
  names(ENTREZID2)=ENTREZID
  
  tmp=read.csv(fout,sep='\t',header=T)
  tmp$symb=apply(tmp,1,function(x) paste(ENTREZID2[strsplit(as.character(unlist(x['Genes'][1])),', ')[[1]]],collapse=','))
  
  write.table(tmp,fout,sep='\t',row.names=F,col.names=T,quote=F)
  return(tmp)
}




# r1=ExtractCell(conA,nname)   nname is selected cell name

#  conos subset graph
ExtractCell=function(conA,cname){
  g=conA$graph
  g2=igraph::subgraph(g,cname)
  method <- igraph::multilevel.community
  km <- method(g2)
  mutil <- km$membership
  names(mutil) <- km$names
  
  method <- igraph::label.propagation.community
  km <- method(g2)
  label <- km$membership
  names(label) <- km$names
  
  method <- igraph::walktrap.community
  km <- method(g2)
  walktrap <- km$membership
  names(walktrap) <- km$names
  
  tmp=leiden.community(g2)
  leiden <- tmp$membership
  
  tmp=leiden.community(g2,resolution = 2, n.iterations = 8)
  leiden2 <- tmp$membership
  print(table(leiden2))
  wij2 <- as_adj(g2, attr = "weight")
  coords2 <- conos:::projectKNNs(wij = wij2, dim = 2, verbose = TRUE, 
                                 sgd_batches = 1e+08, gamma = 1, M = 1, seed = 1, 
                                 alpha = 0.2, rho = 1, threads = 8)
  colnames(coords2) <- colnames(wij2)
  res=list('graph'=g2,'mutil'=mutil,'label'=label,'walktrap'=walktrap,'largvis'=t(coords2))
  
} 





larvisTotSNE=function(g2,ndim=20){
  coords20 <- conos:::projectKNNs(wij = wij2, dim = ndim, verbose = TRUE, 
                                  sgd_batches = 1e+08, gamma = 1, M = 1, seed = 1, 
                                  alpha = 0.2, rho = 1, threads = 8)
  colnames(coords20) <- colnames(wij2)
  
  emb <-Rtsne.multicore::Rtsne.multicore(t(coords20),  perplexity=50, num_threads=10)$Y
  rownames(emb) <- colnames(coords20)
  return(emb)
}











factorEmbedding=function(emb,anoSample,appname,alpha=0.1,lcol=5,panel.size=400,mark.cluster.cex = 1){
  lcol=lcol
  lrow=ceiling(length(unique(anoSample))/lcol)
  filename=paste(appname,'.indivisual.emb.png',sep='')
  mar = c(0.5, 0.5, 0.5, 0.5)
  mgp = c(2, 0.65, 0)
  panel.size=panel.size
  cex=1.2
  mark.cluster.cex = mark.cluster.cex
  mark.clusters = TRUE
  
  png(file = filename, height = lrow * panel.size, 
      width = lcol * panel.size)
  par(mfrow = c(lrow, lcol), mar = mar, mgp = mgp, 
      cex = cex)
  
  
  for (i in unique(anoSample)){
    
    y2=anoSample
    y2[y2!=i]='other'
    plotEmbedding(emb,groups=y2, mark.clusters=TRUE, alpha=alpha,main=i)
    
  }
  
  dev.off()

}






mutiEmbedding=function(emb,gps,appname,alpha=0.1,lcol=2,panel.size=500,mark.cluster.cex = 1){
  lcol=lcol
  lrow=ceiling(length(unique(gps))/lcol)
  filename=paste(appname,'.emb.png',sep='')
  mar = c(0.5, 0.5, 0.5, 0.5)
  mgp = c(2, 0.65, 0)
  panel.size=panel.size
  cex=1.2
  mark.cluster.cex = mark.cluster.cex
  mark.clusters = TRUE
  
  png(file = filename, height = lrow * panel.size, 
      width = lcol * panel.size)
  par(mfrow = c(lrow, lcol), mar = mar, mgp = mgp, 
      cex = cex)
  
  for (i in names(gps)){
    
    y2=gps[[i]][rownames(emb)]
    plotEmbedding(emb,groups=y2, mark.clusters=TRUE, alpha=alpha,main=i,mark.cluster.cex=mark.cluster.cex,show.legend=TRUE)
    
  }
  
  dev.off()
  
}
# factorEmbedding(emb,anoSample,paste(appname,'.sample',sep=''),alpha=0.05)
#factorEmbedding(emb,anoCell,paste(appname,'.cell',sep=''),alpha=0.05)

#







basicP2proc2=function (cd, n.cores = 1, batch = NULL, n.odgenes = 3000, nPcs = 100, 
                       k = 30, perplexity = 50, log.scale = TRUE, trim = 10, keep.genes = NULL, 
                       min.cells.per.gene = 0, min.transcripts.per.cell = 100, get.largevis = TRUE, 
                       get.tsne = TRUE, make.geneknn = TRUE,alphaf=0.05,gs=NULL,rm=NULL) 
{
  rownames(cd) <- make.unique(rownames(cd))
  p2 <- Pagoda2$new(cd, n.cores = n.cores, batch = batch, keep.genes = keep.genes, 
                    trim = trim, log.scale = log.scale, min.cells.per.gene = min.cells.per.gene, 
                    min.transcripts.per.cell = min.transcripts.per.cell)
  p2$adjustVariance(plot = F, gam.k = 10,alpha = alphaf)
  
  p2$misc$odgenes.old=p2$misc$odgenes
  if (!is.null(gs)){
    p2$misc$odgenes=intersect(p2$misc$odgenes,gs)
  }
  
  if (!is.null(rm)){
    p2$misc$odgenes=setdiff(p2$misc$odgenes,rm)
  }
  
  print(length(p2$misc$odgenes))
  p2$calculatePcaReduction(nPcs = nPcs, n.odgenes = n.odgenes, 
                           maxit = 1000)
  p2$makeKnnGraph(k = k, type = "PCA", center = TRUE, weight.type = "none", 
                  n.cores = n.cores, distance = "cosine")
  p2$getKnnClusters(method = igraph::multilevel.community, 
                    type = "PCA", name = "multilevel")
  if (get.largevis) {
    M <- 30
    p2$getEmbedding(type = "PCA", embeddingType = "largeVis", 
                    M = M, perplexity = perplexity, gamma = 1/M, alpha = 1)
  }
  if (get.tsne) {
    if (perplexity > nrow(p2$counts)/5) {
      perplexity <- floor((nrow(p2$counts) - 1)/3)
      cat("perplexity is too large, reducing to", perplexity, 
          "\n")
    }
    p2$getEmbedding(type = "PCA", embeddingType = "tSNE", 
                    perplexity = perplexity, distance = "L2")
  }
  if (make.geneknn) {
    p2$makeGeneKnnGraph()
  }
  invisible(p2)
}




DEheatmap=function(p2myeloid,conosCluster,appname,removeGene=NULL,num=20){
  de1 <- p2myeloid$getDifferentialGenes(groups=conosCluster)
  ggg=NULL
  for(n in names(de1)){
    x=de1[[n]]
    z <- x[order(-x$Z),]
    if(!is.null(removeGene)){
      index=grepl('^RP[LKS]',rownames(z))
      z=z[!index,]
      index=grepl('^MT-',rownames(z))
      z=z[!index,]
      index=grepl('^HSP',rownames(z))
      z=z[!index,]
    }
    markers=rownames(z)[1:num]
    #markers=markers[!is.na(markers)]
    
    #allg=cbind(allg,markers)
    ggg=c(ggg,markers)
    
  }
  
  markers=unique(ggg)
  
  
  
  
  
  x <- as.matrix(p2myeloid$counts[names(conosCluster),markers])
  ## trim outliers
  x <- apply(x, 2, function(xp) {
    qs <- quantile(xp,c(0.01,0.98))
    xp[xp<qs[1]] <- qs[1]
    xp[xp>qs[2]] <- qs[2]
    xp
  })
  x <- x[,(apply(x,2,sd) != 0)]
  x <- t(scale(x))
  ## sample 2000 cells for plotting
  #x <- x[,sample(ncol(x),2000)]
  o <- order(conosCluster[colnames(x)])
  
  
  x=x[,o]
  ss=apply(x,2,function(x) sum(x))
  
  x=x[,abs(ss)>0.1]
  
  
  annot <- data.frame(CellType=conosCluster[colnames(x)],row.names = colnames(x))
  
  
  annot2=annot
  
  #annot2 = data.frame(ID = as.factor(as.character(annot[,1])))
  rownames(annot2)=colnames(x)
  
  
  cellA=annot2[,1]
  names(cellA)=rownames(annot2)
  
  o <- order(cellA)

  x=x[,o]
  
  pal <- colorRampPalette(c('navy','white','firebrick3'))(50)
  ## draw heatmap
  
  fout=paste(appname,'.DiffG.png')
  rgb.palette <- colorRampPalette(c("blue","white","red"), space = "rgb" )
  heat=pheatmap(x,cluster_cols=FALSE,annotation_col = annot2,show_colnames = F,annotation_legend = TRUE, #,gaps_col =  1,
                cluster_rows = FALSE,color=rgb.palette(100),filename=fout,fontsize_row =5,width=8,height=8.8,   #3.5*0.02*length(markers),
                breaks = c(seq(min(x),-0.01,length.out = 50),0.01,seq(0.1,2,length.out = 48),max(x)))
 return(heat)
}












DEcaculate2.term=function(p2,appname,conosCluster,removeGene=NULL,cutoff=3,num=100,GO=NULL){
  allg=NULL  
  de1 <- p2$getDifferentialGenes(groups=conosCluster,z.threshold = cutoff)
  
  folderN=paste(appname,'.DE',sep='')
  pwd=getwd()
  pwd2=paste(pwd,'/',folderN,'/',sep='')
  
  system(paste('mkdir ',folderN))
  
  setwd(pwd2)
  
  
  
  
  f1=paste(appname,'_diffGene.rds',sep='')
  saveRDS(de1,f1)
  for(n in names(de1)){
    x=de1[[n]]
    z <- x[order(-x$Z),]
    if(!is.null(removeGene)){
      index=grepl('^RP[LKS]',rownames(z))
      z=z[!index,]
    }
    markers=rownames(z)[1:num]
    markers=markers[!is.na(markers)]
    allg=cbind(allg,markers)
    x <- as.matrix(p2$counts[names(conosCluster),markers])
    ## trim outliers
    x <- apply(x, 2, function(xp) {
      qs <- quantile(xp,c(0.01,0.98))
      xp[xp<qs[1]] <- qs[1]
      xp[xp>qs[2]] <- qs[2]
      xp
    })
    x <- x[,(apply(x,2,sd) != 0)]
    x <- t(scale(x))
    ## sample 2000 cells for plotting
    #x <- x[,sample(ncol(x),2000)]
    o <- order(conosCluster[colnames(x)])
    annot <- data.frame(cluster=conosCluster[colnames(x)],row.names = colnames(x))
    
    annot$dtype='other'
    annot[as.character(annot[,1])==n,'dtype']=n
    annot$dtype=as.factor(annot$dtype)
    
    
    pal <- colorRampPalette(c('navy','white','firebrick3'))(50)
    ## draw heatmap
    
    fout=paste(appname,'_',n,'_marker.heatmap.png',sep='')
    rgb.palette <- colorRampPalette(c("blue","white","red"), space = "rgb" )
    pheatmap(x[,o],labels_col=FALSE,cluster_cols=FALSE,annotation_col = annot,show_colnames = F,
             cluster_rows = FALSE,color=rgb.palette(100),filename=fout,fontsize_row =4,width=5,height=8,   #3.5*0.02*length(markers),
             breaks = c(seq(min(x),-0.01,length.out = 50),seq(0.01,max(x),length.out = 50)))
    
    
    iid= paste(appname,'_cluster_',n,sep='')
    if(!is.null(GO)){
      
      figGO=GOanalysis.term(markers,iid) 
      figGO=figGO[['BP']]
      
      for (jj in seq(nrow(figGO))){
        
        gs=figGO[jj,'genes']
        gs=strsplit(gs,',')[[1]]
        
        if (length(gs) >= 4 ){
        
        fout=paste(appname,'_',n,'_',figGO[jj,'Term'],'_BP.heatmap.png',sep='')
        rgb.palette <- colorRampPalette(c("blue","white","red"), space = "rgb" )
        pheatmap(x[gs,o],labels_col=FALSE,cluster_cols=FALSE,annotation_col = annot,show_colnames = F,
                 cluster_rows = FALSE,color=rgb.palette(100),filename=fout,fontsize_row =4,width=5,height=5*0.02*length(gs),
                 breaks = c(seq(min(x),-0.01,length.out = 50),seq(0.01,max(x),length.out = 50)))
        
      }}
      
      
      }
  }
  colnames(allg)=names(de1)
  
  
  write.table(allg,paste(appname,'_Diff_gene.xls',sep=''),sep='\t',col.names=T,row.names=F,quote=F)
  setwd(pwd)
}  






GOanalysis.term=function(markers,n){
  
  ENTREZID=unlist(mget(markers, org.Hs.egSYMBOL2EG, ifnotfound=NA))
  ENTREZID=ENTREZID[!is.na(ENTREZID)]
  
  allr=list()
  
  for(function_type in c("BP", "CC", "MF")){
    
    param <- new("GOHyperGParams", geneIds=ENTREZID,
                 #universe=universe,
                 annotation="org.Hs.eg.db", ontology=function_type,pvalueCutoff=0.1,
                 conditional=FALSE, testDirection="over")
    hyp <- hyperGTest(param)
    sumTable <- summary(hyp)
    
    
    
    david=sumTable[1:20,]
    david$Pvalue=-log(david[,2])
    termNumber=nrow(david)
    
    david$genes=apply(david,1,function(x) { paste(names(ENTREZID[ENTREZID %in% get(as.character(x[1]),org.Hs.egGO2ALLEGS)]),collapse = ',' ) } )
    
    allr[[function_type]]=david
    write.table(david,paste(n,'_' ,function_type, ".xls", sep=""),sep='\t',col.names=T,row.names=F,quote=F)
    
    
    library(ggplot2)
    p1 <- ggplot(data=david, aes(x=Pvalue, y=Term, size=Count, colour = factor(david$Count)))
    p1 <- p1 + geom_point()
    p1 <- p1 + guides(color = FALSE)
    p1 <- p1 + theme(panel.grid.major = element_line(colour='blue'),
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank())
    p1 <- p1 + xlab(paste("-log10(Pvalue)", sep="")) + ylab("")
    p1 <- p1 + labs(title=paste("DAVID:", function_type, sep=""))
    p1 <- p1 + theme(axis.text.x=element_text(size=10, face="plain", colour ='black'))
    p1 <- p1 + theme(axis.text.y=element_text(size=6, face="plain", colour ='black'))
    #p1=p1+theme(axis.title.y=element_text(size=10,face="plain",colour ='black'))
    #p1=p1+theme(legend.position="bottom")
    p1 <- p1 + xlim(min(david$Pvalue), max(david$Pvalue))
    #p1=p1+ylim(-10,+15)
    print(p1)
    ggsave(file=paste(n,'_' ,function_type, ".png", sep=""), scale=0.8, dpi=600, width = 7, height=1+0.25*termNumber) 
  } 
  
  saveRDS(allr,paste(n,'.GOterm.rds',sep=''))
  return(allr)
}






prepare_rawList2=function(bigM2,group4,atype,anoSample,ncut=15){
  cell_ano=group4[group4==atype]
  cname=names(cell_ano)
  
  
  #index=grepl('BMET1-',names(anoSample))
  #anoSample=anoSample[!index] 
  gname=names(anoSample)
  
  t.bigM2=bigM2[,gname]
  
  tab=table(anoSample[gname])
  nnames=names(tab[tab>ncut])
  
  raw.mats=list()
  
  n.anoSample=anoSample[gname]
  p2.lis=list()
  for (n in nnames){
    cd=t.bigM2[,n.anoSample==n]
    n1=n
    raw.mats[[n1]]=cd
  }
  nn=unlist(lapply(raw.mats, function(x) ncol(x)))
  print(atype)
  print(nn)
  return(raw.mats)
  
}




GETnomrlizedCount2=function(raw.mats,appname){
  
  count1=do.call(cbind,lapply(raw.mats,function(x) rowSums(as.matrix(x))))
  
  
  cm.tmp=count1
  cm.tmp <- as.matrix(cm.tmp)
  rownames(cm.tmp) <- rownames(cm)
  ## calculate cpm
  cpm <- sweep(cm.tmp, 2, apply(cm.tmp,2, sum), FUN='/')
  dat_cpm <- log10(cpm * 1e6 + 1)
  rownames(dat_cpm)=rownames(count1)
  
  r=list('dat_cpm'=dat_cpm,'raw'=raw.mats)
  
  fo=paste(appname,'.normlized.count.rds',sep='')
  saveRDS(r,fo)
  
  return(dat_cpm)
}


#  function of ligand and receptor correlation 
# @ list of row counts 
# cell annotation 

#r=runLigand_Receptor(bigM,'Macrophage1','Macrophage2',p2,anoCell)




fsubset=function(raw,anoCell,input_cell,cutoff=30){

  cname=names(anoCell[anoCell %in% input_cell])
  
  raw2=lapply(raw,function(x) x[,intersect(cname,colnames(x))])
  
  num=unlist(lapply(raw2,function(x) ncol(x)))
  
  
  print(num)
  
  num=num[ num>cutoff ]
  raw2=raw2[names(num) ]
  
  print(unlist(lapply(raw2,function(x) ncol(x))))
  
  return(raw2)
}







Toch=function(conosCluster){
  cluster2=as.character(conosCluster)
  names(cluster2)=names(conosCluster)
  return(cluster2)
}




prepare_rawList3=function(bigM,group4,atype,anoSample,ncut=15){
  gname=names(group4[group4==atype])
  print(length(gname))
  
  tab=table(anoSample[gname])
  nnames=names(tab[tab>15])
  
  raw.mats=list()
  
  n.anoSample=anoSample[gname]
  for (n in nnames){
    cd=bigM[,names(n.anoSample[n.anoSample==n])]
    n1=n
    raw.mats[[n1]]=cd
  }
  
  nn=unlist(lapply(raw.mats, function(x) ncol(x)))
  print(atype)
  print(nn)
  return(raw.mats)
  
  
}



prepare_rawList3=function(bigM,group4,atype,anoSample,ncut=15){
  gname=names(group4[group4==atype])
  print(length(gname))
  
  tab=table(anoSample[gname])
  nnames=names(tab[tab>15])
  
  raw.mats=list()
  
  n.anoSample=anoSample[gname]
  for (n in nnames){
    cd=bigM[,names(n.anoSample[n.anoSample==n])]
    n1=n
    raw.mats[[n1]]=cd
  }
  
  nn=unlist(lapply(raw.mats, function(x) ncol(x)))
  print(atype)
  print(nn)
  return(raw.mats)
  
  
}




#r=runLigand_Receptor(bigM,'Macrophage1','Macrophage2',p2,anoCell)
runLigand_Receptor2=function(bigM,c1,c2,p2,anoCell,anoSample){
  
  
  cell1=prepare_rawList3(bigM,anoCell,c1,anoSample,ncut=15)
  cell2=prepare_rawList3(bigM,anoCell,c2,anoSample,ncut=15)
  
  
  dcell1=GETnomrlizedCount2(cell1,'cell1')
  dcell2=GETnomrlizedCount2(cell2,'cell2')
  
  inter=intersect(colnames(dcell1),colnames(dcell2))
  
  dcell1=dcell1[,inter]
  dcell2=dcell2[,inter]
  table(colnames(dcell1)==colnames(dcell2))
  
  
  RL=readRDS('/d0/home/meisl/bin/data/Ligand_Receptor/Ligand_recptor.rds')
  
  index1= as.character(RL[,1]) %in% rownames(dcell1)
  RL=RL[index1,]
  
  index1= as.character(RL[,2]) %in% rownames(dcell1)
  RL=RL[index1,]
  
  selected=unique(c(as.character(RL[,1]),as.character(RL[,2])))
  
  
  ## 
  de1 <- p2$getDifferentialGenes(groups=anoCell[anoCell %in% c(c1,c2)],z.threshold = -100)
  names(de1)
  pairsL=RL[as.character(RL[,1]) %in% selected , ]
  pairsL$cell1=de1[[c1]][as.character(pairsL[,1]),'Z']
  pairsL$cell2=de1[[c2 ]][as.character(pairsL[,2]),'Z']
  Ligand=apply(pairsL,1,function(x) cor(dcell1[as.character(x[1]),],dcell2[as.character(x[2]),]))
  pairsL$correlation=Ligand
  #pairsL=pairsL[order(abs(pairsL$correlation),decreasing=T),]
  
  colnames(pairsL)[3]=paste('Ligand_',c1,sep='')
  colnames(pairsL)[4]=paste('Receptor_',c2,sep='')
  
  
  
  pairsR=RL[as.character(RL[,2]) %in% selected , ]
  Receptor=apply(pairsR,1,function(x) cor(dcell1[as.character(x[2]),],dcell2[as.character(x[1]),]))
  pairsR$cell1=de1[[c1 ]][as.character(pairsR[,2]),'Z']
  pairsR$cell2=de1[[c2 ]][as.character(pairsR[,1]),'Z']
  pairsR$correlation=Receptor
  #pairsR=pairsR[order(abs(pairsR$correlation),decreasing=T),]
  
  colnames(pairsR)[3]=paste('Receptor_',c1,sep='')
  colnames(pairsR)[4]=paste('Ligand_',c2,sep='')
  
  res=list('pairsL'=pairsL,'pairsR'=pairsR)
  return(res)
}


#r=runLigand_Receptor(bigM,'Macrophage1','Macrophage2',p2,anoCell)


read10xMatrix2=function (path) {
  matrixFile <- paste(path, "/matrix.mtx.gz", sep = "")
  genesFile <- paste(path, "/features.tsv.gz", sep = "")
  barcodesFile <- paste(path, "/barcodes.tsv.gz", sep = "")
  if (!file.exists(matrixFile)) {
    stop("Matrix file does not exist")
  }
  if (!file.exists(genesFile)) {
    stop("Genes file does not exist")
  }
  if (!file.exists(barcodesFile)) {
    stop("Barcodes file does not exist")
  }
  x <- as(Matrix::readMM(gzfile(matrixFile)), "dgCMatrix")
  genes <- read.table(gzfile(genesFile))
  rownames(x) <- genes[, 2]
  barcodes <- read.table(gzfile(barcodesFile))
  colnames(x) <- barcodes[, 1]
  invisible(x)
}

