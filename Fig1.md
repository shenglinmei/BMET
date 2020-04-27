---
title: "R Notebook Fig1"
output:
  rmarkdown::github_document:
    toc: true
---

#load data and cell annotation 

```{r fig.width=5, fig.height=5}

require(tidyr)
require(dplyr)



source('../utils.R')
source('../Lib.r')

scon=readRDS('..//BMET_all.conos.rds')

samples=scon$samples
xf=names(samples)
xf <- gsub("Whole","Benign",xf);
xf <- gsub("Noninvolved","Distal",xf);
names(samples)=xf
scon$samples=samples


load('../data/ano.RData')


# color palette for cell type 
cellcol=readRDS('../data/cellType.color.rds')




```


Joint emebdding 

```{r fig.width=5, fig.height=5}

pf <- function(n) { cellcol[1:n]}
p0 <- scon$plotGraph(raster=T,groups=typefc2,font.size=c(3.4,5),alpha=0.08,plot.na=F,palette=pf,size=0.8)
p0
ggsave('F1C.png',p0,height=5,width=5)


```


emebdding of indivisual sample

```{r fig.width=5, fig.height=5}

typefc=as.factor(typefc2)




get.fractions <- function(x,correct=T) {
  xn <- names(x);
  xf <- gsub(".*-","",as.character(x));
  xf <- gsub("Whole","Benign",xf);
  xf <- gsub("Noninvolved","Distal",xf);
  
  xf <- factor(xf,levels=c('Benign','Distal','Involved','Tumor'))
  if(!is.null(xn)) names(xf) <- xn;
  xf
}
get.patients <- function(x) { gsub("-.*","",x) }

## FAILING



#order 
x <- names(scon$samples)
do.call(cbind,tapply(x,get.fractions(x),get.patients))
tapply(x,get.fractions(x),function(z) sort(get.patients(z)))

#pl <- lapply(scon$samples,function(x) conos:::embeddingPlot(x$embeddings$PCA$tSNE,groups=typefc,mark.groups=F,alpha=0.2,plot.na=F)+theme_void() + theme(legend.position="none",panel.border = element_rect(colour = "black", fill=NA, size=0.2)) )


pl <- lapply(scon$samples,function(x) embeddingPlot(x$embeddings$PCA$tSNE,groups=typefc,palette=pf,mark.groups=F,alpha=0.2,plot.na=F,raster=T,raster.width=1,raster.height=1,size=0.3)+theme_void() + theme(legend.position="none",panel.border = element_rect(colour = "black", fill=NA, size=0.3)) )

pl <- lapply(sn(names(scon$samples)),function(x) embeddingPlot(scon$samples[[x]]$embeddings$PCA$tSNE,groups=typefc,palette=pf,mark.groups=F,alpha=0.2,plot.na=F,raster=T,raster.width=1,raster.height=1,size=0.3)+theme_void() + 
               annotate("text", Inf, -Inf, label = x,vjust = -2,hjust = 1.3,size=4.4 )+
                theme(legend.position="none",panel.border = element_rect(colour = "black", fill=NA, size=0.3)) )



base.emb=scon$embedding
emblist=lapply(scon$samples, function(x) {
  tmp=x$embeddings$PCA$tSNE
  inter=intersect(rownames(tmp),names(typefc))
  base.emb[inter,]
})

#pl <- lapply(emblist,function(x) conos:::embeddingPlot(x,groups=typefc,mark.groups=F,alpha=0.3,plot.na=F,raster=T,raster.width=1,raster.height=1,size=0.05)+theme_void() + theme(legend.position="none",panel.border = element_rect(colour = "black", fill=NA, size=0.3)) )

xlims=c(min(base.emb[,1])+0.3,max(base.emb[,1])-0.3)
ylims=c(min(base.emb[,2])+0.3,max(base.emb[,2])-0.3)

pl <- lapply(sn(names(emblist)),function(x) embeddingPlot(emblist[[x]],groups=typefc,palette=pf,mark.groups=F,alpha=0.3,plot.na=F,raster=T,raster.width=1.5,raster.height=1.5,size=1.2)+theme_void() + theme(legend.position="none",panel.border = element_rect(colour = "black", fill=NA, size=0.3))+ xlim(xlims) + ylim(ylims)+ 
               annotate("text", Inf, Inf, label =strsplit(x,'-')[[1]][1] ,hjust = 1.1,vjust=1.3,size=2) )
#pl[[2]]




pats <- unique(get.patients(x));
patient.palette <- setNames(rainbow(length(pats),v=0.5,s=0.5),sample(pats))

#pl <- lapply(names(scon$samples),function(xn) { x <- scon$samples[[xn]]; conos:::embeddingPlot(x$embeddings$PCA$tSNE,groups=typefc,mark.groups=F,alpha=0.2,plot.na=F)+theme_void() + theme(legend.position="none",panel.border = element_rect(colour = patient.palette[get.patients(xn)], fill=NA, size=2)) } )

# manual ordering ..
pats <- unique(get.patients(x));
ben.pats <- grep("^BMM",pats,value=T)
pats <- c(setdiff(pats,c(ben.pats,'BMET10','BMET7')),c('BMET10','BMET7'))
plo <- unlist(lapply(1:length(pats),function(i) match(c(paste(ben.pats[i],'Benign',sep='-'),paste(pats[i],c('Distal','Involved','Tumor'),sep='-')),x)))

all.panels <- cowplot::plot_grid(plotlist=pl[plo],ncol=4)

pdf(file='all.panels.pdf',width=3,height=3/4*length(pats))
print(all.panels)
dev.off();




library(cowplot)

sample.spacing <- 0.005;
column.spacing <- 0.01;
all.panels <- ggdraw();
for(i in 1:length(plo)) {
  ir <- ceiling(i/4); ic <- (i-1) %% 4;
  if(!is.na(plo[i])) { all.panels <- all.panels+draw_plot(pl[[plo[i]]],ic*1/4+column.spacing/2,1-ir/length(pats)+sample.spacing/2,1/4-column.spacing,1/length(pats)-sample.spacing) }
  
}
#all.panels



```



```{r fig.width=6, fig.height=9}

base.width <- 3;
#pdf(file='all.panels2.pdf',width=base.width,height=base.width*(1/4-column.spacing)/(1/length(pats) - sample.spacing))
print(all.panels)
#dev.off();


```


# fraction plots
# bening vs. non-cancer
```{r fig.width=5, fig.height=5}



x <- setNames(ifelse(fractionf=='Benign','Benign','Cancer'),names(typefc))
cancer.vs.bening.emb <- scon$plotGraph(groups=x,font.size=c(4,6),alpha=0.1,raster=T,plot.na=F,mark.groups=F,palette=function(n) { as.character(rev(c('gray70',fraction.palette['Benign']))) },size=0.5,show.legend=T)+theme(legend.position=c(0.8,0.2),legend.title=element_text(size=0), legend.text=element_text(size=18), legend.key.size = unit(1.7, 'lines'))+ guides(color = guide_legend(override.aes = list(size=2,alpha=1)))
cancer.vs.bening.emb

pdf(file='cancer.vs.bening.pdf',width=4,height=4)
print(cancer.vs.bening.emb)
dev.off();

# tumor vs. involved+non-involved
x <- setNames(fractionf=='Tumor',names(typefc))
x <- setNames(ifelse(fractionf=='Tumor','Tumor','Involved/Non'),names(typefc))
x[fractionf=='Benign'] <- NA;

x[x=='Involved/Non']='Involved/Distal'

tumor.vs.invnon.emb <- scon$plotGraph(groups=x,font.size=c(4,6),alpha=0.1,raster=T,plot.na=F,mark.groups=F,palette=function(n) { as.character(c('gray70',fraction.palette['Tumor'])) },size=0.5,show.legend=T)+theme(legend.position=c(0.72,0.32),legend.title=element_text(size=0), legend.text=element_text(size=18), legend.key.size = unit(1.7, 'lines'))+ guides(color = guide_legend(override.aes = list(size=2,alpha=1)))
tumor.vs.invnon.emb


pdf(file='tumor.vs.invnon.pdf',width=4,height=4)
print(tumor.vs.invnon.emb)
dev.off();


```

# fraction change boxplots


```{r fig.width=7, fig.height=5}


fraction=fractionf

#table(names(samplef)==names(typefc))

# how about fractional changes?
library(dplyr)
xt <- table(typefc,samplef)
xf <- t(t(xt)/colSums(xt))
df <- melt(xf)
df$patient <- gsub("-.*","",df$samplef); df$fraction <- gsub(".*-","",df$samplef);
#df <- df[df$fraction %in% c("Tumor","Noninvolved"),]
#df <- df[df$fraction %in% c("Involved","Noninvolved"),]
#df <- df[df$fraction %in% c("Tumor","Involved"),]
df$samplef <- NULL



# probably should jsut compare between cancer fractions and Bening controls
df <- melt(xf)
df$patient <- gsub("-.*","",df$samplef); df$fraction <- gsub(".*-","",df$samplef);
df$fraction <- get.fractions(df$samplef)
#df <- df[df$fraction %in% c("Tumor","Whole"),]
#df <- df[df$fraction %in% c("Involved","Noninvolved"),]
#df <- df[df$fraction %in% c("Tumor","Involved"),]
df <- df %>% arrange(typefc,fraction)


x <- as_tibble(df) %>% group_by(typefc,fraction) %>% summarise_at(vars(value),funs(mean(.,na.rm=T))) %>% ungroup() %>% spread(fraction,value) %>% mutate(diff=Tumor-Benign) %>% arrange(desc(diff)) %>% as.data.frame()

# sort cell types
df$typefc <- factor(df$typefc,levels=as.character(x$typefc))
#df$fraction <- factor(df$fraction,levels=c("Benign","Noninvolved","Involved","Tumor"))

# full plot (for the supp)
p <- ggplot(na.omit(df),aes(x=typefc,y=value,dodge=fraction,fill=fraction))+geom_boxplot(notch=FALSE,outlier.shape=NA) +scale_fill_manual(values=fraction.palette) +  geom_point(position = position_jitterdodge(jitter.width=0.1),color=adjustcolor(1,alpha=0.3),aes(pch=fraction),size=0.8)+ theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y = element_text(angle = 90, hjust = 0.5))  +xlab("") +ylab("Fraction of total cells")+ theme(legend.position="top")
p

pdf(file='fraction.all.pdf',width=7.5,height=4.5)
print(p)
dev.off();







# reduced plot for the main panel
selt <- c("Macrophage","Monocytes","NK","Helper T","Memory T","Mature B","Immature B")
df <- df[df$typefc %in% selt,]
df$typefc <- factor(as.character(df$typefc),levels=selt)
p <- ggplot(na.omit(df),aes(x=typefc,y=value,dodge=fraction,fill=fraction))+geom_boxplot(notch=FALSE,outlier.shape=NA) +scale_fill_manual(values=fraction.palette) +  geom_point(position = position_jitterdodge(jitter.width=0.1),color=adjustcolor(1,alpha=0.3),aes(pch=fraction))+ theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y = element_text(angle = 90, hjust = 0.5))  +xlab("") +ylab("fraction of total cells")+ theme(legend.position="top")
p

pdf(file='fraction.short.pdf',width=6,height=5)
print(p)
dev.off();


```



```{r fig.width=5, fig.height=5}


cellT=names(typefc[typefc=='Tumor'])
T.samp=samplef[cellT]
tab=table(T.samp)
T.samp=T.samp[T.samp %in% names(tab[tab>5])]  # at lest 5 tumor cells 
tab=table(T.samp)

dat=data.frame('name'=names(tab),'num'=tab)
dat$samp=apply(dat,1,function(x) strsplit(x,'-')[[1]][1])
dat$Fraction=apply(dat,1,function(x) strsplit(x,'-')[[1]][2])


dat$samp=ordered(as.factor(dat$samp),levels=c('BMET5','BMET1','BMET11','BMET2','BMET10','BMET6','BMET3'))

ts.within.fraction.all <- ggplot(na.omit(dat),aes(x=samp,y=num.Freq,fill=Fraction))+
  geom_bar(stat="identity") +scale_fill_manual(values=fraction.palette) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y = element_text(angle = 90, hjust = 0.5))  +xlab("") +ylab("Number of tumor cells")+ theme(legend.position="top")



dat$ratio=dat$num.Freq/table(samplef)[dat$name]

ts.within.fraction.all <- ggplot(na.omit(dat),aes(x=samp,y=ratio,fill=Fraction))+
  geom_bar(stat="identity") +scale_fill_manual(values=fraction.palette) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y = element_text(angle = 90, hjust = 0.5))  +xlab("") +ylab("Fraction of total cells")+ theme(legend.position="top")

ts.within.fraction.all

ggsave('Tumor.barplot.ratio.pdf',ts.within.fraction.all,width=5,height=4)


```


#  changes in transcriptional states 

```{r fig.width=5, fig.height=5}


library(Rtsne)
library(abind)

samplef <- gsub("Whole","Benign",samplef);


typefc2=as.factor(typefc2)
cm <- scon$getClusterCountMatrices(groups=typefc)

table(names(typefc2)==names(samplef))


cct <- table(typefc2,samplef)

ctdm <- lapply(sn(colnames(cm[[1]])),function(ct) {
  print(ct)
  tcm <- do.call(rbind,lapply(cm,function(x) x[,ct]))
  tcm <- t(tcm/pmax(1,rowSums(tcm)))
  tcd <- pagoda2:::jsDist(tcm); dimnames(tcd) <- list(colnames(tcm),colnames(tcm));
  # calculate how many cells there are
  attr(tcd,'cc') <- cct[ct,colnames(tcm)]
  tcd
})







min.cells <- 10
# VP_1: expression distances between samples, plotted for each cell type
pdf(file='cell_type_total_tSNE3.pdf',height=4,width=5.5)
x <- lapply(names(ctdm),function(ct) {
 # print(ct)
  xd <- ctdm[[ct]]
  nc <- attr(ctdm[[ct]],'cc');
  vi <- nc[rownames(xd)]>=min.cells;
  xd <- xd[vi,vi]
  if(ncol(xd)>10){  # at lest 10 samples 
    xde <- Rtsne(xd,is_distance=TRUE, perplexity=min(5,floor(nrow(xd)/4)))$Y;
    df <- data.frame(xde); rownames(df) <- rownames(xd); colnames(df) <- c("x","y");
    df$fraction <- gsub(".*-","",rownames(df))
    df$patient <- gsub("-.*","",rownames(df))
    df$ncells <- nc[rownames(df)]
    p1 <- ggplot(df,aes(x,y,color=patient,shape=fraction,size=log10(ncells))) + geom_point() + theme_bw() + ggtitle(ct) + theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank())
    print(p1)
  }
})
dev.off();





# a cube across all cell types and sample pairs
# weights of individual cell types determined by the minimal number of cells on each side of the pairwise comparison
x <- abind(lapply(ctdm,function(x) {
  nc <- attr(x,'cc');
  #wm <- (outer(nc,nc,FUN='pmin'))
  wm <- sqrt(outer(nc,nc,FUN='pmin'))
  return( x*wm )
}),along=3)
# just the weights (for total sum of weights normalization)
y <- abind(lapply(ctdm,function(x) {
  nc <- attr(x,'cc');
  sqrt(outer(nc,nc,FUN='pmin'))
}),along=3)

# normalize by total weight sums
xd <- apply(x,c(1,2),sum)/apply(y,c(1,2),sum)



xde <- Rtsne(xd,is_distance=TRUE, perplexity=4)$Y
xde <- Rtsne(xd,is_distance=TRUE, perplexity=4,max_iter=1e4)$Y

df <- data.frame(xde); rownames(df) <- rownames(xd); colnames(df) <- c("x","y");



df$fraction <- gsub(".*-","",rownames(df))
df$patient <- gsub("-.*","",rownames(df))
df$ncells <- colSums(cct)[rownames(df)]

df$patient <- gsub("BMM.*","Benign",df$patient)
df$patient <- gsub("^[A-Z]$","Normal",df$patient)

p1 <- ggplot(df,aes(x,y,color=patient,shape=fraction,size=log10(ncells))) + geom_point() + theme_bw();
p1



```


# inter-sample variation 

```{r fig.width=5, fig.height=5}



# distance magnitude comparisons
# first, on combined distance matrix
x <- xd; x[upper.tri(x)] <- NA; diag(x) <- NA;
df2 <- na.omit(melt(x))

df2$patient1 <- gsub("-.*","",df2$Var1)
df2$patient2 <- gsub("-.*","",df2$Var2)
df2$fraction1 <- gsub(".*-","",df2$Var1)
df2$fraction2 <- gsub(".*-","",df2$Var2)

df2$samePatient <- df2$patient1==df2$patient2;
df2$sameFraction <- df2$fraction1==df2$fraction2;

df2$withTumor <- df2$fraction1=='Tumor' | df2$fraction2=='Tumor'
df2$withInvolved <- df2$fraction1=='Involved' | df2$fraction2=='Involved'
df2$withNoninvolved <- df2$fraction1=='Distal' | df2$fraction2=='Distal'



df2$type <- NA
df2$type[df2$sameFraction & df2$fraction1=='Benign'] <- 'Benign'
df2$type[df2$sameFraction & df2$fraction1=='Involved'] <- 'Involved'
df2$type[df2$sameFraction & df2$fraction1=='Distal'] <- 'Distal'
df2$type[df2$sameFraction & df2$fraction1=='Tumor'] <- 'Tumor'
#df2$type[df2$sameFraction & df2$fraction1=='Healthy'] <- 'Healthy'
df2$type <- factor(df2$type,levels=c('Benign','Distal','Involved','Tumor'))
#df2$type <- factor(df2$type,levels=c('Healthy','Bening','Noninvolved','Involved','Tumor'))

#fraction.palette <- c('green','dodgerblue1','orange','indianred1')
#fraction.palette <- c('white','green','dodgerblue1','orange','indianred1')

type.palette <- cellcol
#pie(1:length(type.palette),labels=names(type.palette),col=type.palette)



# VP_2: shows distances of datasets within each category of a factor (in this case, fraction), sorting cell types based on between-factor difference magnitude
p <- ggplot(na.omit(df2[df2$sameFraction,]),aes(x=type,y=value))+geom_boxplot(notch=TRUE,outlier.shape=NA,aes(fill=type))+scale_fill_manual(values=fraction.palette)+geom_jitter(position=position_jitter(0.2),color=adjustcolor('black',alpha=0.2))+ theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y = element_text(angle = 90, hjust = 0.5)) +  guides(fill=FALSE) + xlab('') + ylab('expression distance')
p
pdf(file='within.fraction.expression.distances.pdf',width=2,height=4)
print(p)
dev.off();



```


