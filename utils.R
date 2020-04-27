
sn=function(x) { names(x) <- x; return(x); }





getMarkers=function(){
  gg='INOS|IL12|FCGR1A|FCGR1B|FCGR1C|CD80|CXCR10|IL23|CXCL9|CXCL10|CXCL11|CD86|IL1A|IL1B|IL6|TNF|MHCII|CCL5|IRF5|IRF1|CD40|IDO1|KYNU|CCR7'
  M1=strsplit(gg,split='|', fixed=TRUE)[[1]]

  gg='ARG1|ARG2|IL10|CD32|CD163|CD23|FCER2|CD200R1|MARCO|CSF1R|CD206|IL1RA|IL14R|CCL4|CCL13|CCL20|CCl17|CCL18|CCl22|CCL24|LYVE1|VEGFA|VEGFB|VEGFC|VEGFD|EGF|CTSA|CTSB|CTSC|CTSD|TGFB1|TGFB2|TGFB3|MMP14|MMP19|MMP9|CLEC7A|WNT7B|FASL|TNSF12|TNSF8CD276|VTCN1|MSR1|FN1|IRF4'
  M2=strsplit(gg,split='|', fixed=TRUE)[[1]]

  gs='ATF3,BCL2A1,CCL20,CXCL2,DUSP2,EREG,IL1B,PDE4B,PTGS2,PTX3,TNF,TNFAIP3,IL6,G0S2,S100A8,IL8,IL1RN,TREM1,OSM,NLRP3,IL10,AQP9'
  TIM=strsplit(gs,split=',', fixed=TRUE)[[1]]

  f1='/d0-mendel/home/meisl/bin/data/score/GSE5099_Macropage_VS_Mono0h_downMono.txt'
  f2='/d0-mendel/home/meisl/bin/data/score/GSE5099_Macropage_VS_Mono0h_upMacro.txt'
  f3='/d0-mendel/home/meisl/bin/data/score/GSE8286_Macropage_VS_Mono0h_downMono.txt'
  f4='/d0-mendel/home/meisl/bin/data/score/GSE8286_Macropage_VS_Mono0h_upMacro.txt'

  gs=read.csv(f1,sep='\t',header=F)
  GSE5099_Mono=as.character(gs[,1])

  gs=read.csv(f2,sep='\t',header=F)
  GSE5099_Macro=as.character(gs[,1])

  gs=read.csv(f3,sep='\t',header=F)
  GSE8286_Mono=as.character(gs[,1])

  gs=read.csv(f4,sep='\t',header=F)
  GSE8286_Macro=as.character(gs[,1])


  cytotoxicity=c('GZMA','GZMB','GZMM','GZMK','GZMH','PRF1','CD8A','CD8B')


  gs=read.csv('/home/meisl/bin/data/Treg.activity.txt',sep='\t',header=F)
  TregActivity=as.character(gs[,1])


  ehaust3=read.csv('/d0-mendel/home/meisl/bin/data/exhausted.gene.sig',sep='\t',header=F)
  Exhaustion=as.character(ehaust3[,1])



  dat=list('M1'=M1,
           'M2'=M2,
           'TIM'=TIM,
           'GSE5099_Mono'=GSE5099_Mono,
           'GSE5099_Macro'=GSE5099_Macro,
           'GSE8286_Mono'=GSE8286_Mono,
           'GSE8286_Macro'=GSE8286_Macro,
           'cytotoxicity'=cytotoxicity,
           'Exhaustion'=Exhaustion,
           'TregActivity'=TregActivity
           )
  return(dat)
}

