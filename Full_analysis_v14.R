# Meta Deg manuscript ####
rm(list=ls())
require(limma);require(Biobase);require(BiocGenerics)
wd= setwd('C:\\Users/bryan.linggi/Box Sync/PMED_Data(bryan.linggi@robartsinc.com 2)/Projects/')
source('C:/Users/bryan.linggi/Box Sync/PMED_Data(bryan.linggi@robartsinc.com 2)/Code_General/DEGfunctions_v01.R')
# 
# # include methods like descriptions (markerd with &) 
# & We performed a search in NCBI Gene expressed Omnibus (GEO) for microarray datasets using the search terms[ (((ulcerative colitis) AND "Homo sapiens"[porgn] AND "gse"[Filter])) AND (((("expression profiling by array"[DataSet Type]))) AND "gse"[Filter]) AND ("gse"[Filter])] 
# & which resulted in 85 datasets. We excluded datasets that did not contain normal and uc patients,those only containing onlyediatric patietns,  < 10 total samples, was uninflamed samples only, or  those that were from the same cohort, or those that did not contain data in intesity values (exclude thoset that have been Z transformed), resulting in 14 datasets.
# & Each dataset was dowloaded from GEO using GEOquery to create and ExpressionSet containing the intensity values, phenotype data (metadata), and gene feature information. Intensity values
# & were plotted to visualize the expression range using boxplots. If intensity values were not provided from GEO in log2 tranformed, log2 transformation was performed
# & If range values between samples were not similar, quantile normalization was performed.
# & genes were annotated in R using the annotation databases indicated in the GEO description. All probesets were retained (did not exclude those without gene symbols)

####get DEG datasets, output from DEGout codes#####
gse16=readRDS('GSE16879/Output/GSE16879_coloncontrol_DEG.rds')#redid 23mar2020 to include all probes, (before had removed na symbol, but had not done same in other datasets)
gse22= readRDS('GSE22619/Output/GSE22619_DEG.RDS')
 gse38=readRDS('GSE38713/Output/GSE38713_DEG.rds')
 gse59=readRDS('GSE59071/Output/GSE59071_DEG.RDS')
 gse73=readRDS('GSE73661/Output/GSE73661_DEG.rds')
 gse87=readRDS('GSE87466/Output/GSE87466_DEG.RDS')
 gse9452=readRDS('GSE9452/Output/GSE9452_DEG.RDS')
gse97=readRDS('GSE97012/Output/GSE97012_DEG.rds')
gse53=readRDS('GSE53306/Output/GSE53306_DEGrds')
#gse48=readRDS('GSE48634/Output/GSE48634_DEGrds')#uninflamed
gse47=readRDS('GSE47908/Output/GSE47908_DEGrds')
gse42=readRDS('GSE42911/Output/GSE42911_DEGrds')
#gse36=readRDS('GSE36807/Output/GSE36807_DEGrds') # is uninflamed acc to paper
gse133= readRDS('GSE13367/Output/GSE13367_DEGrds')
gse114=readRDS('GSE114527/Output/GSE114527_DEGrds')# , is agilent, cant get in in , tried a lot, # is illumina looks ok as geo download 23mar2020

# 
# #gse105= readRDS('GSE105074/Output/GSE105074_DEG_rect_rcrds')# is uninflamed, do not use
# 
# 
 gse=list('GSE16879'=gse16,'GSE22619'=gse22,'GSE38713'=gse38,'GSE59071'=gse59, 'GSE73661'=gse73,'GSE87466'=gse87 , 'GSE9452'=gse9452,# 7 only
           'GSE97012'=gse97, GSE53306=gse53, GSE47908=gse47, GSE42911=gse42,
           GSE13367=gse133,GSE114527=gse114)
# # # NOTE: GSE105074, GSE48634 is uninflamed, removed from V9 to get GSE14
# # #
#  saveRDS(gse,'../MetaProjects/MetaDEG_v01/Output/GSE14.rds')
#gse=readRDS('../MetaProjects/MetaDEG_v01/Output/GSE14.rds')# get saved 14 datasets (has toptable, topconfect output, and pdata in each list) (saved previously as list)
gsett=sapply(gse, function(x) x$toptable) # 

## info on datasets#####
pdata= print(sapply(gse,function(x) summary(pData(x$input))))#used either Disease or INDC, use this to create summary table
fdata= print(sapply(gse,function(x) summary(fData(x$input))))#used either Disease or INDC, use this to create summary table
adata= print(sapply(gse,function(x) annotation(x$input)))#used either Disease or INDC, use this to create summary table

# pdata= print(sapply(gse,function(x) summary(x$input@phenoData$INDC)))
# td= print(summary(pData(gse$GSE73661$input))) # this one is comparintg on visit not disease

#2##stats on individual by tt####
# input is gse (list of gse lists)
# output is table of tc deg
#& Differential gene expression was determine using limma with abs(log2(FC) > 1.5 and adjust pvalue (benjamini) < .05
lfc=log2(1.5)  #log2(1) is 0
apv=.05
#
DEGtable=data.frame(matrix(nrow=length(gse)))
for (i in 1:length(gse))   {
  In=gse[[i]]$toptable
  apval=    In[which(In$adj.P.Val< apv),]
  up=    apval[which(apval$logFC > lfc ),]
  down= apval [which(apval$logFC < -lfc ),]

  rownames(DEGtable)[i]=names(gse[i])
  colnames(DEGtable)[1]='Up'
  ifelse( nrow(up)   >0,DEGtable$Up[i] <-   nrow(up),   DEGtable$Up[i]   <- 0)
  ifelse( nrow(down) >0,DEGtable$Down[i] <- nrow(down), DEGtable$Down[i] <- 0)
}
#DEGtable$Sum=DEGtable$Up+DEGtable$Down
  print(DEGtable) ## fig 1
# #write.csv(DEGtable,'../MetaProjects/MetaDEG_v01/Output/DEGtable.csv')
# 

#3 Volcano plots####
  #& volcano plots were produced with the R package 'EnhancedVolcano'
require('EnhancedVolcano')

for (i in 1:length(gse)){
dev.new()
res1=gse[[i]]$toptable
symix=grep('symbol',tolower(colnames(res1)))
pl=EnhancedVolcano(res1,
                lab = res1[,symix],
                x = 'logFC',
                y = 'P.Value',
                xlim = c(-5, 8), title = names(gse)[[i]],subtitle = '')
print(pl)
}
 # # 6 Enrichment on individ, using toptable ####
  #& Gene set enrichment analysis was performed with the Reactome Pathways database (reactome.org) using the the R packages
  #... ReactomePA and clusterProfiler using genes that have abs(log2(FC) > 1.5 and adjusted pvalue (benjamini) < .05
 #... the top 10 (lowest adjusted pvalue) pathways for each datasets (many were shared between datasets) were graphed  rm(res1)
   require(ReactomePA);require(clusterProfiler)
  # options************************************************
   lfc=log2(1.5)
   apv= .05
  #  *******************************************************

   enrupTt=vector(mode = "list", length = length(gse))
   enrdownTt=vector(mode = "list", length = length(gse))
   k=1
    require(dplyr)
   for (k in 1:length(gse)){
     res1=gse[[k]]$toptable
     symix=grep('symbol',tolower(colnames(res1)))
     apvix= grep('adj',tolower(colnames(res1)))
     lfcix= grep('logfc',tolower(colnames(res1)))

     ixfcup= res1[,lfcix] > lfc
     ixfcdown= res1[,lfcix] < - lfc
     ixp= res1[,apvix] <   apv
    plot.new()
       if (   sum(ixfcup & ixp)    >0){
     enrupTt[[k]]= enrichuniv(res1[ ixfcup&ixp,symix],univ=res1[,symix]   )
     tp=enrupTt[[k]][1:10,]
     require(ggplot2)
    # p=   ggplot(tp, aes(x=reorder(Description,p.adjust), y=-log10(p.adjust)))+geom_col()+
    #      theme(axis.text.x = element_text(face="bold", color="993333",
    #                                       size=6, angle=45,vjust = 1, hjust=1))+ggtitle(paste('Up;', names(gse)[k]))
    #   plot(p)

     enrdownTt[[k]]=enrichuniv(res1[ ixfcdown &ixp,symix],univ=res1[,symix]   )
    tp=  enrdownTt[[k]][1:10,]
   #     q=ggplot(tp, aes(x=reorder(Description,p.adjust), y=-log10(p.adjust)))+geom_col()+
   #     theme(axis.text.x = element_text(face="bold", color="993333",
   #                                      size=6, angle=45,vjust = 1, hjust=1))+ggtitle(paste('Down;', names(gse)[k]))
   # plot(q)
       }}

    #plot together


   require(plyr)
  a=sapply(enrupTt, length)
   ix=which(a!=0)# some with 0, need to skip
   enrupTt2=enrupTt[ix]#*** Numbers%@@$
   names(enrupTt2)=names(gse)[ix]#*** Numbers%@@$
   allenrTtup = Reduce(function(x,y) merge(x,y,all.x=T,all.y=T,by='Description'),enrupTt2)
   allenrupTtsmall=allenrTtup[,c(1,grep('p.adjust',tolower(colnames(allenrTtup))))]
   colnames(allenrupTtsmall)[2:ncol(allenrupTtsmall)]= names(gse)[ix]#*** Numbers%@@$

   a=lapply(enrupTt2,function(x) x$Description[1:10])
   b=unique(unlist((a)))
   c=allenrupTtsmall[allenrupTtsmall$Description %in% b,]
   require(reshape2)
   mc=melt(c)
   ggplot(mc, aes(x=Description,y=-log10(value), fill=(variable))) +geom_bar(position = 'dodge',stat='identity')+
     theme(axis.text.x = element_text(
         color='black',size=6, angle=45,vjust = 1, hjust=1))+ggtitle('Up')

   ggplot(mc, aes(x=variable,y=-log10(value), fill=(Description))) +geom_bar(position = 'dodge',stat='identity')+
     theme(axis.text.x = element_text(
       color="black",size=6, angle=45,vjust = 1, hjust=1))+ggtitle('Up')
  #down
   a=sapply(enrdownTt, length)
   ix=which(a!=0)# some with 0, need to skip
   enrTt2=enrdownTt[ix]#*** Numbers%@@$
   names(enrTt2)=names(gse)[ix]#*** Numbers%@@$
   allenrTt = Reduce(function(x,y) merge(x,y,all.x=T,all.y=T,by='Description'),enrTt2)
   allenrTtsmall=allenrTt[,c(1,grep('p.adjust',tolower(colnames(allenrTt))))]
   colnames(allenrTtsmall)[2:ncol(allenrTtsmall)]= names(gse)[ix]#*** Numbers%@@$

   a=lapply(enrdownTt,function(x) x$Description[1:10])
   b=unique(unlist((a)))
   c=allenrTtsmall[allenrTtsmall$Description %in% b,]
   require(reshape2)
   mc=melt(c)
   ggplot(mc, aes(x=Description,y=-log10(value), fill=(variable))) +geom_bar(position = 'dodge',stat='identity')+
     theme(axis.text.x = element_text(
       color='black',size=6, angle=45,vjust = 1, hjust=1))+ggtitle('Down')

   ggplot(mc, aes(x=variable,y=-log10(value), fill=(Description))) +geom_bar(position = 'dodge',stat='identity')+
     theme(axis.text.x = element_text(
       color="black",size=6, angle=45,vjust = 1, hjust=1))+ggtitle('Down')

#7 Run MetMetaDegnalysis ######

#### meta ###
  ##run###
#check colnames
 sapply(gse,function(x) colnames(x$toptable))
for (j in 1:length(gsett)){
colnames(gsett[[j]])[grep('symbol',tolower(colnames(gsett[[j]])))]='SYMBOL'
}
sapply(gsett,function(x) colnames(x))

# & DEG metanalysis was perfomred using the R. Package 'MetavolcanoR' using the paramaters cvar=T, collaps=T (yes, is mispelled in code) , and metathr =0.1
 pcriteria="P.Value"
 foldchangecol='logFC'
 genenamecol1="SYMBOL"
 geneidcol=NULL
 collaps=T
 names(gsett)
 require(MetaVolcanoR)
 require(metafor)
 # meta_degs_rem <- rem_mv(diffexp=gsett,
 #                        pcriteria=pcriteria,
 #                         foldchangecol=foldchangecol,
 #                         genenamecol=genenamecol1,
 #                         geneidcol=NULL,
 #                         collaps=T,
 #                         llcol='CI.L',
 #                         rlcol='CI.R',
 #                         vcol = NULL,
 #                         cvar=T,
 #                         metathr=0.1,
 #                         jobname="MetaVolcano",
 #                         outputfolder=".",
 #                         draw='HTML',
 #                         ncores=1)
#saveRDS(meta_degs_rem,'../MetaProjects/MetaDEG_v01/Output/mrem14_2.rds')
# #load previous run
meta_degs_rem=readRDS('../Output/mrem14_2.rds')
#meta_degs_rem@MetaVolcano #metavolcano plot

#result tables
mr=meta_degs_rem@metaresult
write.csv(mr, '../Output/MetaDEG_allgenes.csv')
#write.csv(mr,'../MetaProjects/MetaDEG_v01/Output/mr.csv')
#& the adjusted p value for meta results was calculated using the r p.adjust function with the parameter ('BH')
mrinput=meta_degs_rem@input
mr$adjP=p.adjust(mr$randomP,'BH')
#& meta deg genes were retained that were regulated in the same direction in at least 75% of the datasets (at least 11)
sc=.75*length(gse)
#use other volcano plotter
require(EnhancedVolcano)
mrsc=mr[abs(mr$signcon)>sc,]
EnhancedVolcano(mrsc,
                lab = mrsc$SYMBOL,
                x = 'randomSummary',
                y = 'randomP',
                xlim = c(-5, 8), title = 'Meta14',subtitle = '',transcriptLabSize = 2)
#diagnostics
par(mfrow=c(3,1))
hist(mr$signcon,breaks = 12)
hist(log10(mr$randomSummary),breaks = 12)
hist(mr$adjP)
combo=merge(mrinput,mr)
par(mfrow=c(1,1))

#8 plot forests ####
# & forest plots were created in using the metavolcanor function with  a adjusted p value of <0.05 and at least 11 with the same direction

#options**************************************************************************
 # number of sets to agree on sign
pthresh=.05#*** Numbers%@#@$
mrup=mr[mr$signcon >= sc& mr$adjP < pthresh & mr$randomSummary > 0,]
mrdown=mr[mr$signcon <=- sc & mr$adjP < pthresh & mr$randomSummary < 0,]
#*********************************************************************************

#plot.new();meta_degs_rem@MetaVolcano
pcriteria="P.Value"
foldchangecol='logFC'
genenamecol1="SYMBOL"
getwd()
##up
for (i in 1:2){
  require(MetaVolcanoR)
  pl=draw_forest(remres=meta_degs_rem,
                 gene=mrup$SYMBOL[i],
                 genecol=genenamecol1,
                 foldchangecol=foldchangecol,
                 llcol="CI.L",
                 rlcol="CI.R",
                 jobname="MetaVolcano",
                 outputfolder="../MetaProjects/MetaDEG_v01/Output/Forests",
                 draw="PDF")
  print(pl)
}
##down
for (i in 1:2){#*** Numbers%@#@$
  pl=draw_forest(remres=meta_degs_rem,
                 gene=mrdown$SYMBOL[i],
                 genecol=genenamecol1,
                 foldchangecol=foldchangecol,
                 llcol="CI.L",
                 rlcol="CI.R",
                 jobname="MetaVolcano",
                 outputfolder="../MetaProjects/MetaDEG_v01/Output/Forests",
                 draw="PDF")
  print(pl)
}



# reactome enrichment of different sets#####
require(vctrs)
require(clusterProfiler)
require(ReactomePA)
#&Gene set enrichment analysis was performed with the Reactome Pathways database (reactome.org) using the the R packages
#... ReactomePA and clusterProfiler using genes that have abs(log2(FC) > 1.5 and adjusted pvalue (benjamini) < .05

menrUp=enrichuniv(mrup$SYMBOL,univ = mr$SYMBOL)# mrup and mrdown defined above (line 247)
#& tables were created with the top 10 (lowest p.adjust) for enriched reactome pathways

# ggplot(menrUp %>% filter(p.adjust<.005), aes(y=-log10(p.adjust),x=Description))+ geom_bar(stat='identity')
# m2=melt(menrUp[menrUp$p.adjust<.005,c('Description','p.adjust')])
# ggplot(m2 , aes(x=Description,y=-log10(value), fill=(variable))) +geom_bar(position = 'dodge',stat='identity')+
#   theme(axis.text.x = element_text(
#     color='black',size=6, angle=45,vjust = 1, hjust=1))+ggtitle('Up')
menrDown=enrichuniv(mrdown$SYMBOL, univ=mr$SYMBOL)

#ne
# ## more forests#####
# for (i in d1[1:3]){
#   pl=draw_forest(remres=meta_degs_rem,
#                  gene=i,
#                  genecol=genenamecol1,
#                  foldchangecol=foldchangecol,
#                  llcol="CI.L",
#                  rlcol="CI.R",
#                  jobname="MetaVolcano",
#                  outputfolder="../MetaProjects/MetaDEG_v01/Output/Forests",
#                  draw="PDF")
#   print(pl)}
# 
# # 
# # ### zoom in on pathways ######
#  pcriteria="P.Value"
#  foldchangecol='logFC'
#  genenamecol1="SYMBOL"
#  getwd()
# #use meta analysis above as input
# #find pathway nubmer
# #pick row
#  grep('TCA',tolower(metenr$Description))
# ro=25
# pathname= metenr$Description[ro]
# print(pathname)
# genes1= metenr$geneID[ro]
# genes2=unlist(strsplit(genes1,'[/]'))#entrezid, only regulated genes
#  ent=bitr(genes2,toType='SYMBOL',fromType='ENTREZID', OrgDb = 'org.Hs.eg.db')[,2]#convert to symbol
#  for (i in 1:3){
#    pl=draw_forest(remres=meta_degs_rem,
#                   gene=ent[i],
#                   genecol=genenamecol1,
#                   foldchangecol=foldchangecol,
#                   llcol="CI.L",
#                   rlcol="CI.R",
#                   jobname="MetaVolcano",
#                   outputfolder="../MetaProjects/MetaDEG_v01/Output/Forests",
#                   draw="PDF")
#    print(pl)
#  }
#  #  pathway
#  
#  inp2=as.list(mr$randomSummary[mr$SYMBOL %in% ent])
#  names(inp2)=genes2
# inp2
# 
#  viewPathway(pathname, readable=TRUE,foldChange = inp2)
# 
# 
# 
#  
#  
#  rm(list=ls())
#  
#  
#  setwd('C://Users/bryan.linggi/Box Sync/PMED_Data(bryan.linggi@robartsinc.com 2)/')
#  # 
#  ## 
#  #mr=read.csv('MetaProjects/MetaDEG_v01/Output/mr.csv')# metadeg output from 14
#  ### test overlap of meta deg vs smillie deg in each subsets
#  #pseudocode
#  # 1. see if meta genes are present in smillie deg , 
#  # if yes, then what which cluster 
#  # if no, call meta unique
#  # 2. see if smillie deg in each cluster are in meta
#  # if yes, mark

 #get smillie dataset######
#& single cell RNA Seq data was analyzed from Smillie et al ref. The differential gene expression for scrNA was compared using all
#... data from Inflamed versus Healthy comparisons (Epithilial, Innate, and Adaptive cell subsets) from supplemental file s2, mmc4
#... overlap and differences between genes Differentially expressed by meta analysis (regulated in same direction 11 or more datasets, with a adjP <.05) and those in the scRNA were compred
 source("https://gist.github.com/schaunwheeler/5825002/raw/3526a15b032c06392740e20b6c9a179add2cee49/xlsxToR.r")

 sdeg=xlsxToR("../Code-Reference/Smillie2020/DEG_s2.0-S0092867419307329-mmc4.xlsx") #DEG ,  per scRNA cell cluster
 sdeg
 
 #fix layout of dataframes
 a1=sdeg[[1]]
 a=sdeg[c(7,8,9)]# get only sheets with inflamed vs healthy.  Epithelial, Innate, Adaptive
 colnames(a[[1]])
 a1=a[[1]]
 for (i in 1:3){
   colnames(a[[i]])=a[[i]][1,]
   a[[i]]=a[[i]][2:nrow(a[[i]]),]
   print(colnames(a[[i]]))
 }
 colnames(a[[3]])
 
 #put in single 
 sdeg2= Reduce(rbind,a[1:3])
 #split to up and down
 str(sdeg2)
 summary(as.numeric(sdeg2$padjD))
 sdegup=sdeg2[sdeg2$log2fc>0,]
 sdegdown=sdeg2[sdeg2$log2fc<0,]
 
 ###
 length(unique(sdeg2$gene))# =2676 only about 1/3 uinque, multiple cell types show deg 
 
 # look
 grep('IL1B',sdeg2$gene)
 #prep meta data 
 mr2=   mr[abs(mr$signcon)>sc & mr$adjP<.01,]      
 mr2up= mr2[mr2$randomSummary>0,]
 mr2down= mr2[mr2$randomSummary<0,]
 
 #visualize
 require(VennDiagram)
 venn.diagram(list(metaup=mr2up$SYMBOL, metadown=mr2down$SYMBOL, scrnaUp=sdegup$gene, scnrnaDown=sdegdown$gene),filename = '../MetaProjects/MetaDEG_v01/Output/scnrnacompdeg1.tiff')
 

 #check what cell cluster the uniquely scrna are from__down
 a=table(factor(sdegdown$ident[sdegdown$gene %in% setdiff(sdegdown$gene, mr2down$SYMBOL)]))
 d=sdegdown[sdegdown$gene %in% setdiff(sdegdown$gene, mr2down$SYMBOL),'ident']
 #and opposite
 b=table(factor(sdegdown$ident[sdegdown$gene %in% intersect(sdegdown$gene, mr2down$SYMBOL)]))
 e=sdegdown[sdegdown$gene %in% intersect(sdegdown$gene, mr2down$SYMBOL),'ident']
 #graph
 cmb=data.frame()
 cmb= data.frame(c(d,e))
 cmb$comparison=  c(rep('sc_only',length(d)),rep('intersect',length(e)))
 colnames(cmb)[1]='ident'
 require(ggplot2)
 ggplot(cmb, aes(y=ident,fill=comparison))+geom_bar()+ggtitle('down')
 
 #check what cell cluster the uniquely scrna are from_up
 a=table(factor(sdegup$ident[sdegup$gene %in% setdiff(sdegup$gene, mr2up$SYMBOL)]))
 d=sdegup[sdegup$gene %in% setdiff(sdegup$gene, mr2up$SYMBOL),'ident']
 #and opposite
 b=table(factor(sdegup$ident[sdegup$gene %in% intersect(sdegup$gene, mr2up$SYMBOL)]))
 e=sdegup[sdegup$gene %in% intersect(sdegup$gene, mr2up$SYMBOL),'ident']
 #graph
 cmb=data.frame()
 cmb= data.frame(c(d,e))
 cmb$comparison=  c(rep('sc_only',length(d)),rep('intersect',length(e)))
 colnames(cmb)[1]='ident'
 require(ggplot2)
 ggplot(cmb, aes(y=ident,fill=comparison))+geom_bar()+ggtitle('up')
 
 #get deg down , then up, for meta and then label with possible cell type using below
 # & for DEG meta genes not regulated in scRNA dataset, we used the data from mmc2 (padj <.05, >2fc)
 #...note that for meta, there is no fc cutoff (not sure how it translates to randomSummary) 
 #...see STAR methods 'Gene Specificity' from that paper), which contains genes and their assocaited cell
 #...type assignment, to assign the the presumptive cell type 
 #... derivation. 
 #get and prep ref 
 metadownonly=setdiff(mr2down$SYMBOL,sdegdown$gene)
 metauponly=setdiff(mr2up$SYMBOL,sdegup$gene)
 
 
 a=xlsxToR("../Code-Reference/Smillie2020/Copy of 1-s2.0-S0092867419307329-mmc2.xlsx") #gene signatures per scRNA cell cluster
 b=a
 #fix layout of dataframes
 a1=a[[1]]
 for (i in 1:3){
   colnames(a[[i]])=a[[i]][1,]
   a[[i]]=a[[i]][2:nrow(a[[i]]),]
   print(colnames(a[[i]]))
 }
 #put in single 
 dat= Reduce(rbind,a[1:3])
 
 dd=setdiff(mr2up$SYMBOL, sdegup$gene)
 sort(dd)
 d3=dat[dat$gene %in% setdiff( mr2up$SYMBOL, sdegup$gene),c('gene','ident')]
 mtchD=data.frame(metadownonly)
 mtchD$ident=NA
 for (i in 1:nrow(mtchD)){
   mtchD$ident[i]= dat$ident[match(mtchD$metadownonly[i], dat$gene)]
 }
 mtchD=mtchD[!is.na(mtchD$ident),]
 #add to cmb to replot , redo prep
 #check what cell cluster the uniquely scrna are from__down
 a=table(factor(sdegdown$ident[sdegdown$gene %in% setdiff(sdegdown$gene, mr2down$SYMBOL)]))
 d=sdegdown[sdegdown$gene %in% setdiff(sdegdown$gene, mr2down$SYMBOL),'ident']
 #and opposite
 b=table(factor(sdegdown$ident[sdegdown$gene %in% intersect(sdegdown$gene, mr2down$SYMBOL)]))
 e=sdegdown[sdegdown$gene %in% intersect(sdegdown$gene, mr2down$SYMBOL),'ident']
 #graph
 cmb=data.frame()
 cmb= data.frame(c(d,e,mtchD$ident))
 cmb$comparison=  c(rep('sc_only',length(d)),rep('intersect',length(e)),rep('metaonly',length(mtch$metadownonly)))
 colnames(cmb)[1]='ident'
 require(ggplot2)
 ggplot(cmb, aes(y=ident,fill=comparison))+geom_bar()+ggtitle('down')
 # 
 
 #for up----
 mtchU=data.frame(metauponly)
 mtchU$ident=NA
 for (i in 1:nrow(mtch)){
   mtchU$ident[i]= dat$ident[match(mtchU$metauponly[i], dat$gene)]
 }
 mtchU=mtchU[!is.na(mtchU$ident),]
 #add to cmb to replot , redo prep
 #check what cell cluster the uniquely scrna are from__up
 a=table(factor(sdegup$ident[sdegup$gene %in% setdiff(sdegup$gene, mr2up$SYMBOL)]))
 d=sdegup[sdegup$gene %in% setdiff(sdegup$gene, mr2up$SYMBOL),'ident']
 #and opposite
 b=table(factor(sdegup$ident[sdegup$gene %in% intersect(sdegup$gene, mr2up$SYMBOL)]))
 e=sdegup[sdegup$gene %in% intersect(sdegup$gene, mr2up$SYMBOL),'ident']
 #graph
 cmb=data.frame()
 cmb= data.frame(c(d,e,mtchU$ident))
 cmb$comparison=  c(rep('sc_only',length(d)),rep('intersect',length(e)),rep('metaonly',length(mtch$metauponly)))
 colnames(cmb)[1]='ident'
 require(ggplot2)
 ggplot(cmb, aes(y=ident,fill=comparison))+geom_bar()+ggtitle('up')
 
 ### compare selected genes to see difference###
 # top 5 from metaup
 t5= mrup %>% top_n(5,randomSummary)
 t5scU=sdeg2[sdeg2$gene %in% t5$SYMBOL,]
 #which not included
 setdiff(t5$SYMBOL,t5scU$gene)# reg1a
 # where is reg1a expressed
 dat$ident[grep('REG1A', dat$gene)]#TA 1, Cycling TA, this is not showing up for sc, is it consistent?
 
 # top 5 from metadown
 t5= mrdown %>% top_n(-5,randomSummary)
 t5scD=sdeg2[sdeg2$gene %in% t5$SYMBOL,]
 #which not included
 a=print(setdiff(t5$SYMBOL,t5scD$gene))# 
 # where is reg1a expressed
 dat$ident[grep(a[1],dat$gene)]#CLDN8,  "E.Secretory" "Goblet" 
 dat$ident[grep(a[2],dat$gene)]#SLC26A2  "TA 2"  "Immature Enterocytes 1" "E.Epithelial" "Enterocyte Progenitors"
 #..."Immature Enterocytes 2" "E.Absorptive_All"       "E.Immature_Enterocytes"
 # general pattern, absorbtive, enterocytes, secretory, those deg not showing for E. sec and Gobl.
 
 ###top 5 up from sc (doesnt have down, list, is pairwise comparison both directions)
 t5sU= sdeg2 %>% top_n(5, log2fc) #"NCS1"     "RBPMS"    "SERPINB5" "TNIP3"    "SAA1"   
print( t5sU[,c('gene','ident')])
m5U=mrup[mrup$SYMBOL %in% t5sU$gene,]# serpinB5 (E. Imm Enter) and TNIP3 (Enterocyets) there, not 
#... the others ,which are cycling TA , E Epithelial, and TA 2
m5D=mrdown[mrdown$SYMBOL %in% t5sU$gene,]# none down 

# ! check to make sure genes are in total list, not just annotaion, microarray/seq differences


 
# others of interest
 # osm, osmr, trem1, il13ra2
 p1=mr[grep('OSM',mr$SYMBOL)[1],]
 p11=sdeg2[grep('OSM',sdeg2$gene)[1],]#not there
 p111=dat$ident[grep('OSM',dat$gene)]#cell type?  "Inflammatory Monocytes"
 
 p2=mr[grep('OSMR',mr$SYMBOL)[1],]
 p22=sdeg2[grep('OSMR',sdeg2$gene),]#not there
 p222=dat$ident[grep('OSMR',dat$gene)]#cell type, none
 
 
 p3=mr[grep('TREM1',mr$SYMBOL),]
 p33=sdeg2[grep('TREM1',sdeg2$gene),]#not there
 p333=dat$ident[grep('TREM1',dat$gene)]#cell type? "DC2"    "Inflammatory Monocytes"  "M.DCs"  
 
 
 