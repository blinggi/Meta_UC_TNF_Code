# Meta Deg  ####
# from UC disease sig analysis
rm(list=ls())
require(limma);require(Biobase);require(BiocGenerics)
wd= setwd('C:\\Users/bryan.linggi/Box Sync/PMED_Data(bryan.linggi@robartsinc.com 2)/Projects/')
source('C:/Users/bryan.linggi/Box Sync/PMED_Data(bryan.linggi@robartsinc.com 2)/Code_General/DEGfunctions_v01.R')

####get DEG datasets, output from DEGout codes#####
gse16=readRDS('GSE16879/Output/GSE16879_coloncontrol_DEG.rds')#redid 23mar2020 to include all probes, (before had removed na symbol, but had not done same in other datasets)

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

