######################## alternative splicing analysis using BW vs DLN as an example ##################################

library("reshape")

#A3SS
A3SS=read.table("A3SS.MATS.JC.txt",header=TRUE)
A3SS=transform(A3SS,IJC_SAMPLE_1=colsplit(IJC_SAMPLE_1,split=",",names=c('BW_NA_10','BW_NA_6','BW_NA_7')))
A3SS=transform(A3SS,SJC_SAMPLE_1=colsplit(SJC_SAMPLE_1,split=",",names=c('BW_NA_10','BW_NA_6','BW_NA_7')))
A3SS=transform(A3SS,IJC_SAMPLE_2=colsplit(IJC_SAMPLE_2,split=",",names=c('DLN_NA_1','DLN_NA_2','DLN_NA_4')))
A3SS=transform(A3SS,SJC_SAMPLE_2=colsplit(SJC_SAMPLE_2,split=",",names=c('DLN_NA_1','DLN_NA_2','DLN_NA_4')))
A3SS=transform(A3SS,IncLevel1=colsplit(IncLevel1,split=",",names=c('BW_NA_10','BW_NA_6','BW_NA_7')))
A3SS=transform(A3SS,IncLevel2=colsplit(IncLevel2,split=",",names=c('DLN_NA_1','DLN_NA_2','DLN_NA_4')))

A3SS$IJC_sum=rowSums(A3SS$IJC_SAMPLE_1)+rowSums(A3SS$IJC_SAMPLE_2)
A3SS$SJC_sum=rowSums(A3SS$SJC_SAMPLE_1)+rowSums(A3SS$SJC_SAMPLE_2)
A3SS_count_over_20=A3SS[A3SS$IJC_sum>20&A3SS$SJC_sum>20,] #confident AS events

A3SS_sig=A3SS_count_over_20[A3SS_count_over_20$FDR<0.05,] #DS events
nrow(A3SS_sig)
write.csv(A3SS_sig,file="A3SS_DS.csv")

#A5SS
A5SS=read.table("A5SS.MATS.JC.txt",header=TRUE)
A5SS=transform(A5SS,IJC_SAMPLE_1=colsplit(IJC_SAMPLE_1,split=",",names=c('BW_NA_10','BW_NA_6','BW_NA_7')))
A5SS=transform(A5SS,SJC_SAMPLE_1=colsplit(SJC_SAMPLE_1,split=",",names=c('BW_NA_10','BW_NA_6','BW_NA_7')))
A5SS=transform(A5SS,IJC_SAMPLE_2=colsplit(IJC_SAMPLE_2,split=",",names=c('DLN_NA_1','DLN_NA_2','DLN_NA_4')))
A5SS=transform(A5SS,SJC_SAMPLE_2=colsplit(SJC_SAMPLE_2,split=",",names=c('DLN_NA_1','DLN_NA_2','DLN_NA_4')))
A5SS=transform(A5SS,IncLevel1=colsplit(IncLevel1,split=",",names=c('BW_NA_10','BW_NA_6','BW_NA_7')))
A5SS=transform(A5SS,IncLevel2=colsplit(IncLevel2,split=",",names=c('DLN_NA_1','DLN_NA_2','DLN_NA_4')))

A5SS$IJC_sum=rowSums(A5SS$IJC_SAMPLE_1)+rowSums(A5SS$IJC_SAMPLE_2)
A5SS$SJC_sum=rowSums(A5SS$SJC_SAMPLE_1)+rowSums(A5SS$SJC_SAMPLE_2)
A5SS_count_over_20=A5SS[A5SS$IJC_sum>20&A5SS$SJC_sum>20,] #confident AS events

A5SS_sig=A5SS_count_over_20[A5SS_count_over_20$FDR<0.05,] #DS events
nrow(A5SS_sig)
write.csv(A5SS_sig,file="A5SS_DS.csv")

#SE
SE=read.table("SE.MATS.JC.txt",header=TRUE)
SE=transform(SE,IJC_SAMPLE_1=colsplit(IJC_SAMPLE_1,split=",",names=c('BW_NA_10','BW_NA_6','BW_NA_7')))
SE=transform(SE,SJC_SAMPLE_1=colsplit(SJC_SAMPLE_1,split=",",names=c('BW_NA_10','BW_NA_6','BW_NA_7')))
SE=transform(SE,IJC_SAMPLE_2=colsplit(IJC_SAMPLE_2,split=",",names=c('DLN_NA_1','DLN_NA_2','DLN_NA_4')))
SE=transform(SE,SJC_SAMPLE_2=colsplit(SJC_SAMPLE_2,split=",",names=c('DLN_NA_1','DLN_NA_2','DLN_NA_4')))
SE=transform(SE,IncLevel1=colsplit(IncLevel1,split=",",names=c('BW_NA_10','BW_NA_6','BW_NA_7')))
SE=transform(SE,IncLevel2=colsplit(IncLevel2,split=",",names=c('DLN_NA_1','DLN_NA_2','DLN_NA_4')))

SE$IJC_sum=rowSums(SE$IJC_SAMPLE_1)+rowSums(SE$IJC_SAMPLE_2)
SE$SJC_sum=rowSums(SE$SJC_SAMPLE_1)+rowSums(SE$SJC_SAMPLE_2)
SE_count_over_20=SE[SE$IJC_sum>20&SE$SJC_sum>20,] #confident AS events

SE_sig=SE_count_over_20[SE_count_over_20$FDR<0.05,] #DS events
nrow(SE_sig)
write.csv(SE_sig,file="SE_DS.csv")

#MXE
MXE=read.table("MXE.MATS.JC.txt",header=TRUE)
MXE=transform(MXE,IJC_SAMPLE_1=colsplit(IJC_SAMPLE_1,split=",",names=c('BW_NA_10','BW_NA_6','BW_NA_7')))
MXE=transform(MXE,SJC_SAMPLE_1=colsplit(SJC_SAMPLE_1,split=",",names=c('BW_NA_10','BW_NA_6','BW_NA_7')))
MXE=transform(MXE,IJC_SAMPLE_2=colsplit(IJC_SAMPLE_2,split=",",names=c('DLN_NA_1','DLN_NA_2','DLN_NA_4')))
MXE=transform(MXE,SJC_SAMPLE_2=colsplit(SJC_SAMPLE_2,split=",",names=c('DLN_NA_1','DLN_NA_2','DLN_NA_4')))
MXE=transform(MXE,IncLevel1=colsplit(IncLevel1,split=",",names=c('BW_NA_10','BW_NA_6','BW_NA_7')))
MXE=transform(MXE,IncLevel2=colsplit(IncLevel2,split=",",names=c('DLN_NA_1','DLN_NA_2','DLN_NA_4')))

MXE$IJC_sum=rowSums(MXE$IJC_SAMPLE_1)+rowSums(MXE$IJC_SAMPLE_2)
MXE$SJC_sum=rowSums(MXE$SJC_SAMPLE_1)+rowSums(MXE$SJC_SAMPLE_2)
MXE_count_over_20=MXE[MXE$IJC_sum>20&MXE$SJC_sum>20,] #confident AS events

MXE_sig=MXE_count_over_20[MXE_count_over_20$FDR<0.05,] #DS events
nrow(MXE_sig)
write.csv(MXE_sig,file="MXE_DS.csv")

#RI
RI=read.table("RI.MATS.JC.txt",header=TRUE)
RI=transform(RI,IJC_SAMPLE_1=colsplit(IJC_SAMPLE_1,split=",",names=c('BW_NA_10','BW_NA_6','BW_NA_7')))
RI=transform(RI,SJC_SAMPLE_1=colsplit(SJC_SAMPLE_1,split=",",names=c('BW_NA_10','BW_NA_6','BW_NA_7')))
RI=transform(RI,IJC_SAMPLE_2=colsplit(IJC_SAMPLE_2,split=",",names=c('DLN_NA_1','DLN_NA_2','DLN_NA_4')))
RI=transform(RI,SJC_SAMPLE_2=colsplit(SJC_SAMPLE_2,split=",",names=c('DLN_NA_1','DLN_NA_2','DLN_NA_4')))
RI=transform(RI,IncLevel1=colsplit(IncLevel1,split=",",names=c('BW_NA_10','BW_NA_6','BW_NA_7')))
RI=transform(RI,IncLevel2=colsplit(IncLevel2,split=",",names=c('DLN_NA_1','DLN_NA_2','DLN_NA_4')))

RI$IJC_sum=rowSums(RI$IJC_SAMPLE_1)+rowSums(RI$IJC_SAMPLE_2)
RI$SJC_sum=rowSums(RI$SJC_SAMPLE_1)+rowSums(RI$SJC_SAMPLE_2)
RI_count_over_20=RI[RI$IJC_sum>20&RI$SJC_sum>20,] #confident AS events

RI_sig=RI_count_over_20[RI_count_over_20$FDR<0.05,] #DS events
nrow(RI_sig)
write.csv(RI_sig,file="RI_DS.csv")

## AS events
A3SS_count_over_20$Type=rep("A3SS",nrow(A3SS_count_over_20))
A5SS_count_over_20$Type=rep("A5SS",nrow(A5SS_count_over_20))
MXE_count_over_20$Type=rep("MXE",nrow(MXE_count_over_20))
SE_count_over_20$Type=rep("SE",nrow(SE_count_over_20))
RI_count_over_20$Type=rep("RI",nrow(RI_count_over_20))

AS=rbind(rbind(rbind(rbind(A3SS_count_over_20[,c(2,20,23,26)]
                           ,A5SS_count_over_20[,c(2,20,23,26)])
                     ,MXE_count_over_20[,c(2,22,25,28)])
               ,SE_count_over_20[,c(2,20,23,26)])
         ,RI_count_over_20[,c(2,20,23,26)])
AS$absIncleveldifference=abs(AS$IncLevelDifference)
write.csv(AS,file="AS_events.csv") # Total confident AS events between two groups 16453
AS_genes=AS[!duplicated(AS$GeneID),]

####DSG
A3SS_sig$Type=rep("A3SS",nrow(A3SS_sig))
A5SS_sig$Type=rep("A5SS",nrow(A5SS_sig))
MXE_sig$Type=rep("MXE",nrow(MXE_sig))
SE_sig$Type=rep("SE",nrow(SE_sig))
RI_sig$Type=rep("RI",nrow(RI_sig))

DSG=rbind(rbind(rbind(rbind(A3SS_sig[,c(2,20,23,26)]
                            ,A5SS_sig[,c(2,20,23,26)])
                      ,MXE_sig[,c(2,22,25,28)])
                ,SE_sig[,c(2,20,23,26)])
          ,RI_sig[,c(2,20,23,26)])
DSG$absIncleveldifference=abs(DSG$IncLevelDifference)

# The DS event with the largest absIncleveldifference of a DSG represent the differential splicing level of the gene.
DSG_unique_gene=DSG[order(DSG[,'GeneID'],-DSG[,'absIncleveldifference']),]
DSG_unique_gene=DSG_unique_gene[!duplicated(DSG_unique_gene$GeneID),]
write.csv(DSG_unique_gene,file="DSG_unique_BW_DLN.csv") # Total unique DS genes between two groups
