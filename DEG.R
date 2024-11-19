library(DESeq2)

#Read in sample table and Trim counts table
counts=read.table("woodfrog_rna.txt", header = T, stringsAsFactors = F)
rownames(counts)=counts$ENSEMBL_GeneID
counts=counts[,-1]
noint = rownames(counts) %in% c("__no_feature","__ambiguous","__too_low_aQual", "__not_aligned","__alignment_not_unique")
samples=read.csv("sample.csv", stringsAsFactors = F)
samples=samples[samples$name %in% colnames(counts),]
samples=samples[order(samples$name),]
samples=samples[-1,] #remove outlier sample BW_F_113_E 
counts=counts[,-1] #remove outlier sample BW_F_113_E

# remove non-informative rows and weakly expressed genes in all samples
counts=counts[!noint,]
colData=data.frame(samples[,c("stage","type")])
rownames(colData)=samples$name
colnames(colData)=c("stage","type")
colnames(counts)==rownames(colData) #check the order of sample table and count table

dds=DESeqDataSetFromMatrix(countData = counts, 
                           colData = colData,
                           design = ~ stage + type + stage * type) #interaction between life stage and type

keep_gene=rowSums(counts(dds)>=1)>=nrow(colData) # remove genes with less than 1 read in all samples
dds=dds[keep_gene,] 

####DE analysis of adult vs. hatching####
dds$stage = factor(dds$stage, levels = c("adult","hatchling"))
dds$type = factor(dds$type, levels=c("roadside","woodland"))

colnames(dds)==rownames(colData) # check if the count table has the same order as sample table

# Test for DE genes
dds_stage=DESeq(dds)
res_stage=results(dds_stage, contrast = c("stage", "adult", "hatchling"))   #test DE for adult versus hatchling while controlling for habitat type differences
head(res_stage)

# Inspect the result tables of significantly upregulated and downregulated genes, at a 5% false discovery rate (FDR) as follows:
resSig_stage=res_stage[which(res_stage$padj<0.05),]
head(resSig_stage[order(resSig_stage$log2FoldChange, decreasing = TRUE), ])
head(resSig_stage[order(resSig_stage$log2FoldChange, decreasing = FALSE), ])

write.csv(resSig_stage[order(resSig_stage$log2FoldChange, decreasing = TRUE), ],file = "DE_stage_23.csv") #control for difference in habitat type
deg_stage=rownames(resSig_stage) 
# Count the number of genes with significant differential expression at a FDR of 5%:
table(res_stage$padj < 0.05) 


#### DE analysis of roadside vs woodland ####
#trim sample table
samples_keep=grep("hatchling",samples$stage)
samples_hatchling=samples[samples_keep,] #remain hatchling only
#trim count table
counts_keep=grep("NA",colnames(counts))
counts_hatchling=counts[,counts_keep] #remain hatchling only

dds_hatchling=dds[,dds$stage=="hatchling"]

# DE analysis between 3 x 3 roadside vs woodland ----------
##### DE analysis of BW vs RO #####
#Trim sample table for BW vs RO
keep_samples_BW=grep("BW", samples_hatchling$pond)
samples_bw=samples_hatchling[keep_samples_BW,]
keep_samples_RO=grep("RO", samples_hatchling$pond)
samples_ro=samples_hatchling[keep_samples_RO,]

samples_bw_ro=rbind(samples_bw,samples_ro)

#Trim counts table for DE analysis between BW vs RO
bw_keep=grep("BW", colnames(counts_hatchling)) 
bw=dds_hatchling[,bw_keep]

ro_keep=grep("RO",colnames(counts_hatchling))
ro=dds_hatchling[,ro_keep]

bw_ro=cbind(bw@assays@data@listData$counts,ro@assays@data@listData$counts) #combine two count table and then make DEseqdataset from matrix
colData_bw_ro=data.frame(samples_bw_ro[,c("pond")])
rownames(colData_bw_ro)=samples_bw_ro$name
colnames(colData_bw_ro)="pond"

colnames(bw_ro)==rownames(colData_bw_ro)  # check if the count table has the same order as sample table

dds_bw_ro=DESeqDataSetFromMatrix(countData = bw_ro, 
                                 colData = colData_bw_ro,
                                 design = ~ pond) #making DESeqDataSet for BW and RO only


dds_bw_ro$pond = factor(dds_bw_ro$pond, levels = c("BW","RO"))

# Test for DE genes
dds_bw_ro_1=DESeq(dds_bw_ro)
res_bw_ro=results(dds_bw_ro_1, contrast = c("pond", "BW", "RO"))   #test DE for BW versus RO
head(res_bw_ro)


# Inspect the result tables of significantly upregulated and downregulated genes, at a 5% false discovery rate (FDR) as follows:
resSig_bw_ro=res_bw_ro[which(res_bw_ro$padj<0.05),]
head(resSig_bw_ro[order(resSig_bw_ro$log2FoldChange, decreasing = TRUE), ])
head(resSig_bw_ro[order(resSig_bw_ro$log2FoldChange, decreasing = FALSE), ])

write.csv(resSig_bw_ro[order(resSig_bw_ro$log2FoldChange, decreasing = TRUE), ],file = "DE_BW_vs_RO_all.csv")
deg_bw_ro=rownames(resSig_bw_ro)

# Count the number of genes with significant differential expression at a FDR of 5%:
table(resSig_bw_ro$padj < 0.05)

##### DE analysis of BW vs DLN####
#Trim sample table for BW vs DLN
keep_samples_BW=grep("BW", samples_hatchling$pond)
samples_BW=samples_hatchling[keep_samples_BW,]
keep_samples_DLN=grep("DLN", samples_hatchling$pond)
samples_DLN=samples_hatchling[keep_samples_DLN,]

samples_BW_DLN=rbind(samples_BW,samples_DLN)

#Trim counts table for DE analysis between BW vs DLN
BW_keep=grep("BW", colnames(counts_hatchling)) 
BW=dds_hatchling[,BW_keep]

DLN_keep=grep("DLN",colnames(counts_hatchling))
DLN=dds_hatchling[,DLN_keep]

BW_DLN=cbind(BW@assays@data@listData$counts,DLN@assays@data@listData$counts) #combine two count table and then make DEseqdataset from matrix
colData_BW_DLN=data.frame(samples_BW_DLN[,c("pond")])
rownames(colData_BW_DLN)=samples_BW_DLN$name
colnames(colData_BW_DLN)="pond"

colnames(BW_DLN)==rownames(colData_BW_DLN)  # check if the count table has the same order as sample table

dds_BW_DLN=DESeqDataSetFromMatrix(countData = BW_DLN, 
                                  colData = colData_BW_DLN,
                                  design = ~ pond) #making DESeqDataSet for BW and DLN only


dds_BW_DLN$pond = factor(dds_BW_DLN$pond, levels = c("BW","DLN"))

# Test for DE genes
dds_BW_DLN_1=DESeq(dds_BW_DLN)
res_BW_DLN=results(dds_BW_DLN_1, contrast = c("pond", "BW", "DLN"))   #test DE for BW versus DLN
head(res_BW_DLN)


# Inspect the result tables of significantly upregulated and downregulated genes, at a 5% false discovery rate (FDR) as follows:
resSig_BW_DLN=res_BW_DLN[which(res_BW_DLN$padj<0.05),]
head(resSig_BW_DLN[order(resSig_BW_DLN$log2FoldChange, decreasing = TRUE), ])
head(resSig_BW_DLN[order(resSig_BW_DLN$log2FoldChange, decreasing = FALSE), ])

write.csv(resSig_BW_DLN[order(resSig_BW_DLN$log2FoldChange, decreasing = TRUE), ],file = "DE_BW_vs_DLN_all.csv")
deg_BW_DLN=rownames(resSig_BW_DLN)

# Count the number of genes with significant differential expression at a FDR of 5%:
table(resSig_BW_DLN$padj < 0.05)


##### DE analysis of BW vs ST####
#Trim sample table for BW vs ST
keep_samples_BW=grep("BW", samples_hatchling$pond)
samples_BW=samples_hatchling[keep_samples_BW,]
keep_samples_ST=grep("ST", samples_hatchling$pond)
samples_ST=samples_hatchling[keep_samples_ST,]

samples_BW_ST=rbind(samples_BW,samples_ST)

#Trim counts table for DE analysis between BW vs ST
BW_keep=grep("BW", colnames(counts_hatchling)) 
BW=dds_hatchling[,BW_keep]

ST_keep=grep("ST",colnames(counts_hatchling))
ST=dds_hatchling[,ST_keep]

BW_ST=cbind(BW@assays@data@listData$counts,ST@assays@data@listData$counts) #combine two count table and then make DEseqdataset from matrix
colData_BW_ST=data.frame(samples_BW_ST[,c("pond")])
rownames(colData_BW_ST)=samples_BW_ST$name
colnames(colData_BW_ST)="pond"

colnames(BW_ST)==rownames(colData_BW_ST)  # check if the count table has the same order as sample table

dds_BW_ST=DESeqDataSetFromMatrix(countData = BW_ST, 
                                 colData = colData_BW_ST,
                                 design = ~ pond) #making DESeqDataSet for BW and ST only


dds_BW_ST$pond = factor(dds_BW_ST$pond, levels = c("BW","ST"))

# Test for DE genes
dds_BW_ST_1=DESeq(dds_BW_ST)
res_BW_ST=results(dds_BW_ST_1, contrast = c("pond", "BW", "ST"))   #test DE for BW versus ST
head(res_BW_ST)


# Inspect the result tables of significantly upregulated and downregulated genes, at a 5% false discovery rate (FDR) as follows:
resSig_BW_ST=res_BW_ST[which(res_BW_ST$padj<0.05),]
head(resSig_BW_ST[order(resSig_BW_ST$log2FoldChange, decreasing = TRUE), ])
head(resSig_BW_ST[order(resSig_BW_ST$log2FoldChange, decreasing = FALSE), ])

write.csv(resSig_BW_ST[order(resSig_BW_ST$log2FoldChange, decreasing = TRUE), ],file = "DE_BW_vs_ST_all.csv")
deg_BW_ST=rownames(resSig_BW_ST)

# Count the number of genes with significant differential expression at a FDR of 5%:
table(resSig_BW_ST$padj < 0.05)


##### DE analysis of RSB vs RO ####
#Trim sample table for RSB vs RO
keep_samples_RSB=grep("RSB", samples_hatchling$pond)
samples_RSB=samples_hatchling[keep_samples_RSB,]
keep_samples_RO=grep("RO", samples_hatchling$pond)
samples_RO=samples_hatchling[keep_samples_RO,]

samples_RSB_RO=rbind(samples_RSB,samples_RO)

#Trim counts table for DE analysis between RSB vs RO
RSB_keep=grep("RSB", colnames(counts_hatchling)) 
RSB=dds_hatchling[,RSB_keep]

RO_keep=grep("RO",colnames(counts_hatchling))
RO=dds_hatchling[,RO_keep]

RSB_RO=cbind(RSB@assays@data@listData$counts,RO@assays@data@listData$counts) #combine two count table and then make DEseqdataset from matrix
colData_RSB_RO=data.frame(samples_RSB_RO[,c("pond")])
rownames(colData_RSB_RO)=samples_RSB_RO$name
colnames(colData_RSB_RO)="pond"

colnames(RSB_RO)==rownames(colData_RSB_RO)  # check if the count table has the same order as sample table

dds_RSB_RO=DESeqDataSetFromMatrix(countData = RSB_RO, 
                                  colData = colData_RSB_RO,
                                  design = ~ pond) #making DESeqDataSet for RSB and RO only


dds_RSB_RO$pond = factor(dds_RSB_RO$pond, levels = c("RO","RSB"))

# Test for DE genes
dds_RSB_RO_1=DESeq(dds_RSB_RO)
res_RSB_RO=results(dds_RSB_RO_1, contrast = c("pond", "RSB", "RO"))   #test DE for RSB versus RO
head(res_RSB_RO)


# Inspect the result tables of significantly upregulated and downregulated genes, at a 5% false discovery rate (FDR) as follows:
resSig_RSB_RO=res_RSB_RO[which(res_RSB_RO$padj<0.05),]
head(resSig_RSB_RO[order(resSig_RSB_RO$log2FoldChange, decreasing = TRUE), ])
head(resSig_RSB_RO[order(resSig_RSB_RO$log2FoldChange, decreasing = FALSE), ])

write.csv(resSig_RSB_RO[order(resSig_RSB_RO$log2FoldChange, decreasing = TRUE), ],file = "DE_RSB_vs_RO_all.csv")
deg_RSB_RO=rownames(resSig_RSB_RO)

# Count the number of genes with significant differential expression at a FDR of 5%:
table(resSig_RSB_RO$padj < 0.05)


##### DE analysis of RSB vs DLN #####
#Trim sample table for RSB vs DLN 
keep_samples_RSB=grep("RSB", samples_hatchling$pond)
samples_rsb=samples_hatchling[keep_samples_RSB,]
keep_samples_DLN=grep("DLN", samples_hatchling$pond)
samples_dln=samples_hatchling[keep_samples_DLN,]

samples_rsb_dln=rbind(samples_rsb,samples_dln)

#Trim counts table for DE analysis between RSB vs DLN 

rsb_keep=grep("RSB", colnames(counts_hatchling)) 
rsb=dds_hatchling[,rsb_keep]

dln_keep=grep("DLN",colnames(counts_hatchling))
dln=dds_hatchling[,dln_keep]

rsb_dln=cbind(rsb@assays@data@listData$counts,dln@assays@data@listData$counts) #combine two count table and then make DEseqdataset from matrix
colData_rsb_dln=data.frame(samples_rsb_dln[,c("pond")])
rownames(colData_rsb_dln)=samples_rsb_dln$name
colnames(colData_rsb_dln)="pond"

colnames(rsb_dln)==rownames(colData_rsb_dln)  # check if the count table has the same order as sample table

dds_rsb_dln=DESeqDataSetFromMatrix(countData = rsb_dln, 
                                   colData = colData_rsb_dln,
                                   design = ~ pond) #making DESeqDataSet for RSB and DLN only


dds_rsb_dln$pond = factor(dds_rsb_dln$pond, levels = c("DLN","RSB"))

# Test for DE genes
dds_rsb_dln_1=DESeq(dds_rsb_dln)
res_rsb_dln=results(dds_rsb_dln_1, contrast = c("pond", "RSB", "DLN"))   #test DE for RSB versus DLN
head(res_rsb_dln)


# Inspect the result tables of significantly upregulated and downregulated genes, at a 5% false discovery rate (FDR) as follows:
resSig_rsb_dln=res_rsb_dln[which(res_rsb_dln$padj<0.05),]
head(resSig_rsb_dln[order(resSig_rsb_dln$log2FoldChange, decreasing = TRUE), ])
head(resSig_rsb_dln[order(resSig_rsb_dln$log2FoldChange, decreasing = FALSE), ])

write.csv(resSig_rsb_dln[order(resSig_rsb_dln$log2FoldChange, decreasing = TRUE), ],file = "DE_RSB_vs_DLN_all.csv")
deg_rsb_dln=rownames(resSig_rsb_dln)

# Count the number of genes with significant differential expression at a FDR of 5%:
table(resSig_rsb_dln$padj < 0.05)


##### DE analysis of RSB vs ST ####
#Trim sample table for RSB vs ST
keep_samples_RSB=grep("RSB", samples_hatchling$pond)
samples_RSB=samples_hatchling[keep_samples_RSB,]
keep_samples_ST=grep("ST", samples_hatchling$pond)
samples_ST=samples_hatchling[keep_samples_ST,]

samples_RSB_ST=rbind(samples_RSB,samples_ST)

#Trim counts table for DE analysis between RSB vs ST
RSB_keep=grep("RSB", colnames(counts_hatchling)) 
RSB=dds_hatchling[,RSB_keep]

ST_keep=grep("ST",colnames(counts_hatchling))
ST=dds_hatchling[,ST_keep]

RSB_ST=cbind(RSB@assays@data@listData$counts,ST@assays@data@listData$counts) #combine two count table and then make DEseqdataset from matrix
colData_RSB_ST=data.frame(samples_RSB_ST[,c("pond")])
rownames(colData_RSB_ST)=samples_RSB_ST$name
colnames(colData_RSB_ST)="pond"

colnames(RSB_ST)==rownames(colData_RSB_ST)  # check if the count table has the same order as sample table

dds_RSB_ST=DESeqDataSetFromMatrix(countData = RSB_ST, 
                                  colData = colData_RSB_ST,
                                  design = ~ pond) #making DESeqDataSet for RSB and ST only


dds_RSB_ST$pond = factor(dds_RSB_ST$pond, levels = c("RSB","ST"))

# Test for DE genes
dds_RSB_ST_1=DESeq(dds_RSB_ST)
res_RSB_ST=results(dds_RSB_ST_1, contrast = c("pond", "RSB", "ST"))   #test DE for RSB versus ST
head(res_RSB_ST)


# Inspect the result tables of significantly upregulated and downregulated genes, at a 5% false discovery rate (FDR) as follows:
resSig_RSB_ST=res_RSB_ST[which(res_RSB_ST$padj<0.05),]
head(resSig_RSB_ST[order(resSig_RSB_ST$log2FoldChange, decreasing = TRUE), ])
head(resSig_RSB_ST[order(resSig_RSB_ST$log2FoldChange, decreasing = FALSE), ])

write.csv(resSig_RSB_ST[order(resSig_RSB_ST$log2FoldChange, decreasing = TRUE), ],file = "DE_RSB_vs_ST_all.csv")
deg_RSB_ST=rownames(resSig_RSB_ST)

# Count the number of genes with significant differential expression at a FDR of 5%:
table(resSig_RSB_ST$padj < 0.05)



##### DE analysis of EC vs RO ####
#Trim sample table for EC vs RO
keep_samples_EC=grep("EC", samples_hatchling$pond)
samples_EC=samples_hatchling[keep_samples_EC,]
keep_samples_RO=grep("RO", samples_hatchling$pond)
samples_RO=samples_hatchling[keep_samples_RO,]

samples_EC_RO=rbind(samples_EC,samples_RO)

#Trim counts table for DE analysis between EC vs RO
EC_keep=grep("EC", colnames(counts_hatchling)) 
EC=dds_hatchling[,EC_keep]

RO_keep=grep("RO",colnames(counts_hatchling))
RO=dds_hatchling[,RO_keep]

EC_RO=cbind(EC@assays@data@listData$counts,RO@assays@data@listData$counts) #combine two count table and then make DEseqdataset from matrix
colData_EC_RO=data.frame(samples_EC_RO[,c("pond")])
rownames(colData_EC_RO)=samples_EC_RO$name
colnames(colData_EC_RO)="pond"

colnames(EC_RO)==rownames(colData_EC_RO)  # check if the count table has the same order as sample table

dds_EC_RO=DESeqDataSetFromMatrix(countData = EC_RO, 
                                 colData = colData_EC_RO,
                                 design = ~ pond) #making DESeqDataSet for EC and RO only


dds_EC_RO$pond = factor(dds_EC_RO$pond, levels = c("EC","RO"))

# Test for DE genes
dds_EC_RO_1=DESeq(dds_EC_RO)
res_EC_RO=results(dds_EC_RO_1, contrast = c("pond", "EC", "RO"))   #test DE for EC versus RO
head(res_EC_RO)


# Inspect the result tables of significantly upregulated and downregulated genes, at a 5% false discovery rate (FDR) as follows:
resSig_EC_RO=res_EC_RO[which(res_EC_RO$padj<0.05),]
head(resSig_EC_RO[order(resSig_EC_RO$log2FoldChange, decreasing = TRUE), ])
head(resSig_EC_RO[order(resSig_EC_RO$log2FoldChange, decreasing = FALSE), ])

write.csv(resSig_EC_RO[order(resSig_EC_RO$log2FoldChange, decreasing = TRUE), ],file = "DE_EC_vs_RO_all.csv")
deg_EC_RO=rownames(resSig_EC_RO)

# Count the number of genes with significant differential expression at a FDR of 5%:
table(resSig_EC_RO$padj < 0.05)


##### DE analysis of EC vs DLN #####
#Trim sample table for EC vs DLN
keep_samples_EC=grep("EC", samples_hatchling$pond)
samples_EC=samples_hatchling[keep_samples_EC,]
keep_samples_DLN=grep("DLN", samples_hatchling$pond)
samples_DLN=samples_hatchling[keep_samples_DLN,]

samples_EC_DLN=rbind(samples_EC,samples_DLN)

#Trim counts table for DE analysis between EC vs DLN
EC_keep=grep("EC", colnames(counts_hatchling)) 
EC=dds_hatchling[,EC_keep]

DLN_keep=grep("DLN",colnames(counts_hatchling))
DLN=dds_hatchling[,DLN_keep]

EC_DLN=cbind(EC@assays@data@listData$counts,DLN@assays@data@listData$counts) #combine two count table and then make DEseqdataset from matrix
colData_EC_DLN=data.frame(samples_EC_DLN[,c("pond")])
rownames(colData_EC_DLN)=samples_EC_DLN$name
colnames(colData_EC_DLN)="pond"

colnames(EC_DLN)==rownames(colData_EC_DLN)  # check if the count table has the same order as sample table

dds_EC_DLN=DESeqDataSetFromMatrix(countData = EC_DLN, 
                                  colData = colData_EC_DLN,
                                  design = ~ pond) #making DESeqDataSet for EC and DLN only


dds_EC_DLN$pond = factor(dds_EC_DLN$pond, levels = c("DLN","EC"))

# Test for DE genes
dds_EC_DLN_1=DESeq(dds_EC_DLN)
res_EC_DLN=results(dds_EC_DLN_1, contrast = c("pond", "EC", "DLN"))   #test DE for EC versus DLN
head(res_EC_DLN)


# Inspect the result tables of significantly upregulated and downregulated genes, at a 5% false discovery rate (FDR) as follows:
resSig_EC_DLN=res_EC_DLN[which(res_EC_DLN$padj<0.05),]
head(resSig_EC_DLN[order(resSig_EC_DLN$log2FoldChange, decreasing = TRUE), ])
head(resSig_EC_DLN[order(resSig_EC_DLN$log2FoldChange, decreasing = FALSE), ])

write.csv(resSig_EC_DLN[order(resSig_EC_DLN$log2FoldChange, decreasing = TRUE), ],file = "DE_EC_vs_DLN_all.csv")
deg_EC_DLN=rownames(resSig_EC_DLN)

# Count the number of genes with significant differential expression at a FDR of 5%:
table(resSig_EC_DLN$padj < 0.05)



##### DE analysis of EC vs ST#####
#Trim sample table for EC vs ST 
keep_samples_EC=grep("EC", samples_hatchling$pond)
samples_EC=samples_hatchling[keep_samples_EC,]
keep_samples_ST=grep("ST", samples_hatchling$pond)
samples_ST=samples_hatchling[keep_samples_ST,]

samples_EC_ST=rbind(samples_EC,samples_ST)

#Trim counts table for DE analysis between EC vs ST 

EC_keep=grep("EC", colnames(counts_hatchling)) 
EC=dds_hatchling[,EC_keep]

ST_keep=grep("ST",colnames(counts_hatchling))
ST=dds_hatchling[,ST_keep]

EC_ST=cbind(EC@assays@data@listData$counts,ST@assays@data@listData$counts) #combine two count table and then make DEseqdataset from matrix
colData_EC_ST=data.frame(samples_EC_ST[,c("pond")])
rownames(colData_EC_ST)=samples_EC_ST$name
colnames(colData_EC_ST)="pond"

colnames(EC_ST)==rownames(colData_EC_ST)  # check if the count table has the same order as sample table

dds_EC_ST=DESeqDataSetFromMatrix(countData = EC_ST, 
                                 colData = colData_EC_ST,
                                 design = ~ pond) #making DESeqDataSet for EC and ST only


dds_EC_ST$pond = factor(dds_EC_ST$pond, levels = c("EC","ST"))

# Test for DE genes
dds_EC_ST_1=DESeq(dds_EC_ST)
res_EC_ST=results(dds_EC_ST_1, contrast = c("pond", "EC", "ST"))   #test DE for EC versus ST
head(res_EC_ST)


# Inspect the result tables of significantly upregulated and downregulated genes, at a 5% false discovery rate (FDR) as follows:
resSig_EC_ST=res_EC_ST[which(res_EC_ST$padj<0.05),]
head(resSig_EC_ST[order(resSig_EC_ST$log2FoldChange, decreasing = TRUE), ])
head(resSig_EC_ST[order(resSig_EC_ST$log2FoldChange, decreasing = FALSE), ])

write.csv(resSig_EC_ST[order(resSig_EC_ST$log2FoldChange, decreasing = TRUE), ],file = "DE_EC_vs_ST_all.csv")
deg_EC_ST=rownames(resSig_EC_ST)

# Count the number of genes with significant differential expression at a FDR of 5%:
table(resSig_EC_ST$padj < 0.05)


# DE analysis within roadside locations ------
##### DE analysis of BW vs RSB ####
#Trim sample table for BW vs RSB
keep_samples_BW=grep("BW", samples_hatchling$pond)
samples_BW=samples_hatchling[keep_samples_BW,]
keep_samples_RSB=grep("RSB", samples_hatchling$pond)
samples_RSB=samples_hatchling[keep_samples_RSB,]

samples_BW_RSB=rbind(samples_BW,samples_RSB)

#Trim counts table for DE analysis between BW vs RSB
BW_keep=grep("BW", colnames(counts_hatchling)) 
BW=dds_hatchling[,BW_keep]

RSB_keep=grep("RSB",colnames(counts_hatchling))
RSB=dds_hatchling[,RSB_keep]

BW_RSB=cbind(BW@assays@data@listData$counts,RSB@assays@data@listData$counts) #combine two count table and then make DEseqdataset from matrix
colData_BW_RSB=data.frame(samples_BW_RSB[,c("pond")])
rownames(colData_BW_RSB)=samples_BW_RSB$name
colnames(colData_BW_RSB)="pond"

colnames(BW_RSB)==rownames(colData_BW_RSB)  # check if the count table has the same order as sample table

dds_BW_RSB=DESeqDataSetFromMatrix(countData = BW_RSB, 
                                  colData = colData_BW_RSB,
                                  design = ~ pond) #making DESeqDataSet for BW and RSB only


dds_BW_RSB$pond = factor(dds_BW_RSB$pond, levels = c("BW","RSB"))

# Test for DE genes
dds_BW_RSB_1=DESeq(dds_BW_RSB)
res_BW_RSB=results(dds_BW_RSB_1, contrast = c("pond", "BW", "RSB"))   #test DE for BW versus RSB
head(res_BW_RSB)


# Inspect the result tables of significantly upregulated and downregulated genes, at a 5% false discovery rate (FDR) as follows:
resSig_BW_RSB=res_BW_RSB[which(res_BW_RSB$padj<0.05),]
head(resSig_BW_RSB[order(resSig_BW_RSB$log2FoldChange, decreasing = TRUE), ])
head(resSig_BW_RSB[order(resSig_BW_RSB$log2FoldChange, decreasing = FALSE), ])

write.csv(resSig_BW_RSB[order(resSig_BW_RSB$log2FoldChange, decreasing = TRUE), ],file = "DE_BW_vs_RSB_all.csv")
deg_BW_RSB=rownames(resSig_BW_RSB)

# Count the number of genes with significant differential expression at a FDR of 5%:
table(resSig_BW_RSB$padj < 0.05)


##### DE analysis of BW vs EC #####
#Trim sample table for BW vs EC
keep_samples_BW=grep("BW", samples_hatchling$pond)
samples_BW=samples_hatchling[keep_samples_BW,]
keep_samples_EC=grep("EC", samples_hatchling$pond)
samples_EC=samples_hatchling[keep_samples_EC,]

samples_BW_EC=rbind(samples_BW,samples_EC)

#Trim counts table for DE analysis between BW vs EC
BW_keep=grep("BW", colnames(counts_hatchling)) 
BW=dds_hatchling[,BW_keep]

EC_keep=grep("EC",colnames(counts_hatchling))
EC=dds_hatchling[,EC_keep]

BW_EC=cbind(BW@assays@data@listData$counts,EC@assays@data@listData$counts) #combine two count table and then make DEseqdataset from matrix
colData_BW_EC=data.frame(samples_BW_EC[,c("pond")])
rownames(colData_BW_EC)=samples_BW_EC$name
colnames(colData_BW_EC)="pond"

colnames(BW_EC)==rownames(colData_BW_EC)  # check if the count table has the same order as sample table

dds_BW_EC=DESeqDataSetFromMatrix(countData = BW_EC, 
                                 colData = colData_BW_EC,
                                 design = ~ pond) #making DESeqDataSet for BW and EC only


dds_BW_EC$pond = factor(dds_BW_EC$pond, levels = c("BW","EC"))

# Test for DE genes
dds_BW_EC_1=DESeq(dds_BW_EC)
res_BW_EC=results(dds_BW_EC_1, contrast = c("pond", "BW", "EC"))   #test DE for BW versus EC
head(res_BW_EC)


# Inspect the result tables of significantly upregulated and downregulated genes, at a 5% false discovery rate (FDR) as follows:
resSig_BW_EC=res_BW_EC[which(res_BW_EC$padj<0.05),]
head(resSig_BW_EC[order(resSig_BW_EC$log2FoldChange, decreasing = TRUE), ])
head(resSig_BW_EC[order(resSig_BW_EC$log2FoldChange, decreasing = FALSE), ])

write.csv(resSig_BW_EC[order(resSig_BW_EC$log2FoldChange, decreasing = TRUE), ],file = "DE_BW_vs_EC_all.csv")
deg_BW_EC=rownames(resSig_BW_EC)

# Count the number of genes with significant differential expression at a FDR of 5%:
table(resSig_BW_EC$padj < 0.05)


##### DE analysis of RSB vs EC #####
#Trim sample table for RSB vs EC
keep_samples_RSB=grep("RSB", samples_hatchling$pond)
samples_RSB=samples_hatchling[keep_samples_RSB,]
keep_samples_EC=grep("EC", samples_hatchling$pond)
samples_EC=samples_hatchling[keep_samples_EC,]

samples_RSB_EC=rbind(samples_RSB,samples_EC)

#Trim counts table for DE analysis between RSB vs EC
RSB_keep=grep("RSB", colnames(counts_hatchling)) 
RSB=dds_hatchling[,RSB_keep]

EC_keep=grep("EC",colnames(counts_hatchling))
EC=dds_hatchling[,EC_keep]

RSB_EC=cbind(RSB@assays@data@listData$counts,EC@assays@data@listData$counts) #combine two count table and then make DEseqdataset from matrix
colData_RSB_EC=data.frame(samples_RSB_EC[,c("pond")])
rownames(colData_RSB_EC)=samples_RSB_EC$name
colnames(colData_RSB_EC)="pond"

colnames(RSB_EC)==rownames(colData_RSB_EC)  # check if the count table has the same order as sample table

dds_RSB_EC=DESeqDataSetFromMatrix(countData = RSB_EC, 
                                  colData = colData_RSB_EC,
                                  design = ~ pond) #making DESeqDataSet for RSB and EC only


dds_RSB_EC$pond = factor(dds_RSB_EC$pond, levels = c("EC","RSB"))

# Test for DE genes
dds_RSB_EC_1=DESeq(dds_RSB_EC)
res_RSB_EC=results(dds_RSB_EC_1, contrast = c("pond", "RSB", "EC"))   #test DE for RSB versus EC
head(res_RSB_EC)


# Inspect the result tables of significantly upregulated and downregulated genes, at a 5% false discovery rate (FDR) as follows:
resSig_RSB_EC=res_RSB_EC[which(res_RSB_EC$padj<0.05),]
head(resSig_RSB_EC[order(resSig_RSB_EC$log2FoldChange, decreasing = TRUE), ])
head(resSig_RSB_EC[order(resSig_RSB_EC$log2FoldChange, decreasing = FALSE), ])

write.csv(resSig_RSB_EC[order(resSig_RSB_EC$log2FoldChange, decreasing = TRUE), ],file = "DE_RSB_vs_EC_all.csv")
deg_RSB_EC=rownames(resSig_RSB_EC)

# Count the number of genes with significant differential expression at a FDR of 5%:
table(resSig_RSB_EC$padj < 0.05)




# DE analysis within woodland locations ------
##### DE analysis of RO vs DLN #####
#Trim sample table for RO vs DLN
keep_samples_RO=grep("RO", samples_hatchling$pond)
samples_RO=samples_hatchling[keep_samples_RO,]
keep_samples_DLN=grep("DLN", samples_hatchling$pond)
samples_DLN=samples_hatchling[keep_samples_DLN,]

samples_RO_DLN=rbind(samples_RO,samples_DLN)

#Trim counts table for DE analysis between RO vs DLN
RO_keep=grep("RO", colnames(counts_hatchling)) 
RO=dds_hatchling[,RO_keep]

DLN_keep=grep("DLN",colnames(counts_hatchling))
DLN=dds_hatchling[,DLN_keep]

RO_DLN=cbind(RO@assays@data@listData$counts,DLN@assays@data@listData$counts) #combine two count table and then make DEseqdataset from matrix
colData_RO_DLN=data.frame(samples_RO_DLN[,c("pond")])
rownames(colData_RO_DLN)=samples_RO_DLN$name
colnames(colData_RO_DLN)="pond"

colnames(RO_DLN)==rownames(colData_RO_DLN)  # check if the count table has the same order as sample table

dds_RO_DLN=DESeqDataSetFromMatrix(countData = RO_DLN, 
                                  colData = colData_RO_DLN,
                                  design = ~ pond) #making DESeqDataSet for RO and DLN only


dds_RO_DLN$pond = factor(dds_RO_DLN$pond, levels = c("DLN","RO"))

# Test for DE genes
dds_RO_DLN_1=DESeq(dds_RO_DLN)
res_RO_DLN=results(dds_RO_DLN_1, contrast = c("pond", "RO", "DLN"))   #test DE for RO versus DLN
head(res_RO_DLN)


# Inspect the result tables of significantly upregulated and downregulated genes, at a 5% false discovery rate (FDR) as follows:
resSig_RO_DLN=res_RO_DLN[which(res_RO_DLN$padj<0.05),]
head(resSig_RO_DLN[order(resSig_RO_DLN$log2FoldChange, decreasing = TRUE), ])
head(resSig_RO_DLN[order(resSig_RO_DLN$log2FoldChange, decreasing = FALSE), ])

write.csv(resSig_RO_DLN[order(resSig_RO_DLN$log2FoldChange, decreasing = TRUE), ],file = "DE_RO_vs_DLN_all.csv")
deg_RO_DLN=rownames(resSig_RO_DLN)

# Count the number of genes with significant differential expression at a FDR of 5%:
table(resSig_RO_DLN$padj < 0.05)


##### DE analysis of RO vs ST #####
#Trim sample table for RO vs ST
keep_samples_RO=grep("RO", samples_hatchling$pond)
samples_RO=samples_hatchling[keep_samples_RO,]
keep_samples_ST=grep("ST", samples_hatchling$pond)
samples_ST=samples_hatchling[keep_samples_ST,]

samples_RO_ST=rbind(samples_RO,samples_ST)

#Trim counts table for DE analysis between RO vs ST
RO_keep=grep("RO", colnames(counts_hatchling)) 
RO=dds_hatchling[,RO_keep]

ST_keep=grep("ST",colnames(counts_hatchling))
ST=dds_hatchling[,ST_keep]

RO_ST=cbind(RO@assays@data@listData$counts,ST@assays@data@listData$counts) #combine two count table and then make DEseqdataset from matrix
colData_RO_ST=data.frame(samples_RO_ST[,c("pond")])
rownames(colData_RO_ST)=samples_RO_ST$name
colnames(colData_RO_ST)="pond"

colnames(RO_ST)==rownames(colData_RO_ST)  # check if the count table has the same order as sample table

dds_RO_ST=DESeqDataSetFromMatrix(countData = RO_ST, 
                                 colData = colData_RO_ST,
                                 design = ~ pond) #making DESeqDataSet for RO and ST only


dds_RO_ST$pond = factor(dds_RO_ST$pond, levels = c("RO","ST"))

# Test for DE genes
dds_RO_ST_1=DESeq(dds_RO_ST)
res_RO_ST=results(dds_RO_ST_1, contrast = c("pond", "RO", "ST"))   #test DE for RO versus ST
head(res_RO_ST)


# Inspect the result tables of significantly upregulated and downregulated genes, at a 5% false discovery rate (FDR) as follows:
resSig_RO_ST=res_RO_ST[which(res_RO_ST$padj<0.05),]
head(resSig_RO_ST[order(resSig_RO_ST$log2FoldChange, decreasing = TRUE), ])
head(resSig_RO_ST[order(resSig_RO_ST$log2FoldChange, decreasing = FALSE), ])

write.csv(resSig_RO_ST[order(resSig_RO_ST$log2FoldChange, decreasing = TRUE), ],file = "DE_RO_vs_ST_all.csv")
deg_RO_ST=rownames(resSig_RO_ST)

# Count the number of genes with significant differential expression at a FDR of 5%:
table(resSig_RO_ST$padj < 0.05)


##### DE analysis of DLN vs ST #####
#Trim sample table for DLN vs ST
keep_samples_DLN=grep("DLN", samples_hatchling$pond)
samples_DLN=samples_hatchling[keep_samples_DLN,]
keep_samples_ST=grep("ST", samples_hatchling$pond)
samples_ST=samples_hatchling[keep_samples_ST,]

samples_DLN_ST=rbind(samples_DLN,samples_ST)

#Trim counts table for DE analysis between DLN vs ST
DLN_keep=grep("DLN", colnames(counts_hatchling)) 
DLN=dds_hatchling[,DLN_keep]

ST_keep=grep("ST",colnames(counts_hatchling))
ST=dds_hatchling[,ST_keep]

DLN_ST=cbind(DLN@assays@data@listData$counts,ST@assays@data@listData$counts) #combine two count table and then make DEseqdataset from matrix
colData_DLN_ST=data.frame(samples_DLN_ST[,c("pond")])
rownames(colData_DLN_ST)=samples_DLN_ST$name
colnames(colData_DLN_ST)="pond"

colnames(DLN_ST)==rownames(colData_DLN_ST)  # check if the count table has the same order as sample table

dds_DLN_ST=DESeqDataSetFromMatrix(countData = DLN_ST, 
                                  colData = colData_DLN_ST,
                                  design = ~ pond) #making DESeqDataSet for DLN and ST only


dds_DLN_ST$pond = factor(dds_DLN_ST$pond, levels = c("DLN","ST"))

# Test for DE genes
dds_DLN_ST_1=DESeq(dds_DLN_ST)
res_DLN_ST=results(dds_DLN_ST_1, contrast = c("pond", "DLN", "ST"))   #test DE for DLN versus ST
head(res_DLN_ST)


# Inspect the result tables of significantly upregulated and downregulated genes, at a 5% false discovery rate (FDR) as follows:
resSig_DLN_ST=res_DLN_ST[which(res_DLN_ST$padj<0.05),]
head(resSig_DLN_ST[order(resSig_DLN_ST$log2FoldChange, decreasing = TRUE), ])
head(resSig_DLN_ST[order(resSig_DLN_ST$log2FoldChange, decreasing = FALSE), ])

write.csv(resSig_DLN_ST[order(resSig_DLN_ST$log2FoldChange, decreasing = TRUE), ],file = "DE_DLN_vs_ST_all.csv")
deg_DLN_ST=rownames(resSig_DLN_ST)

# Count the number of genes with significant differential expression at a FDR of 5%:
table(resSig_DLN_ST$padj < 0.05)


####parallel gene analysis#####
gene1=read.csv("DE_BW_vs_EC_all.csv",header=T)
gene2=read.csv("DE_BW_vs_RSB_all.csv",header=T)
gene3=read.csv("DE_DLN_vs_ST_all.csv",header=T)
gene4=read.csv("DE_RO_vs_DLN_all.csv",header=T)
gene5=read.csv("DE_RO_vs_ST_all.csv",header=T)
gene6=read.csv("DE_RSB_vs_EC_all.csv",header=T)

gene_within=rbind(rbind(rbind(rbind(rbind(gene1,gene2),gene3),gene4),gene5),gene6)
gene_within_unique=gene_within[!duplicated(gene_within$X),]                  
gene_within_unique1=gene_within_unique$X

gene7=read.csv("DE_BW_vs_DLN_all.csv",header=T)
gene8=read.csv("DE_BW_vs_RO_all.csv",header=T)
gene9=read.csv("DE_BW_vs_ST_all.csv",header=T)
gene10=read.csv("DE_EC_vs_DLN_all.csv",header=T)
gene11=read.csv("DE_EC_vs_RO_all.csv",header=T)
gene12=read.csv("DE_EC_vs_ST_all.csv",header=T)
gene13=read.csv("DE_RSB_vs_DLN_all.csv",header=T)
gene14=read.csv("DE_RSB_vs_RO_all.csv",header=T)
gene15=read.csv("DE_RSB_vs_ST_all.csv",header=T)

#########gene parallel########
gene7_within=gene7$X %in% gene_within_unique1
gene7_clean=gene7[!gene7_within,]

gene8_within=gene8$X %in% gene_within_unique1
gene8_clean=gene8[!gene8_within,]

gene9_within=gene9$X %in% gene_within_unique1
gene9_clean=gene9[!gene9_within,]

gene10_within=gene10$X %in% gene_within_unique1
gene10_clean=gene10[!gene10_within,]

gene11_within=gene11$X %in% gene_within_unique1
gene11_clean=gene11[!gene11_within,]

gene12_within=gene12$X %in% gene_within_unique1
gene12_clean=gene12[!gene12_within,]

gene13_within=gene13$X %in% gene_within_unique1
gene13_clean=gene13[!gene13_within,]

gene14_within=gene14$X %in% gene_within_unique1
gene14_clean=gene14[!gene14_within,]

gene15_within=gene15$X %in% gene_within_unique1
gene15_clean=gene15[!gene15_within,]


merge1=merge(gene7_clean,gene8_clean,by="X")
merge2=merge(merge1,gene9_clean,by="X")
merge3=merge(merge2,gene10_clean,by="X")
merge4=merge(merge3,gene11_clean,by="X")
merge5=merge(merge4,gene12_clean,by="X")
merge6=merge(merge5,gene13_clean,by="X")
merge7=merge(merge6,gene14_clean,by="X")
merge8=merge(merge7,gene15_clean,by="X")

write.csv(merge8,file = "deg_parallel.csv")