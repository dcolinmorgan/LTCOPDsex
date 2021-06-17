library(qsmooth)
library(tidyverse)
library(GEOquery)
library(limma)
library(umap)
library(data.table)
library(sva)

##LTCOPD
gset <- getGEO("GSE76925", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL10558", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

load("/proj/regeps/regep00/studies/LTCOPD/analyses/remdj/expression/analysis/expSet.clean_V14_1435701821.RData")
exprs(gset) <- normalizeBetweenArrays(exprs(gset)) # normalize data
LTCOPD_and_ECLIPSE_Generic_Plate_Report_Fixed_062514 <- read_delim("/proj/regeps/regep00/studies/LTCOPD/analyses/remdj/paperMethyLTCOPD_ChanningReview/LTCOPD and ECLIPSE Generic Plate Report Fixed_062514.csv",delim='\t')

j<-LTCOPD_and_ECLIPSE_Generic_Plate_Report_Fixed_062514[LTCOPD_and_ECLIPSE_Generic_Plate_Report_Fixed_062514$COLLECTION=='Lung Tissue',]
t<-j[,c("SUBJECT","PLATE ID","GENDER")]
jeff<-t[t$SUBJECT %in% colnames(exp.clean),]
dd<-unique(jeff)
LTCOPD<-(exprs(gset))
# LTCOPD<-DGEList(LTCOPD)

D2 = qsmooth(object=(LTCOPD), group_factor=dd$`GENDER`)

combat_edata3 = ComBat(dat=LTCOPD, batch=dd$`PLATE ID`)

D3 = ComBat(dat=D2, batch=dd$`PLATE ID`)

rownames(combat_edata3)<-gset@featureData@data$`Gene symbol`
colnames(combat_edata3)<-paste(gset@phenoData@data$`age:ch1`,'_',gset@phenoData@data$`Sex:ch1`,'_',gset@phenoData@data$`copd:ch1`,sep='')
D3<-data.frame(D2@qsmoothData)
rownames(D3)<-gset@featureData@data$`Gene symbol`
colnames(D3)<-paste(gset@phenoData@data$`age:ch1`,'_',gset@phenoData@data$`Sex:ch1`,'_',gset@phenoData@data$`copd:ch1`,sep='')


# write.table(LTCOPD,'LTCOPD_raw.txt',sep='\t')
write.table(round(combat_edata3,3),'LTCOPD_combat.txt',sep='\t')
write.table(round(D3,3),'LTCOPD_qsmooth_combat.txt',sep='\t')
# write.table(gset@featureData@data$`Gene symbol`,'LTCOPD_genes.txt',sep='\t')

ff<-data.frame(gset$`Sex:ch1`)
ff$race<-gset$`race:ch1`
ff$sex<-gset$`Sex:ch1`
ff$py<-gset$`packyears:ch1`
ff$copd<-gset$`copd:ch1`
ff$age<-round(as.numeric(gset$`age:ch1`))

mmdata<-combat_edata3[,ff$sex=="M"]
rownames(mmdata)<-gset@featureData@data$`Gene symbol`
ffdata<-combat_edata3[,ff$sex=="F"]
rownames(ffdata)<-gset@featureData@data$`Gene symbol`

# table(ff$sex,ff$copd)

mm<-ff[ff$sex=="M",]
ff<-ff[ff$sex=="F",]

designA <- model.matrix(~copd+age+py,ff)

fit <- lmFit(ffdata, designA)  # fit linear model
fit2A <- eBayes(fit, 0.01)

tT2<-topTable(fit2A, coef=2,adjust="fdr", sort.by="p", number=1000) ## plate batch variableas surrogate ??
tT2 <- subset(tT2, select=c("ID",'AveExpr','logFC',"P.Value",'adj.P.Val')) ##
write.table(tT2, file=paste('analyses/F_LTCOPD_de_','061521','.txt',sep=''), row.names=T, sep="\t")


designB <- model.matrix(~copd+age+py,mm)

fit <- lmFit(mmdata, designB)  # fit linear model
fitA <- eBayes(fit, 0.01)

tT<-topTable(fitA, coef=2,adjust="fdr", sort.by="p", number=1000) ## plate batch variableas surrogate ??
tT <- subset(tT, select=c("ID",'AveExpr','logFC',"P.Value",'adj.P.Val')) ##
write.table(tT, file=paste('analyses/M_LTCOPD_de_','061521','.txt',sep=''), row.names=T, sep="\t")





##LTRC
counts <- fread(paste0("/proj/regeps/regep00/studies/COPDGene/analyses/reagh/LTRC_ILD/RNA/LTRC_topmed_to3_rnaseq_1.rsem_genes_expected_count.txt"),header = T)
counts <- counts %>% column_to_rownames(var = "gene_id")

analysis.final.modCopd <- fread(paste0('/proj/edith/regeps/regep00/studies/COPDGene/analyses/reagh/LTRC_ILD/modCopd_pathConservRna.pheno_022221.csv'),colClasses = list(character=1))

### Factor binary variables

RNApheno.modCopd <- analysis.final.modCopd %>% mutate(modCopd.control=as.factor(modCopd_pathConserv),
                                                      gender=as.factor(gender),
                                                      race=as.factor(race),
                                                      smoking_current=as.factor(smoking_current),
                                                      batch=as.factor(batch))

table(RNApheno.modCopd$modCopd.control, exclude =F)


## Read in and prepare count data

## Clean and prep counts file  
### Limit samples using rna bam ids from phenotype file
rnaBamId.modCopd <- as.character(RNApheno.modCopd$rnaBamId)
counts.clean <- counts[,rnaBamId.modCopd]
counts.clean <- counts.clean[!grepl("ERCC", rownames(counts.clean)),]

# counts.clean3<-counts.clean[RNApheno.modCopd$gender=="0"]
### Make DGEList object 
counts.modCopd <- qsmooth(object=counts.clean,group_factor=RNApheno.modCopd$gender)

dim(counts.modCopd)

### Remove genes with low expression in >50% of subjects
# cpms = cpm(counts.modCopd)
# keep.modCopd = rowSums(cpms >1) >=322
# counts.modCopd = counts.modCopd[keep.modCopd,]
data1<-counts.modCopd@qsmoothData
colnames(data1)<-paste(RNApheno.modCopd$age,'_',RNApheno.modCopd$gender,'_',RNApheno.modCopd$modCopd.control,sep='')
combat_edata1 = ComBat(dat=data1, batch=RNApheno.modCopd$batch)
# combat_edata4 = ComBat(dat=counts.clean, batch=RNApheno.modCopd$batch)
colnames(combat_edata1)<-paste(RNApheno.modCopd$age,'_',RNApheno.modCopd$gender,'_',RNApheno.modCopd$modCopd.control,sep='')

hist(combat_edata1)# jj<-(str_split(rownames(combat_edata1),'[.]'))
# write.csv(jj,"~/analyses/LTRC/ilmn.txt",sep='\t')

# write.csv(combat_edata1,"~/analyses/LTRC/F_results.modCopd_060521.csv",sep=',', row.names = T)
# 
# write.csv(counts.modCopd$counts,"~/analyses/LTRC/modCopd_060521.csv",sep=',', row.names = T)



write.table(round(combat_edata1,3),'LTRC_qsmooth_combat.txt',sep='\t')
write.table(round(data1,3),'LTRC_qsmooth.txt',sep='\t')
# write.table(counts.clean,'LTRC_raw.txt',sep='\t')



ff<-data.frame(RNApheno.modCopd$gender)
ff$race<-RNApheno.modCopd$race
ff$sex<-RNApheno.modCopd$gender
ff$py<-round(as.numeric(RNApheno.modCopd$smoking_packyears))
ff$copd<-RNApheno.modCopd$modCopd.control
ff$age<-round(as.numeric(RNApheno.modCopd$age))

mmdata<-combat_edata1[,ff$sex=="1"]
# rownames(mmdata)<-gset@featureData@data$`Gene symbol`
ffdata<-combat_edata1[,ff$sex=="0"]
# rownames(ffdata)<-gset@featureData@data$`Gene symbol`

# table(ff$sex,ff$copd)

mm<-ff[ff$sex=="1",]
ff<-ff[ff$sex=="0",]

designA <- model.matrix(~copd+age+py,ff)

fit <- lmFit(ffdata, designA)  # fit linear model
fit2A <- eBayes(fit, 0.01)

tT2<-topTable(fit2A, coef=2,adjust="fdr", sort.by="p", number=1000) ## plate batch variableas surrogate ??
tT2 <- subset(tT2, select=c('AveExpr','logFC',"P.Value",'adj.P.Val')) ##
write.table(tT2, file=paste('analyses/F_LTRC_de_','061521','.txt',sep=''), row.names=T, sep="\t")


designB <- model.matrix(~copd+age+py,mm)

fit <- lmFit(mmdata, designB)  # fit linear model
fitA <- eBayes(fit, 0.01)

tT<-topTable(fitA, coef=2,adjust="fdr", sort.by="p", number=1000) ## plate batch variableas surrogate ??
tT <- subset(tT, select=c('AveExpr','logFC',"P.Value",'adj.P.Val')) ##
write.table(tT, file=paste('analyses/M_LTRC_de_','061521','.txt',sep=''), row.names=T, sep="\t")





