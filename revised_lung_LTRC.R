library(tidyverse)
library(limma)

# bigsheet<-read.table('/proj/regeps/regep00/studies/LTRC/data/phenotype/data/freezes/20210310/ltrcLongTopMedHarm_20210310.csv',sep=',',header=TRUE)

# lung<-readRDS('~/../resiq/meffil/QC/LTRC.lung/old/minfi/results/ltrc.lung.beta.rcp')
lung<-readRDS('~/../resiq/meffil/QC/LTRC.lung/qc.clean/ltrc.lung.beta.rcp')

annot <- read.table("~/analyses/LTRC/min_850kannot.txt",sep='\t',header=FALSE)
colnames(annot)<-c('ensemble','chr','AddressA','GencodeBasicV12_NAME','Relation_to_Island')
# meta2 <- read_csv("/proj/regeps/regep00/studies/LTRC/data/phenotype/data/freezes/20210310/ltrcSamplePhenoIdMap_20210310.csv")
meta2 <- read_csv("/proj/regeps/regep00/studies/LTRC/data/phenotype/data/freezes/20210428/ltrcSamplePhenoIdMap_20210428.csv")

meta <- read.table("analyses/LTRC/ltrc_group_patid_LongFull_merge_phenoMap.csv",sep='\t',header=TRUE)
meta<-meta[meta$Project=='Methylation.lung',]
meta<-meta[!meta$clinCopd==.5,]   ## plot out phenotype variables for repeat visit subjects >> 
                                  ## pull out repeats, characterize lung sample location (upper and lower) >> 
                                  ## check dates, check diagnosis codes >> regional effect included into model??
                                  ## copd upper lung bias due to resectioning / volume reduction
metaA<-meta[meta$IPF.path==1,]
metaB<-meta[meta$IPF.path==0,]
meta<-rbind(metaA,metaB)
meta$smoking_packyears[is.na(meta$smoking_packyears)]<-mean(meta$smoking_packyears,na.rm=TRUE)

jeff2<-meta2 %>% subset(Project=='Methylation.lung')
# jeff2<-jeff2 %>% subset(Qc.Status=='Pass')
jeff2$patid<-as.numeric(jeff2$patid)

# jeff3<-jeff2[str_which(jeff2$lungLobe,'R_Lower'),]
# jeff4<-jeff2[str_which(jeff2$lungLobe,'L_Upper'),]

meta<-merge(jeff2,meta,on='patid')
meta<-meta[!duplicated(meta$patid), ]

# meta[is.na(meta$lungLobe)] <- 0
meta<-meta[meta$lungLobe!=0,]
# meta<-drop_na(meta,'lungLobe') ## drop subjects without lung location

meta$lungLobe<-replace_na(meta$lungLobe,'0')
RNApheno.modCopd <- meta %>% mutate(modCopd.control=as.factor(clinCopd),
                                    ipf.control=as.factor(IPF.path),
                                    lungLobe=(lungLobe),
                                                      gender=as.factor(gender),
                                                      race=(race), ## make sure this going to factor in design matrix
                                                      smoking_current=as.factor((smoking_packyears))) ## this seems to drop out later in code

table(RNApheno.modCopd$modCopd.control, exclude =F)

RNApheno.modCopd$patid<-as.numeric(RNApheno.modCopd$patid)
# pheno_meta$patid<-pheno_meta$ALIAS

# RNApheno.modCopd<-merge(pheno_meta,RNApheno.modCopd,on=patid)
# date<-date()
# date<-str_replace_all(string=date, pattern=" ", repl="")
date<-Sys.Date()
## Read in and prepare count data

## Clean and prep counts file 
### Limit samples using rna bam ids from phenotype file
rnaBamId.modCopd <- as.character(RNApheno.modCopd$topmedId)
ww<-str_split(colnames(lung),'_')#[[1:length(colnames(lung))]][1]
qq<-data.frame(ww)
qq<-data.frame(t(qq))
colnames(lung)<-qq$X1
col.num <- which(colnames(lung) %in% RNApheno.modCopd$topmedId)
counts.clean <- lung[,col.num]



col.num <- which(RNApheno.modCopd$topmedId %in% colnames(counts.clean))
meta <- RNApheno.modCopd[col.num,]
meta$bmi<-NULL
meta$ht_cm<-NULL
meta$wt_kg<-NULL
meta<-unique(meta)




design.modCopd2 <- data.frame(model.matrix(~modCopd.control+lungLobe+age+race+gender+smoking_packyears, data=meta))
fit.modCopd <- lmFit(counts.clean, design.modCopd2)#sv)
fit.modCopd <- eBayes(fit.modCopd) ## check into which fdr >> strange all adj pval behavior 

results.modCopd<-topTable(fit.modCopd, coef=2, number=Inf) ## make sure this isnt (intercept) ## hardcode this aas "copdcont"
results.modCopd <- results.modCopd %>% add_rownames(var = 'ensemble') ## threshold of meaningful difference
results.modCopd$ensemble <- gsub("\\..*","",results.modCopd$ensemble) ##
results0.001FDR.modCopd <- subset(results.modCopd, results.modCopd$adj.P.Val < 0.05,)
results0.001FDR.modCopd<-merge(annot,results0.001FDR.modCopd,by='ensemble')
write.csv(results0.001FDR.modCopd,paste("~/analyses/LTRC/bench/lobe_lung_methyl_p001_",date,"_allcopd_cov.txt",sep=''), row.names = F)

design.modCopd2 <- data.frame(model.matrix(~ipf.control+lungLobe+age+race+gender+smoking_packyears, data=meta))
fit.modCopd <- lmFit(counts.clean, design.modCopd2)#sv)
fit.modCopd <- eBayes(fit.modCopd)
results.modCopd<-topTable(fit.modCopd, coef=2, number=Inf) ## make sure this isnt (intercept) ## hardcode this
results.modCopd <- results.modCopd %>% add_rownames(var = 'ensemble') ## threshold of meaningful difference
results.modCopd$ensemble <- gsub("\\..*","",results.modCopd$ensemble) ##
results0.001FDR.modCopd <- subset(results.modCopd, results.modCopd$adj.P.Val < 0.05,)
results0.001FDR.modCopd<-merge(annot,results0.001FDR.modCopd,by='ensemble')
write.csv(results0.001FDR.modCopd,paste("~/analyses/LTRC/bench/lobe_lung_methyl_p001_",date,"_allipf_cov.txt",sep=''), row.names = F)

design.modCopd2 <- data.frame(model.matrix(~modCopd.control, data=meta))
fit.modCopd <- lmFit(counts.clean, design.modCopd2)#sv)
fit.modCopd <- eBayes(fit.modCopd)
results.modCopd<-topTable(fit.modCopd, coef=2, number=Inf) ## make sure this isnt (intercept) ## hardcode this
results.modCopd <- results.modCopd %>% add_rownames(var = 'ensemble') ## threshold of meaningful difference
results.modCopd$ensemble <- gsub("\\..*","",results.modCopd$ensemble) ##
results0.001FDR.modCopd <- subset(results.modCopd, results.modCopd$adj.P.Val < 0.05,)
results0.001FDR.modCopd<-merge(annot,results0.001FDR.modCopd,by='ensemble')
write.csv(results0.001FDR.modCopd,paste("~/analyses/LTRC/bench/lobe_lung_methyl_p001_",date,"_copd_cov.txt",sep=''), row.names = F)

design.modCopd2 <- data.frame(model.matrix(~ipf.control, data=meta))
fit.modCopd <- lmFit(counts.clean, design.modCopd2)#sv)
fit.modCopd <- eBayes(fit.modCopd)
results.modCopd<-topTable(fit.modCopd, coef=2, number=Inf) ## make sure this isnt (intercept) ## hardcode this
results.modCopd <- results.modCopd %>% add_rownames(var = 'ensemble') ## threshold of meaningful difference
results.modCopd$ensemble <- gsub("\\..*","",results.modCopd$ensemble) ##
results0.001FDR.modCopd <- subset(results.modCopd, results.modCopd$adj.P.Val < 0.05,)
results0.001FDR.modCopd<-merge(annot,results0.001FDR.modCopd,by='ensemble')
write.csv(results0.001FDR.modCopd,paste("~/analyses/LTRC/bench/lobe_lung_methyl_p001_",date,"_ipf_cov.txt",sep=''), row.names = F)

design.modCopd3 <- data.frame(model.matrix(~age, data=meta))
fit.modCopd <- lmFit(counts.clean, design.modCopd3)#sv)
fit.modCopd <- eBayes(fit.modCopd)
results.modCopd<-topTable(fit.modCopd, coef=2, number=Inf) ## make sure this isnt (intercept) ## hardcode this
results.modCopd <- results.modCopd %>% add_rownames(var = 'ensemble') ## threshold of meaningful difference
results.modCopd$ensemble <- gsub("\\..*","",results.modCopd$ensemble) ##
results0.001FDR.modCopd <- subset(results.modCopd, results.modCopd$adj.P.Val < 0.05,)
results0.001FDR.modCopd<-merge(annot,results0.001FDR.modCopd,by='ensemble')
write.csv(results0.001FDR.modCopd,paste("~/analyses/LTRC/bench/love_lung_methyl_p001_",date,"_age.txt",sep=''), row.names = F)

design.modCopd4 <- data.frame(model.matrix(~race, data=meta))
fit.modCopd <- lmFit(counts.clean, design.modCopd4)#sv)
fit.modCopd <- eBayes(fit.modCopd)
results.modCopd<-topTable(fit.modCopd, coef=2, number=Inf) ## make sure this isnt (intercept) ## hardcode this
results.modCopd <- results.modCopd %>% add_rownames(var = 'ensemble') ## threshold of meaningful difference
results.modCopd$ensemble <- gsub("\\..*","",results.modCopd$ensemble) ##
results0.001FDR.modCopd <- subset(results.modCopd, results.modCopd$adj.P.Val < 0.05,)
results0.001FDR.modCopd<-merge(annot,results0.001FDR.modCopd,by='ensemble')
write.csv(results0.001FDR.modCopd,paste("~/analyses/LTRC/bench/lobe_lung_methyl_p001_",date,"_race.txt",sep=''), row.names = F)

design.modCopd5 <- data.frame(model.matrix(~gender, data=meta))
fit.modCopd <- lmFit(counts.clean, design.modCopd5)#sv)
fit.modCopd <- eBayes(fit.modCopd)
results.modCopd<-topTable(fit.modCopd, coef=2, number=Inf) ## make sure this isnt (intercept) ## hardcode this
results.modCopd <- results.modCopd %>% add_rownames(var = 'ensemble') ## threshold of meaningful difference
results.modCopd$ensemble <- gsub("\\..*","",results.modCopd$ensemble) ##
results0.001FDR.modCopd <- subset(results.modCopd, results.modCopd$adj.P.Val < 0.05,)
results0.001FDR.modCopd<-merge(annot,results0.001FDR.modCopd,by='ensemble')
write.csv(results0.001FDR.modCopd,paste("~/analyses/LTRC/bench/love_lung_methyl_p001_",date,"_gender.txt",sep=''), row.names = F)


design.modCopd6 <- data.frame(model.matrix(~smoking_packyears, data=meta))
fit.modCopd <- lmFit(counts.clean, design.modCopd6)#sv)
fit.modCopd <- eBayes(fit.modCopd)
results.modCopd<-topTable(fit.modCopd, coef=2, number=Inf) ## make sure this isnt (intercept) ## hardcode this
results.modCopd <- results.modCopd %>% add_rownames(var = 'ensemble') ## threshold of meaningful difference
results.modCopd$ensemble <- gsub("\\..*","",results.modCopd$ensemble) ##
results0.001FDR.modCopd <- subset(results.modCopd, results.modCopd$adj.P.Val < 0.05,)
results0.001FDR.modCopd<-merge(annot,results0.001FDR.modCopd,by='ensemble')
write.csv(results0.001FDR.modCopd,paste("~/analyses/LTRC/bench/lobe_lung_methyl_p001_",date,"_PY.txt",sep=''), row.names = F)

design.modCopd7 <- data.frame(model.matrix(~lungLobe, data=meta))
fit.modCopd <- lmFit(counts.clean, design.modCopd7)#sv)
fit.modCopd <- eBayes(fit.modCopd)
results.modCopd<-topTable(fit.modCopd, coef=2, number=Inf) ## make sure this isnt (intercept) ## hardcode this
results.modCopd <- results.modCopd %>% add_rownames(var = 'ensemble') ## threshold of meaningful difference
results.modCopd$ensemble <- gsub("\\..*","",results.modCopd$ensemble) ##
results0.001FDR.modCopd <- subset(results.modCopd, results.modCopd$adj.P.Val < 0.05,)
results0.001FDR.modCopd<-merge(annot,results0.001FDR.modCopd,by='ensemble')
write.csv(results0.001FDR.modCopd,paste("~/analyses/LTRC/bench/lobe_lung_methyl_p001_",date,"_lobe.txt",sep=''), row.names = F)

table(meta$lungLobe)



##############
##############

## Create design matrices
design.modCopd <- model.matrix(~modCopd.control+lungLobe+ipf.control+age+race+gender+smoking_packyears, data=meta)

# design.modCopd <- model.matrix(~ipf.clinpath, data=meta)

jj<-colnames(design.modCopd)
jj<-jj[2:8]

AA<-t(combn(jj,7))
AA<-rbind(AA,cbind(t(combn(jj,6)),0))
AA<-rbind(AA,cbind(t(combn(jj,5)),0,0))
AA<-rbind(AA,cbind(t(combn(jj,4)),0,0,0))
AA<-rbind(AA,cbind(t(combn(jj,3)),0,0,0,0))
AA<-rbind(AA,cbind(t(combn(jj,2)),0,0,0,0,0))
AA<-rbind(AA,cbind(t(combn(jj,1)),0,0,0,0,0,0))
# AA<-rbind(AA,cbind(t(combn(jj,3)),0,0,0,0,0,0))
# AA<-rbind(AA,cbind(t(combn(jj,2)),0,0,0,0,0,0,0))
# AA<-rbind(AA,cbind(t(combn(jj,1)),0,0,0,0,0,0,0,0))
AA<-data.frame(AA)
AA$Z<-AA %>% unite("Z", X1:X2:X3:X4:X5:X6:X7, na.rm = TRUE, remove = TRUE,sep='+')

for(i in 1:nrow(AA)) {
  design.modCopd2<-data.frame(design.modCopd)
  row <- AA[i,]
  assign("model",row$Z)
  # do stuff with row
  
  design.modCopd2 <- model.matrix(as.formula(paste("~", paste(model),sep='')), data=design.modCopd2)
  print(paste('fitting',i))
  fit.modCopd <- lmFit(counts.clean, design.modCopd2)#sv)
  print(paste('correcting',i))
  fit.modCopd <- eBayes(fit.modCopd)
  
  ## Examine results
  results.modCopd<-topTable(fit.modCopd, coef=2, number=Inf) ## make sure this isnt (intercept) ## hardcode this
  results.modCopd <- results.modCopd %>% add_rownames(var = 'ensemble') ## threshold of meaningful difference
  results.modCopd$ensemble <- gsub("\\..*","",results.modCopd$ensemble) ##
  results0.001FDR.modCopd <- subset(results.modCopd, results.modCopd$adj.P.Val < 0.05,)
  results0.001FDR.modCopd<-merge(annot,results0.001FDR.modCopd,by='ensemble')
  
  out<-paste("~/analyses/LTRC/bench/lung_methyl_p001_",date,"_",paste(row$Z),".txt",sep='')
  print(paste('writing',i))
  write.table(results0.001FDR.modCopd,out,sep='\t', row.names = F)
  
}

