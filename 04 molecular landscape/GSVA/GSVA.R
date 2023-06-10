table(FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE)
# HER2_0      HER2_low HER2_positive 
# 91           434           182 

##GSVA  RNA_C1
her2_0_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_0"),"PatientCode"]
her2_low_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_low"),"PatientCode"]
her2_pos_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_positive"),"PatientCode"]

group_list<-c(rep("HER2_0",length(her2_0_list)),rep("HER2_low",length(her2_low_list)),rep("HER2_positive",length(her2_pos_list)))
names(group_list)<-c(her2_0_list,her2_low_list,her2_pos_list)
table(group_list)

group_list <- group_list[paste(names(group_list),"_RNA_T",sep = "") %in% colnames(FUSCC_HER2_low_project_RNA_seq_log2FPKM)]


exp<-FUSCC_HER2_low_project_RNA_seq_log2FPKM
exp <- exp[,paste(names(group_list),"_RNA_T",sep = "")]
exp <- as.matrix(exp )

c1_genesets<-clusterProfiler::read.gmt("c1.all.v7.5.1.symbols.gmt")
genesets = split(c1_genesets$gene, c1_genesets$term)

length(genesets)
#299

library(GSVA)

t1 <- Sys.time()
esmicro <-gsva(expr=exp,gset.idx.list=genesets, min.sz=5, max.sz=500,method="gsva",kcdf='Gaussian') 
t2 <- Sys.time()
t2-t1

save(file = "RNA_C1.Rdata",esmicro)



##GSVA SSgsea RNA_c2
her2_0_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_0"),"PatientCode"]
her2_low_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_low"),"PatientCode"]
her2_pos_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_positive"),"PatientCode"]

group_list<-c(rep("HER2_0",length(her2_0_list)),rep("HER2_low",length(her2_low_list)),rep("HER2_positive",length(her2_pos_list)))
names(group_list)<-c(her2_0_list,her2_low_list,her2_pos_list)
table(group_list)

group_list <- group_list[paste(names(group_list),"_RNA_T",sep = "") %in% colnames(FUSCC_HER2_low_project_RNA_seq_log2FPKM)]


exp<-FUSCC_HER2_low_project_RNA_seq_log2FPKM
exp <- exp[,paste(names(group_list),"_RNA_T",sep = "")]
exp <- as.matrix(exp )

c2_genesets<-clusterProfiler::read.gmt("c2.all.v7.5.symbols.gmt")
genesets = split(c2_genesets$gene, c2_genesets$term)
genesets <- genesets[substr(names(genesets),1,8) == "REACTOME"]

length(genesets)
#1615

library(GSVA)

t1 <- Sys.time()
esmicro <-gsva(expr=exp,gset.idx.list=genesets, min.sz=5, max.sz=500,method="gsva",kcdf='Gaussian') 
t2 <- Sys.time()
t2-t1

save(file = "RNA_C2.Rdata",esmicro)
