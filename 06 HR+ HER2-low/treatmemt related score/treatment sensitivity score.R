table(FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE)
# HER2_0      HER2_low HER2_positive 
# 91           434           182 



###################
##Endocrine sensitivity score
######################
endocrime_resistence_score <- as.numeric((0.8*FUSCC_HER2_low_project_RNA_seq_log2FPKM["ESR1",]+1.2*FUSCC_HER2_low_project_RNA_seq_log2FPKM["PGR",]+FUSCC_HER2_low_project_RNA_seq_log2FPKM["BCL2",]+FUSCC_HER2_low_project_RNA_seq_log2FPKM["SCUBE2",])/4)
names(endocrime_resistence_score) <- colnames(FUSCC_HER2_low_project_RNA_seq_log2FPKM)

her2_0_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_0")  & FUSCC_HER2_low_project_Cohort.Info$Adjuvant_endocrine_therapy %in% c("Yes","YES")& FUSCC_HER2_low_project_Cohort.Info$Clinical_Subtype %in% c("HR+HER2-"),"PatientCode"]
her2_low_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_low") & FUSCC_HER2_low_project_Cohort.Info$Adjuvant_endocrine_therapy %in% c("Yes","YES")& FUSCC_HER2_low_project_Cohort.Info$Clinical_Subtype %in% c("HR+HER2-"),"PatientCode"]

group_list<-c(rep("HER2-0",length(her2_0_list)),rep("HER2-low",length(her2_low_list)))
names(group_list)<-c(her2_0_list,her2_low_list)
table(group_list)
# HER2-0 HER2-low 
# 57      305

my_comparisons <- list( c("HER2-0", "HER2-low"))

dataset=data.frame(Endocrine_sensitivity_score=endocrime_resistence_score[paste(names(group_list),"_RNA_T",sep = "")],group=group_list)
write.csv(file = "Endocrine_sensitivity_score.csv",dataset,quote = F)
library(ggplot2)
library(ggpubr)
#dataset<-na.omit(dataset)
pdf(file="Endocrine_sensitivity_score_her2_0_low.pdf",width = 2,height = 5)
ggboxplot(dataset, x = "group", y = "Endocrine_sensitivity_score",
          color = "group", palette = "aaas",
          add = "jitter",title = "Endocrine sensitivity score",ylab = "Endocrine sensitivity score")+stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")
dev.off()

wilcox.test(dataset[dataset$group == "HER2-0","Endocrine_sensitivity_score"],dataset[dataset$group == "HER2-low","Endocrine_sensitivity_score"])

###################
##setep
######################

setep <- read.csv(file="SETep.csv")
rownames(setep) <- setep[,1]

her2_0_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_0")  & FUSCC_HER2_low_project_Cohort.Info$Adjuvant_endocrine_therapy %in% c("Yes","YES")& FUSCC_HER2_low_project_Cohort.Info$Clinical_Subtype %in% c("HR+HER2-"),"PatientCode"]
her2_low_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_low") & FUSCC_HER2_low_project_Cohort.Info$Adjuvant_endocrine_therapy %in% c("Yes","YES")& FUSCC_HER2_low_project_Cohort.Info$Clinical_Subtype %in% c("HR+HER2-"),"PatientCode"]

group_list<-c(rep("HER2-0",length(her2_0_list)),rep("HER2-low",length(her2_low_list)))
names(group_list)<-c(her2_0_list,her2_low_list)
table(group_list)
# HER2-0 HER2-low 
# 57      305

my_comparisons <- list( c("HER2-0", "HER2-low"))

dataset=data.frame(SETep=setep[names(group_list),"SETep"],group=group_list)

library(ggplot2)
library(ggpubr)

pdf(file="SETep_her2_0_low.pdf",width = 2,height = 5)
ggboxplot(dataset, x = "group", y = "SETep",
          color = "group", palette = "aaas",
          add = "jitter",title = "SETep",ylab = "SETep")+stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")
dev.off()

wilcox.test(dataset[dataset$group == "HER2-0","SETep"],dataset[dataset$group == "HER2-low","SETep"])

