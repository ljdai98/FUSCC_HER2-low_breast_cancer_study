table(FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE)
# HER2_0      HER2_low HER2_positive 
# 91           434           182 

load("../../04 molecular landscape/GSVA/RNA_C1.Rdata")
esmicro <- data.frame(esmicro)

her2_0_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_0")  & FUSCC_HER2_low_project_Cohort.Info$Clinical_Subtype %in% c("HR+HER2-"),"PatientCode"]
her2_low_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_low") & FUSCC_HER2_low_project_Cohort.Info$Clinical_Subtype %in% c("HR+HER2-"),"PatientCode"]

group_list<-c(rep("HER2-0",length(her2_0_list)),rep("HER2-low",length(her2_low_list)))
names(group_list)<-c(her2_0_list,her2_low_list)
table(group_list)

group_list <- group_list[paste(names(group_list),"_RNA_T",sep = "") %in% colnames(esmicro)]
names(group_list) <- paste(names(group_list),"_RNA_T",sep = "")

my_comparisons <- list( c("HER2-0", "HER2-low"))

library(ggplot2)
library(ggpubr)

pdf("GSVA_C1_HR_pos_HER2_0_vs_HER2_low.pdf",width = 2.5,height = 11)

gene <- "chr17q21"
dataset=data.frame(mRNA=as.numeric(esmicro[gene,names(group_list)]),group=group_list)
p1 <- ggboxplot(dataset, x = "group", y = "mRNA",
                color = "group", palette = "aaas",
                add = "jitter",title = gene,ylab ="log2(FPKM+1)")+stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")


gene <- "chr17q11"
dataset=data.frame(mRNA=as.numeric(esmicro[gene,names(group_list)]),group=group_list)
p2 <- ggboxplot(dataset, x = "group", y = "mRNA",
                color = "group", palette = "aaas",
                add = "jitter",title = gene,ylab ="log2(FPKM+1)")+stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")


gene <- "chr17q12"
dataset=data.frame(mRNA=as.numeric(esmicro[gene,names(group_list)]),group=group_list)
p3 <- ggboxplot(dataset, x = "group", y = "mRNA",
                color = "group", palette = "aaas",
                add = "jitter",title = gene,ylab ="log2(FPKM+1)")+stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")
library(patchwork)
print(p1/p2/p3)

dev.off()


