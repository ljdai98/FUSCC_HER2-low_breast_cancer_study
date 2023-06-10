table(FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE)
# HER2_0      HER2_low HER2_positive 
# 91           434           182 


her2_0_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_0")  & FUSCC_HER2_low_project_Cohort.Info$Clinical_Subtype %in% c("HR+HER2-"),"PatientCode"]
her2_low_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_low") & FUSCC_HER2_low_project_Cohort.Info$Clinical_Subtype %in% c("HR+HER2-"),"PatientCode"]


group_list<-c(rep("HER2-0",length(her2_0_list)),rep("HER2-low",length(her2_low_list)))
names(group_list)<-c(her2_0_list,her2_low_list)
table(group_list)

group_list <- group_list[paste(names(group_list),"_RNA_T",sep = "") %in% colnames(FUSCC_HER2_low_project_RNA_seq_log2FPKM)]
names(group_list) <- paste(names(group_list),"_RNA_T",sep = "")


my_comparisons <- list( c("HER2-0", "HER2-low"))


genes <- c("ESR1","BCL2","PGR","GPR160","FOXA1","BAG1","AR")


library(ggplot2)
library(ggpubr)

pdf(paste("luminal related genes_RNA_0-low(HR+).pdf",sep=""),width = 2,height = 5)
for (gene in genes){
  dataset=data.frame(mRNA=as.numeric(FUSCC_HER2_low_project_RNA_seq_log2FPKM[gene,names(group_list)]),group=group_list)
  p1 <- ggboxplot(dataset, x = "group", y = "mRNA",
                  color = "group", palette = "aaas",
                  add = "jitter",title = gene,ylab ="log2(FPKM+1)")+stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")
  print(p1)
}
dev.off()


for (gene in genes){
  dataset=data.frame(mRNA=as.numeric(FUSCC_HER2_low_project_RNA_seq_log2FPKM[gene,names(group_list)]),group=group_list)
  p1 <- wilcox.test(dataset[dataset$group == "HER2-0","mRNA"],dataset[dataset$group == "HER2-low","mRNA"])
  print(paste0(gene," ",round(p1$p.value,5)))
}
