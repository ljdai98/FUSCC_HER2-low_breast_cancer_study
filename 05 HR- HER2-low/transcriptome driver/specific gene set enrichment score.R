table(FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE)
# HER2_0      HER2_low HER2_positive 
# 91           434           182 

HR_neg_HER2_0_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HR_status %in% c("Negative")& FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_0"),"PatientCode"]
HR_neg_HER2_low_list_basal<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HR_status %in% c("Negative") & FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_low") & FUSCC_HER2_low_project_Cohort.Info$PAM50 %in% c("Basal"),"PatientCode"]
HR_neg_HER2_low_list_nonbasal<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HR_status %in% c("Negative")& FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_low") & FUSCC_HER2_low_project_Cohort.Info$PAM50 %in% c("Her2","LumA","LumB","Normal"),"PatientCode"]
HR_neg_HER2_low_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HR_status %in% c("Negative")& FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_low"),"PatientCode"]
HR_neg_HER2_pos_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HR_status %in% c("Negative")& FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_positive"),"PatientCode"]


group_list<-c(rep("HER2_0",length(HR_neg_HER2_0_list)),
              rep("HER2_low_basal",length(HR_neg_HER2_low_list_basal)),
              rep("HER2_low_nonbasal",length(HR_neg_HER2_low_list_nonbasal)),
              rep("HER2_pos",length(HR_neg_HER2_pos_list)))

names(group_list)<-c(HR_neg_HER2_0_list,HR_neg_HER2_low_list_basal,HR_neg_HER2_low_list_nonbasal,HR_neg_HER2_pos_list)
group_list<-factor(group_list,levels=c("HER2_0","HER2_low_basal","HER2_low_nonbasal","HER2_pos"))
table(group_list)

load("../../04 molecular landscape/GSVA/RNA_C2.Rdata")

#REACTOME_PI3K_AKT_ACTIVATION
gene<-c("REACTOME_PI3K_AKT_ACTIVATION")
exprSet<-esmicro
plot_list<-group_list
names(plot_list)<-paste(names(plot_list),"_RNA_T",sep = "")
plot_list<-factor(plot_list)
plot_list<-plot_list[names(plot_list) %in% colnames(exprSet)]
my_comparisons <- list( c("HER2_0", "HER2_low_nonbasal"), c("HER2_low_basal", "HER2_low_nonbasal"), c("HER2_low_nonbasal", "HER2_pos"))
col <- c("#E41A1C","#4DAF4A","#984EA3","#377EB8")
names(col) <- c("HER2_0","HER2_low_basal", "HER2_low_nonbasal","HER2_pos")
dataset=data.frame(mRNA=as.matrix(exprSet)[gene,names(plot_list)],group=plot_list)
library(ggplot2)
library(ggpubr)
pdf(paste("boxplot_GSVA_",gene,"_HR_neg_4_SUBTYPE.pdf",sep=""),width = 3,height = 5)
ggboxplot(dataset, x = "group", y = "mRNA",ylab = "GSVA score",
          color = "group", palette = col,
          add = "jitter", title = gene,xlab  = FALSE)+stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")
dev.off()

#REACTOME_PI3K_EVENTS_IN_ERBB2_SIGNALING
gene<-c("REACTOME_PI3K_EVENTS_IN_ERBB2_SIGNALING")
exprSet<-esmicro
plot_list<-group_list
names(plot_list)<-paste(names(plot_list),"_RNA_T",sep = "")
plot_list<-factor(plot_list)
plot_list<-plot_list[names(plot_list) %in% colnames(exprSet)]
my_comparisons <- list( c("HER2_0", "HER2_low_nonbasal"), c("HER2_low_basal", "HER2_low_nonbasal"), c("HER2_low_nonbasal", "HER2_pos"))
col <- c("#E41A1C","#4DAF4A","#984EA3","#377EB8")
names(col) <- c("HER2_0","HER2_low_basal", "HER2_low_nonbasal","HER2_pos")
dataset=data.frame(mRNA=as.matrix(exprSet)[gene,names(plot_list)],group=plot_list)
library(ggplot2)
library(ggpubr)
pdf(paste("boxplot_GSVA_",gene,"_HR_neg_4_SUBTYPE.pdf",sep=""),width = 3,height = 5)
ggboxplot(dataset, x = "group", y = "mRNA",ylab = "GSVA score",
          color = "group", palette = col,
          add = "jitter", title = gene,xlab  = FALSE)+stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")
dev.off()


#REACTOME_ERBB2_ACTIVATES_PTK6_SIGNALING
gene<-c("REACTOME_ERBB2_ACTIVATES_PTK6_SIGNALING")
exprSet<-esmicro
plot_list<-group_list
names(plot_list)<-paste(names(plot_list),"_RNA_T",sep = "")
plot_list<-factor(plot_list)
plot_list<-plot_list[names(plot_list) %in% colnames(exprSet)]
my_comparisons <- list( c("HER2_0", "HER2_low_nonbasal"), c("HER2_low_basal", "HER2_low_nonbasal"), c("HER2_low_nonbasal", "HER2_pos"))
col <- c("#E41A1C","#4DAF4A","#984EA3","#377EB8")
names(col) <- c("HER2_0","HER2_low_basal", "HER2_low_nonbasal","HER2_pos")
dataset=data.frame(mRNA=as.matrix(exprSet)[gene,names(plot_list)],group=plot_list)
library(ggplot2)
library(ggpubr)
pdf(paste("boxplot_GSVA_",gene,"_HR_neg_4_SUBTYPE.pdf",sep=""),width = 3,height = 5)
ggboxplot(dataset, x = "group", y = "mRNA",ylab = "GSVA score",
          color = "group", palette = col,
          add = "jitter", title = gene,xlab  = FALSE)+stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")
dev.off()


