table(FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE)
# HER2_0      HER2_low HER2_positive 
# 91           434           182 

HR_neg_HER2_0_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HR_status %in% c("Negative")& FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_0"),"PatientCode"]
HR_neg_HER2_low_list_basal<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HR_status %in% c("Negative") & FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_low") & FUSCC_HER2_low_project_Cohort.Info$PAM50 %in% c("Basal"),"PatientCode"]
HR_neg_HER2_low_list_nonbasal<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HR_status %in% c("Negative")& FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_low") & FUSCC_HER2_low_project_Cohort.Info$PAM50 %in% c("Her2","LumA","LumB","Normal"),"PatientCode"]
HR_neg_HER2_low_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HR_status %in% c("Negative")& FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_low"),"PatientCode"]
HR_neg_HER2_pos_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HR_status %in% c("Negative")& FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_positive"),"PatientCode"]


#gene correlation
library(ggplot2)
library(ggpubr)
dat=data.frame(ERBB4=as.numeric(FUSCC_HER2_low_project_RNA_seq_log2FPKM["ERBB4",paste(HR_neg_HER2_low_list_nonbasal,"_RNA_T",sep = "")]),
               PTK6=as.numeric(FUSCC_HER2_low_project_RNA_seq_log2FPKM["PTK6",paste(HR_neg_HER2_low_list_nonbasal,"_RNA_T",sep = "")]))
rownames(dat) <- paste(HR_neg_HER2_low_list_nonbasal,"_RNA_T",sep = "")

pdf(file="correlation_dotplot_of_ERBB4_PTK6_in_HR_neg_HER2_low_list_nonbasal.pdf",width = 3,height = 3)
ggplot(data=dat, aes(x=ERBB4, y=PTK6),main="LOH")+geom_point(color="red")+stat_smooth(method="lm",se=FALSE)+stat_cor(data=dat, method = "spearman",cor.coef.name="rho") #有序等级资料用spearman，连续资料用pearson
dev.off()

#gene correlation
library(ggplot2)
library(ggpubr)
dat=data.frame(ERBB4=as.numeric(FUSCC_HER2_low_project_RNA_seq_log2FPKM["ERBB4",paste(HR_neg_HER2_low_list_nonbasal,"_RNA_T",sep = "")]),
               FGFR4=as.numeric(FUSCC_HER2_low_project_RNA_seq_log2FPKM["FGFR4",paste(HR_neg_HER2_low_list_nonbasal,"_RNA_T",sep = "")]))
rownames(dat) <- paste(HR_neg_HER2_low_list_nonbasal,"_RNA_T",sep = "")

pdf(file="correlation_dotplot_of_ERBB4_FGFR4_in_HR_neg_HER2_low_list_nonbasal.pdf",width = 3,height = 3)
ggplot(data=dat, aes(x=ERBB4, y=FGFR4),main="LOH")+geom_point(color="red")+stat_smooth(method="lm",se=FALSE)+stat_cor(data=dat, method = "spearman",cor.coef.name="rho") #有序等级资料用spearman，连续资料用pearson
dev.off()

#gene expression
group_list<-c(rep("HER2_0",length(HR_neg_HER2_0_list)),
              rep("HER2_low_basal",length(HR_neg_HER2_low_list_basal)),
              rep("HER2_low_nonbasal",length(HR_neg_HER2_low_list_nonbasal)),
              rep("HER2_pos",length(HR_neg_HER2_pos_list)))

names(group_list)<-c(HR_neg_HER2_0_list,HR_neg_HER2_low_list_basal,HR_neg_HER2_low_list_nonbasal,HR_neg_HER2_pos_list)
group_list<-factor(group_list,levels=c("HER2_0","HER2_low_basal","HER2_low_nonbasal","HER2_pos"))
table(group_list)

#ERBB2
gene<-c("ERBB2")
exprSet<-FUSCC_HER2_low_project_RNA_seq_log2FPKM
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
pdf(paste("boxplot_mRNA_",gene,"_HR_neg_4_SUBTYPE.pdf",sep=""),width = 3,height = 5)
ggboxplot(dataset, x = "group", y = "mRNA",ylab = "log2(FPKM+1)",
          color = "group", palette = col,
          add = "jitter", title = gene,xlab  = FALSE)+stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")
dev.off()

#ERBB2
gene<-c("ERBB2")
exprSet<-FUSCC_HER2_low_project_TMT_PRO_normalized
plot_list<-group_list
names(plot_list)<-paste(names(plot_list),"_PRO_T",sep = "")
plot_list<-factor(plot_list)
plot_list<-plot_list[names(plot_list) %in% colnames(exprSet)]
my_comparisons <- list( c("HER2_0", "HER2_low_nonbasal"), c("HER2_low_basal", "HER2_low_nonbasal"), c("HER2_low_nonbasal", "HER2_pos"))
col <- c("#E41A1C","#4DAF4A","#984EA3","#377EB8")
names(col) <- c("HER2_0","HER2_low_basal", "HER2_low_nonbasal","HER2_pos")
dataset=data.frame(PRO=as.matrix(exprSet)[gene,names(plot_list)],group=plot_list)
library(ggplot2)
library(ggpubr)
pdf(paste("boxplot_PRO_",gene,"_HR_neg_4_SUBTYPE.pdf",sep=""),width = 3,height = 5)
ggboxplot(dataset, x = "group", y = "PRO",ylab = "Protein level",
          color = "group", palette = col,
          add = "jitter", title = gene,xlab  = FALSE)+stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")
dev.off()

#PTK6
gene<-c("PTK6")
exprSet<-FUSCC_HER2_low_project_RNA_seq_log2FPKM
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
pdf(paste("boxplot_mRNA_",gene,"_HR_neg_4_SUBTYPE.pdf",sep=""),width = 2,height = 5)
ggboxplot(dataset, x = "group", y = "mRNA",ylab = "log2(FPKM+1)",
          color = "group", palette = col,
          add = "jitter", title = gene,xlab  = FALSE)+stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")
dev.off()

#FGFR4
gene<-c("FGFR4")
exprSet<-FUSCC_HER2_low_project_RNA_seq_log2FPKM
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
pdf(paste("boxplot_mRNA_",gene,"_HR_neg_4_SUBTYPE.pdf",sep=""),width = 2,height = 5)
ggboxplot(dataset, x = "group", y = "mRNA",ylab = "log2(FPKM+1)",
          color = "group", palette = col,
          add = "jitter", title = gene,xlab  = FALSE)+stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")
dev.off()

#ERBB4
gene<-c("ERBB4")
exprSet<-FUSCC_HER2_low_project_RNA_seq_log2FPKM
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
pdf(paste("boxplot_mRNA_",gene,"_HR_neg_4_SUBTYPE.pdf",sep=""),width = 2,height = 5)
ggboxplot(dataset, x = "group", y = "mRNA",ylab = "log2(FPKM+1)",
          color = "group", palette = col,
          add = "jitter", title = gene,xlab  = FALSE)+stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")
dev.off()


