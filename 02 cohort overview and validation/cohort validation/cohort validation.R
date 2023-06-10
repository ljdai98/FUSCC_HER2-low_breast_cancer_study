table(FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE)
# HER2_0      HER2_low HER2_positive 
# 91           434           182 


## HR Not stratified ----------------
list_a_samplenames<-rownames(FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_0") & 
                                                                         FUSCC_HER2_low_project_Cohort.Info$HR_status %in% c("Positive","Negative"),])
list_a_typename<-c("HER2_0")

list_b_samplenames<-rownames(FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_low") & 
                                                                         FUSCC_HER2_low_project_Cohort.Info$HR_status %in% c("Positive","Negative"),])
list_b_typename<-c("HER2_low")

list_c_samplenames<-rownames(FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_positive") & 
                                                                         FUSCC_HER2_low_project_Cohort.Info$HR_status %in% c("Positive","Negative"),])
list_c_typename<-c("HER2_pos")

group_list <- c(rep(list_a_typename,length(list_a_samplenames)),
                rep(list_b_typename,length(list_b_samplenames)),
                rep(list_c_typename,length(list_c_samplenames)))

names(group_list) <- c(list_a_samplenames,list_b_samplenames,list_c_samplenames)
table(group_list)

##CNV
gene <- "ERBB2"
Subtype = t(FUSCC_HER2_low_project_GISTICgene.thre)[names(group_list[names(group_list) %in% colnames(FUSCC_HER2_low_project_GISTICgene.thre)]),] ; Subtype<-as.data.frame(Subtype);Subtype$cluster = as.character(group_list[names(group_list) %in% colnames(FUSCC_HER2_low_project_GISTICgene.thre)])
ERBB2_Table <- xtabs(~ ERBB2 + cluster, data = Subtype)

Color_ERBB2<-c("#C2C2C4","#EC7878","#E21A21","#88A3D8")
names(Color_ERBB2)<-c("0","1","2","-1")

pdf(file="CNA_ERBB2_0_low_pos_subgroups.pdf",width = 4,height = 6)
barplot(apply(ERBB2_Table, 2, function(x){x/sum(x)}), col = Color_ERBB2[row.names(ERBB2_Table)] ,border = F,space = 1.25)
dev.off()

fisher.test(cbind(ERBB2_Table[,2],ERBB2_Table[,1]))
# p-value = 2.264e-07
fisher.test(cbind(ERBB2_Table[,2],ERBB2_Table[,3]))
# p-value < 2.2e-16


##ERBB2  PRO
exprSet<-FUSCC_HER2_low_project_TMT_PRO_normalized
plot_list<-group_list
names(plot_list)<-paste(names(plot_list),"_PRO_T",sep = "")
plot_list<-plot_list[names(plot_list) %in% colnames(exprSet)]


my_comparisons <- list(c(list_a_typename,list_b_typename),c(list_a_typename,list_c_typename),c(list_b_typename,list_c_typename))

dataset=data.frame(mRNA=as.matrix(exprSet)[gene,names(plot_list)],group=plot_list)

library(ggplot2)
library(ggpubr)
pdf(file="PRO_ERBB2_0_low_pos_subgroups.pdf",width = 3,height = 5)
ggboxplot(dataset, x = "group", y = "mRNA",
          color = "group", palette = "aaas",
          add = "jitter",title = gene,ylab = "Protein abundance")+stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")
dev.off()


##ERBB2 mRNA
gene <- "ERBB2"
exprSet<-FUSCC_HER2_low_project_RNA_seq_log2FPKM
plot_list<-group_list
names(plot_list)<-paste(names(plot_list),"_RNA_T",sep = "")
plot_list<-plot_list[names(plot_list) %in% colnames(exprSet)]


my_comparisons <- list(c(list_a_typename,list_b_typename),c(list_a_typename,list_c_typename),c(list_b_typename,list_c_typename))

dataset=data.frame(mRNA=as.matrix(exprSet)[gene,names(plot_list)],group=plot_list)

library(ggplot2)
library(ggpubr)
pdf(file="mRNA_ERBB2_0_low_pos_subgroups.pdf",width = 3,height = 5)
ggboxplot(dataset, x = "group", y = "mRNA",
          color = "group", palette = "aaas",
          add = "jitter",title = gene,ylab = "log2(FPKM+1)")+stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")
dev.off()


## HR stratified ------------------------
list_a_samplenames<-rownames(FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_0") & 
                                                                  FUSCC_HER2_low_project_Cohort.Info$HR_status %in% c("Positive"),])
list_a_typename<-c("HR_pos_HER2_0")

list_b_samplenames<-rownames(FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_low") & 
                                                                  FUSCC_HER2_low_project_Cohort.Info$HR_status %in% c("Positive"),])
list_b_typename<-c("HR_pos_HER2_low")

list_c_samplenames<-rownames(FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_positive") & 
                                                                  FUSCC_HER2_low_project_Cohort.Info$HR_status %in% c("Positive"),])
list_c_typename<-c("HR_pos_HER2_pos")

list_d_samplenames<-rownames(FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_0") & 
                                                                  FUSCC_HER2_low_project_Cohort.Info$HR_status %in% c("Negative"),])
list_d_typename<-c("HR_neg_HER2_0")

list_e_samplenames<-rownames(FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_low") & 
                                                                  FUSCC_HER2_low_project_Cohort.Info$HR_status %in% c("Negative"),])
list_e_typename<-c("HR_neg_HER2_low")

list_f_samplenames<-rownames(FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_positive") & 
                                                                  FUSCC_HER2_low_project_Cohort.Info$HR_status %in% c("Negative"),])
list_f_typename<-c("HR_neg_HER2_pos")


group_list <- c(rep(list_a_typename,length(list_a_samplenames)),
                rep(list_b_typename,length(list_b_samplenames)),
                rep(list_c_typename,length(list_c_samplenames)),
                rep(list_d_typename,length(list_d_samplenames)),
                rep(list_e_typename,length(list_e_samplenames)),
                rep(list_f_typename,length(list_f_samplenames)))

names(group_list) <- c(list_a_samplenames,list_b_samplenames,list_c_samplenames,
                       list_d_samplenames,list_e_samplenames,list_f_samplenames)
table(group_list)


##ERBB2 CNV
gene <- "ERBB2"
Subtype = t(FUSCC_HER2_low_project_GISTICgene.thre)[names(group_list[names(group_list) %in% colnames(FUSCC_HER2_low_project_GISTICgene.thre)]),] ; Subtype<-as.data.frame(Subtype);Subtype$cluster = as.character(group_list[names(group_list) %in% colnames(FUSCC_HER2_low_project_GISTICgene.thre)])
ERBB2_Table <- xtabs(~ ERBB2 + cluster, data = Subtype)
ERBB2_Table <- ERBB2_Table[,c(4:6,1:3)]

Color_ERBB2<-c("#C2C2C4","#EC7878","#E21A21","#88A3D8")
names(Color_ERBB2)<-c("0","1","2","-1")
pdf(file="CNA_ERBB2_HR_stratification_0_low_pos_subgroups.pdf",width = 6,height = 6)
barplot(apply(ERBB2_Table, 2, function(x){x/sum(x)}),space = 1.5, col = Color_ERBB2[row.names(ERBB2_Table)] ,border = F,cex.names  = 0.5)
dev.off()

fisher.test(cbind(ERBB2_Table[,1],ERBB2_Table[,2]))
#p-value = 5.454e-07
fisher.test(cbind(ERBB2_Table[,2],ERBB2_Table[,3]))
#p-value < 2.2e-16
fisher.test(cbind(ERBB2_Table[,4],ERBB2_Table[,5]))
#p-value = 0.4162
fisher.test(cbind(ERBB2_Table[,5],ERBB2_Table[,6]))
#p-value < 2.2e-16
fisher.test(cbind(ERBB2_Table[,1],ERBB2_Table[,4]))
#p-value = 0.3452
fisher.test(cbind(ERBB2_Table[,2],ERBB2_Table[,5]))
#p-value = 0.002035
fisher.test(cbind(ERBB2_Table[,3],ERBB2_Table[,6]))
#p-value = 0.2364



##ERBB2 pro
exprSet<-FUSCC_HER2_low_project_TMT_PRO_normalized

plot_list<-group_list
names(plot_list)<-paste(names(plot_list),"_PRO_T",sep = "")
plot_list<-plot_list[names(plot_list) %in% colnames(exprSet)]


my_comparisons <- list(c(list_a_typename,list_b_typename),c(list_b_typename,list_c_typename),
                       c(list_d_typename,list_e_typename),c(list_e_typename,list_f_typename),
                       c(list_a_typename,list_d_typename),c(list_b_typename,list_e_typename),c(list_b_typename,list_f_typename))
dataset=data.frame(mRNA=as.matrix(exprSet)[gene,names(plot_list)],group=plot_list)

library(ggplot2)
library(ggpubr)
pdf(file="PRO_ERBB2_HR_stratification_0_low_pos_subgroups.pdf",width = 5,height = 6)
ggboxplot(dataset, x = "group", y = "mRNA",
          color = "group", palette = "aaas",
          add = "jitter",title = gene,ylab = "Protein abundance",order = c(list_a_typename,list_b_typename,list_c_typename,list_d_typename,list_e_typename,list_f_typename))+stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")
dev.off()


##ERBB2 mRNA
gene <- "ERBB2"

exprSet<-FUSCC_HER2_low_project_RNA_seq_log2FPKM

plot_list<-group_list
names(plot_list)<-paste(names(plot_list),"_RNA_T",sep = "")
plot_list<-plot_list[names(plot_list) %in% colnames(exprSet)]


my_comparisons <- list(c(list_a_typename,list_b_typename),c(list_b_typename,list_c_typename),
                       c(list_d_typename,list_e_typename),c(list_e_typename,list_f_typename),
                       c(list_a_typename,list_d_typename),c(list_b_typename,list_e_typename),c(list_b_typename,list_f_typename))

dataset=data.frame(mRNA=as.matrix(exprSet)[gene,names(plot_list)],group=plot_list)

library(ggplot2)
library(ggpubr)
pdf(file="mRNA_ERBB2_HR_stratification_0_low_pos_subgroups.pdf",width = 5,height = 6)
ggboxplot(dataset, x = "group", y = "mRNA",
          color = "group", palette = "aaas",
          add = "jitter",title = gene,ylab = "log2(FPKM+1)",order = c(list_a_typename,list_b_typename,list_c_typename,list_d_typename,list_e_typename,list_f_typename))+stat_compare_means(comparisons = my_comparisons)
dev.off()

