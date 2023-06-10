table(FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE)
# HER2_0      HER2_low HER2_positive 
# 91           434           182 

HR_neg_HER2_0_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HR_status %in% c("Negative")& FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_0"),"PatientCode"]
HR_neg_HER2_low_list_basal<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HR_status %in% c("Negative") & FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_low") & FUSCC_HER2_low_project_Cohort.Info$PAM50 %in% c("Basal"),"PatientCode"]
HR_neg_HER2_low_list_nonbasal<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HR_status %in% c("Negative")& FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_low") & FUSCC_HER2_low_project_Cohort.Info$PAM50 %in% c("Her2","LumA","LumB","Normal"),"PatientCode"]
HR_neg_HER2_low_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HR_status %in% c("Negative")& FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_low"),"PatientCode"]
HR_neg_HER2_pos_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HR_status %in% c("Negative")& FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_positive"),"PatientCode"]


##############################################################
group_list<-c(rep("HER2_low_basal",length(HR_neg_HER2_low_list_basal)),
              rep("HER2_low_nonbasal",length(HR_neg_HER2_low_list_nonbasal)))
names(group_list)<-c(HR_neg_HER2_low_list_basal,HR_neg_HER2_low_list_nonbasal)
group_list<-factor(group_list,levels=c("HER2_low_basal","HER2_low_nonbasal"))
table(group_list)
# HER2_low_basal HER2_low_nonbasal 
# 46                20

c2_genesets<-clusterProfiler::read.gmt("../../04 molecular landscape/GSVA/c2.all.v7.5.symbols.gmt")
genesets = split(c2_genesets$gene, c2_genesets$term)

######
genesets_pick1<-genesets["REACTOME_PI3K_AKT_ACTIVATION"]
genesets_pick1[["REACTOME_PI3K_AKT_ACTIVATION"]]

genesets_pick2<-genesets["REACTOME_PI3K_AKT_SIGNALING_IN_CANCER"]
genesets_pick2[["REACTOME_PI3K_AKT_SIGNALING_IN_CANCER"]]

genesets_pick3<-genesets["REACTOME_SIGNALING_BY_ERBB2"]
genesets_pick3[["REACTOME_SIGNALING_BY_ERBB2"]]

choose_gene<-unique(c(unlist(genesets_pick1),unlist(genesets_pick2),unlist(genesets_pick3)))

choose_gene <- choose_gene[choose_gene %in% rownames(FUSCC_HER2_low_project_RNA_seq_log2FPKM)]

##RNA
compare_matrix <- matrix(NA,nrow = length(choose_gene),ncol = 6)
rownames(compare_matrix) <- choose_gene
colnames(compare_matrix) <- c("Median_nonbasal","Median_basal","P","P.adj","log2FC","For_size")
compare_matrix <- as.data.frame(compare_matrix)

HR_neg_HER2_low_list_nonbasal_RNA <- HR_neg_HER2_low_list_nonbasal[paste(HR_neg_HER2_low_list_nonbasal,"_RNA_T",sep = "") %in% colnames(FUSCC_HER2_low_project_RNA_seq_log2FPKM)]
HR_neg_HER2_low_list_basal_RNA <- HR_neg_HER2_low_list_basal[paste(HR_neg_HER2_low_list_basal,"_RNA_T",sep = "") %in% colnames(FUSCC_HER2_low_project_RNA_seq_log2FPKM)]


for (i in choose_gene){
  compare_matrix[i,1] <- median(as.numeric(FUSCC_HER2_low_project_RNA_seq_log2FPKM[i,paste(HR_neg_HER2_low_list_nonbasal_RNA,"_RNA_T",sep = "")]),na.rm = T)
  compare_matrix[i,2] <- median(as.numeric(FUSCC_HER2_low_project_RNA_seq_log2FPKM[i,paste(HR_neg_HER2_low_list_basal_RNA,"_RNA_T",sep = "")]),na.rm = T)
  compare_matrix[i,3] <- wilcox.test(as.numeric(FUSCC_HER2_low_project_RNA_seq_log2FPKM[i,paste(HR_neg_HER2_low_list_nonbasal_RNA,"_RNA_T",sep = "")]),
                                     as.numeric(FUSCC_HER2_low_project_RNA_seq_log2FPKM[i,paste(HR_neg_HER2_low_list_basal_RNA,"_RNA_T",sep = "")]))[["p.value"]]
  compare_matrix[i,5] <- mean(as.numeric(FUSCC_HER2_low_project_RNA_seq_log2FPKM[i,paste(HR_neg_HER2_low_list_nonbasal_RNA,"_RNA_T",sep = "")]),na.rm = T) - mean(as.numeric(FUSCC_HER2_low_project_RNA_seq_log2FPKM[i,paste(HR_neg_HER2_low_list_basal_RNA,"_RNA_T",sep = "")]),na.rm = T)
}
compare_matrix[,4] <- p.adjust(compare_matrix[,3],method = "fdr")


compare_matrix_RNA <- compare_matrix
compare_matrix_RNA$For_size <- abs(compare_matrix_RNA$Median_nonbasal - compare_matrix_RNA$Median_basal)

write.csv(file = "compare_matrix_RNA_basal_nonbasal.csv",compare_matrix_RNA,quote = F)

##vol
pdf('vol_plot_nonbasal_basal_basal(blue).pdf',width = 7,height = 5) 
library(ggplot2)
library(ggrepel)
library(dplyr)
data<-compare_matrix_RNA 
logFCLableFilter=1 
creteria_logFC=0
creteria_adj_p=0.05 
data$significant <- as.factor(data$P.adj<creteria_adj_p & abs(data$log2FC) > creteria_logFC)
data$gene <- rownames(data)
ggplot(data=data, aes(x = log2FC,y = -log10(P.adj), colour=significant, size=For_size)) +
  geom_point(alpha=0.6,col="#d2dae2")+
  geom_point(data=subset(data, log2FC > 1 & P.adj < 0.05),alpha=0.6,col="#D6604D")+
  geom_point(data=subset(data, log2FC < -1 & P.adj < 0.05),alpha=0.6,col="#4393C3")+
  labs(x="log2(FoldChange)",
       y="-log10(P.adj)")+xlim(c(-1.6,1.6))+ylim(c(0,3))+
  theme_bw()+geom_text_repel(data=subset(data, abs(log2FC) > creteria_logFC & -log10(P.adj)>-log10(creteria_adj_p)),
                             aes(label=gene),col="black",alpha = 0.8,size=3)+scale_size(trans="exp")
dev.off()



