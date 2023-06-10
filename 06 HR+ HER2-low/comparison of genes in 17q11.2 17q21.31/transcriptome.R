table(FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE)
# HER2_0      HER2_low HER2_positive 
# 91           434           182 


##
her2_0_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_0"),"PatientCode"]
her2_low_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_low"),"PatientCode"]
her2_pos_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_positive"),"PatientCode"]

group_list<-c(rep("HER2_0",length(her2_0_list)),rep("HER2_low",length(her2_low_list)),rep("HER2_positive",length(her2_pos_list)))
names(group_list)<-c(her2_0_list,her2_low_list,her2_pos_list)
table(group_list)

# 17q21.31
HR_pos_names <- FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_low","HER2_0") &
                                                     FUSCC_HER2_low_project_Cohort.Info$Clinical_Subtype %in% c("HR+HER2-","HR+HER2+"),"PatientCode"]

HR_pos_names <- HR_pos_names[HR_pos_names %in% colnames(FUSCC_HER2_low_project_GISTICpeaks.del.thre)]
group_list_temp <- as.character(FUSCC_HER2_low_project_GISTICpeaks.del.thre["Del_Peak.43_17q21.31",HR_pos_names])
names(group_list_temp) <- HR_pos_names
names(group_list_temp) <- paste0(names(group_list_temp),"_RNA_T")
group_list_temp <- group_list_temp[names(group_list_temp) %in% colnames(FUSCC_HER2_low_project_RNA_seq_count)]
group_list_temp <- gsub("0","Others",group_list_temp)
group_list_temp <- gsub("1|2","Loss/del",group_list_temp)
table(group_list_temp)

gene_list <- FUSCC_HER2_low_project_GISTICgene_Cytoband[FUSCC_HER2_low_project_GISTICgene_Cytoband$Cytoband %in% c("17q21.31"),"Gene"]
gene_list <- gene_list[gene_list %in% rownames(FUSCC_HER2_low_project_RNA_seq_log2FPKM)]

choose_gene <- gene_list

compare_matrix <- matrix(NA,nrow = length(choose_gene),ncol = 6)
rownames(compare_matrix) <- choose_gene
colnames(compare_matrix) <- c("Median_Others","Median_Loss/del","P","P.adj","log2FC","For_size")
compare_matrix <- as.data.frame(compare_matrix)

group1 <- names(group_list_temp)[group_list_temp %in% "Others"]
group2 <- names(group_list_temp)[group_list_temp %in% "Loss/del"]


for (i in choose_gene){
  compare_matrix[i,1] <- median(as.numeric(FUSCC_HER2_low_project_RNA_seq_log2FPKM[i,group1]),na.rm = T)
  compare_matrix[i,2] <- median(as.numeric(FUSCC_HER2_low_project_RNA_seq_log2FPKM[i,group2]),na.rm = T)
  compare_matrix[i,3] <- wilcox.test(as.numeric(FUSCC_HER2_low_project_RNA_seq_log2FPKM[i,group1]),
                                     as.numeric(FUSCC_HER2_low_project_RNA_seq_log2FPKM[i,group2]))[["p.value"]]
  compare_matrix[i,5] <- mean(as.numeric(FUSCC_HER2_low_project_RNA_seq_log2FPKM[i,group1]),na.rm = T) - mean(as.numeric(FUSCC_HER2_low_project_RNA_seq_log2FPKM[i,group2]),na.rm = T)
}
compare_matrix[,4] <- p.adjust(compare_matrix[,3],method = "fdr")

compare_matrix_RNA <- compare_matrix
compare_matrix_RNA$For_size <- abs(compare_matrix_RNA$`Median_Loss/del` - compare_matrix_RNA$Median_Others)

compare_matrix_RNA <- compare_matrix_RNA[order(compare_matrix_RNA$log2FC,decreasing = T),]

write.csv(file = "compare_matrix_RNA_17q21.31_Lossdel_Others.csv",compare_matrix_RNA,quote = F)



library(ggplot2)
library(ggrepel)
library(dplyr)
data<-compare_matrix_RNA 
logFCLableFilter=0
creteria_logFC=0
creteria_adj_p=1
data$significant <- as.factor(data$P.adj<creteria_adj_p & abs(data$log2FC) > creteria_logFC)
data$gene <- rownames(data)
p1 <- ggplot(data=data, aes(x = log2FC,y = -log10(P.adj), colour=significant, size=For_size)) +
  geom_point(alpha=0.6,col="#d2dae2")+
  geom_point(data=subset(data, log2FC > 0),alpha=0.6,col="#D6604D")+
  geom_point(data=subset(data, log2FC < 0),alpha=0.6,col="#4393C3")+
  labs(x="log2(FoldChange)",
       y="-log10(P.adj)")+ggtitle("17q21.31")+
  theme_bw()+geom_text_repel(data=subset(data, abs(log2FC) > creteria_logFC & -log10(P.adj)>-log10(creteria_adj_p)),
                             aes(label=gene),col="black",alpha = 0.8,size=1.5,max.overlaps = 15)+scale_size(limits = c(0, 0.9),breaks=c(0,0.3,0.6,0.9))



##17q11.2
HR_pos_names <- FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_low","HER2_0") &
                                                     FUSCC_HER2_low_project_Cohort.Info$Clinical_Subtype %in% c("HR+HER2-","HR+HER2+"),"PatientCode"]

HR_pos_names <- HR_pos_names[HR_pos_names %in% colnames(FUSCC_HER2_low_project_GISTICpeaks.del.thre)]
group_list_temp <- as.character(FUSCC_HER2_low_project_GISTICpeaks.del.thre["Del_Peak.42_17q11.2",HR_pos_names])
names(group_list_temp) <- HR_pos_names
names(group_list_temp) <- paste0(names(group_list_temp),"_RNA_T")
group_list_temp <- group_list_temp[names(group_list_temp) %in% colnames(FUSCC_HER2_low_project_RNA_seq_count)]
group_list_temp <- gsub("0","Others",group_list_temp)
group_list_temp <- gsub("1|2","Loss/del",group_list_temp)
table(group_list_temp)

gene_list <- FUSCC_HER2_low_project_GISTICgene_Cytoband[FUSCC_HER2_low_project_GISTICgene_Cytoband$Cytoband %in% c("17q11.2"),"Gene"]
gene_list <- gene_list[gene_list %in% rownames(FUSCC_HER2_low_project_RNA_seq_log2FPKM)]

choose_gene <- gene_list

compare_matrix <- matrix(NA,nrow = length(choose_gene),ncol = 6)
rownames(compare_matrix) <- choose_gene
colnames(compare_matrix) <- c("Median_Others","Median_Loss/del","P","P.adj","log2FC","For_size")
compare_matrix <- as.data.frame(compare_matrix)

group1 <- names(group_list_temp)[group_list_temp %in% "Others"]
group2 <- names(group_list_temp)[group_list_temp %in% "Loss/del"]


for (i in choose_gene){
  compare_matrix[i,1] <- median(as.numeric(FUSCC_HER2_low_project_RNA_seq_log2FPKM[i,group1]),na.rm = T)
  compare_matrix[i,2] <- median(as.numeric(FUSCC_HER2_low_project_RNA_seq_log2FPKM[i,group2]),na.rm = T)
  compare_matrix[i,3] <- wilcox.test(as.numeric(FUSCC_HER2_low_project_RNA_seq_log2FPKM[i,group1]),
                                     as.numeric(FUSCC_HER2_low_project_RNA_seq_log2FPKM[i,group2]))[["p.value"]]
  compare_matrix[i,5] <- mean(as.numeric(FUSCC_HER2_low_project_RNA_seq_log2FPKM[i,group1]),na.rm = T) - mean(as.numeric(FUSCC_HER2_low_project_RNA_seq_log2FPKM[i,group2]),na.rm = T)
}
compare_matrix[,4] <- p.adjust(compare_matrix[,3],method = "fdr")

compare_matrix_RNA <- compare_matrix
compare_matrix_RNA$For_size <- abs(compare_matrix_RNA$`Median_Loss/del` - compare_matrix_RNA$Median_Others)

compare_matrix_RNA <- compare_matrix_RNA[order(compare_matrix_RNA$log2FC,decreasing = T),]

write.csv(file = "compare_matrix_RNA_17q11.2_Lossdel_Others.csv",compare_matrix_RNA,quote = F)


library(ggplot2)
library(ggrepel)
library(dplyr)
data<-compare_matrix_RNA 
logFCLableFilter=0
creteria_logFC=0
creteria_adj_p=1
data$significant <- as.factor(data$P.adj<creteria_adj_p & abs(data$log2FC) > creteria_logFC)
data$gene <- rownames(data)
p2 <- ggplot(data=data, aes(x = log2FC,y = -log10(P.adj), colour=significant, size=For_size)) +
  geom_point(alpha=0.6,col="#d2dae2")+
  geom_point(data=subset(data, log2FC > 0),alpha=0.6,col="#D6604D")+
  geom_point(data=subset(data, log2FC < 0),alpha=0.6,col="#4393C3")+
  labs(x="log2(FoldChange)",
       y="-log10(P.adj)")+ggtitle("17q11.2")+
  theme_bw()+geom_text_repel(data=subset(data, abs(log2FC) > creteria_logFC & -log10(P.adj)>-log10(creteria_adj_p)),
                             aes(label=gene),col="black",alpha = 0.8,size=1.5,max.overlaps = 15)+scale_size(limits = c(0, 0.9),breaks=c(0,0.3,0.6,0.9))

library(patchwork)
pdf('Vol_plot_RNA_combo_Lossdel_Others.pdf',width = 8,height = 4) 
p1+p2
dev.off()
