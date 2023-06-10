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
names(group_list_temp) <- paste0(names(group_list_temp),"_PRO_T")
group_list_temp <- group_list_temp[names(group_list_temp) %in% colnames(FUSCC_HER2_low_project_TMT_PRO_normalized)]
group_list_temp <- gsub("0","Others",group_list_temp)
group_list_temp <- gsub("1|2","Loss/del",group_list_temp)
table(group_list_temp)

gene_list <- FUSCC_HER2_low_project_GISTICgene_Cytoband[FUSCC_HER2_low_project_GISTICgene_Cytoband$Cytoband %in% c("17q21.31"),"Gene"]
gene_list <- gene_list[gene_list %in% rownames(FUSCC_HER2_low_project_TMT_PRO_normalized)]

choose_gene <- gene_list


compare_matrix <- matrix(NA,nrow = length(choose_gene),ncol = 6)
rownames(compare_matrix) <- choose_gene
colnames(compare_matrix) <- c("Median_Others","Median_Loss/del","P","P.adj","log2FC","For_size")
compare_matrix <- as.data.frame(compare_matrix)

group1 <- names(group_list_temp)[group_list_temp %in% "Others"]
group2 <- names(group_list_temp)[group_list_temp %in% "Loss/del"]


for (i in choose_gene){
  compare_matrix[i,1] <- median(as.numeric(FUSCC_HER2_low_project_TMT_PRO_normalized[i,group1]),na.rm = T)
  compare_matrix[i,2] <- median(as.numeric(FUSCC_HER2_low_project_TMT_PRO_normalized[i,group2]),na.rm = T)
  compare_matrix[i,3] <- wilcox.test(as.numeric(FUSCC_HER2_low_project_TMT_PRO_normalized[i,group1]),
                                     as.numeric(FUSCC_HER2_low_project_TMT_PRO_normalized[i,group2]))[["p.value"]]
  compare_matrix[i,5] <- mean(as.numeric(FUSCC_HER2_low_project_TMT_PRO_normalized[i,group1]),na.rm = T) - mean(as.numeric(FUSCC_HER2_low_project_TMT_PRO_normalized[i,group2]),na.rm = T)
}
compare_matrix[,4] <- p.adjust(compare_matrix[,3],method = "fdr")


compare_matrix_RNA <- compare_matrix
compare_matrix_RNA$For_size <- abs(compare_matrix_RNA$`Median_Loss/del` - compare_matrix_RNA$Median_Others)

compare_matrix_RNA <- compare_matrix_RNA[order(compare_matrix_RNA$log2FC,decreasing = T),]

write.csv(file = "compare_matrix_RNA_17q21.31_Lossdel_Others(PRO).csv",compare_matrix_RNA,quote = F)

##17q11.2
HR_pos_names <- FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_low","HER2_0") &
                                                     FUSCC_HER2_low_project_Cohort.Info$Clinical_Subtype %in% c("HR+HER2-","HR+HER2+"),"PatientCode"]

HR_pos_names <- HR_pos_names[HR_pos_names %in% colnames(FUSCC_HER2_low_project_GISTICpeaks.del.thre)]
group_list_temp <- as.character(FUSCC_HER2_low_project_GISTICpeaks.del.thre["Del_Peak.42_17q11.2",HR_pos_names])
names(group_list_temp) <- HR_pos_names
names(group_list_temp) <- paste0(names(group_list_temp),"_PRO_T")
group_list_temp <- group_list_temp[names(group_list_temp) %in% colnames(FUSCC_HER2_low_project_TMT_PRO_normalized)]
group_list_temp <- gsub("0","Others",group_list_temp)
group_list_temp <- gsub("1|2","Loss/del",group_list_temp)
table(group_list_temp)

gene_list <- FUSCC_HER2_low_project_GISTICgene_Cytoband[FUSCC_HER2_low_project_GISTICgene_Cytoband$Cytoband %in% c("17q11.2"),"Gene"]
gene_list <- gene_list[gene_list %in% rownames(FUSCC_HER2_low_project_TMT_PRO_normalized)]

choose_gene <- gene_list


compare_matrix <- matrix(NA,nrow = length(choose_gene),ncol = 6)
rownames(compare_matrix) <- choose_gene
colnames(compare_matrix) <- c("Median_Others","Median_Loss/del","P","P.adj","log2FC","For_size")
compare_matrix <- as.data.frame(compare_matrix)

group1 <- names(group_list_temp)[group_list_temp %in% "Others"]
group2 <- names(group_list_temp)[group_list_temp %in% "Loss/del"]


for (i in choose_gene){
  compare_matrix[i,1] <- median(as.numeric(FUSCC_HER2_low_project_TMT_PRO_normalized[i,group1]),na.rm = T)
  compare_matrix[i,2] <- median(as.numeric(FUSCC_HER2_low_project_TMT_PRO_normalized[i,group2]),na.rm = T)
  compare_matrix[i,3] <- wilcox.test(as.numeric(FUSCC_HER2_low_project_TMT_PRO_normalized[i,group1]),
                                     as.numeric(FUSCC_HER2_low_project_TMT_PRO_normalized[i,group2]))[["p.value"]]
  compare_matrix[i,5] <- mean(as.numeric(FUSCC_HER2_low_project_TMT_PRO_normalized[i,group1]),na.rm = T) - mean(as.numeric(FUSCC_HER2_low_project_TMT_PRO_normalized[i,group2]),na.rm = T)
}
compare_matrix[,4] <- p.adjust(compare_matrix[,3],method = "fdr")


compare_matrix_RNA <- compare_matrix
compare_matrix_RNA$For_size <- abs(compare_matrix_RNA$`Median_Loss/del` - compare_matrix_RNA$Median_Others)

compare_matrix_RNA <- compare_matrix_RNA[order(compare_matrix_RNA$log2FC,decreasing = T),]

write.csv(file = "compare_matrix_RNA_17q11.2_Lossdel_Others(PRO).csv",compare_matrix_RNA,quote = F)


##dot_plot  17q11.2

q11_2_pro <- read.csv("compare_matrix_RNA_17q11.2_Lossdel_Others(PRO).csv",row.names = 1)
q11_2_rna <- read.csv("compare_matrix_RNA_17q11.2_Lossdel_Others.csv",row.names = 1)
common_genes <- intersect(rownames(q11_2_pro),rownames(q11_2_rna))
data=data.frame(`log2FC_RNA`=q11_2_rna[common_genes,"log2FC"],
                `log2FC_PRO`=q11_2_pro[common_genes,"log2FC"])
data$significant <- "0"
data[data$`log2FC_RNA` > 0 & data$`log2FC_PRO` >0,"significant"] <- "1"
data[data$`log2FC_RNA` < 0 & data$`log2FC_PRO` <0,"significant"] <- "2"
data$symbol <- common_genes
data$significant

library(ggplot2)
library(ggrepel)
library(dplyr)

pdf("Dot_plot_log2FC_RNA_log2FC_PRO_17q11.2.pdf",width = 5,height = 4)
ggplot(data=data, aes(x = log2FC_RNA,y = log2FC_PRO, colour=significant)) +
  geom_point()+
  scale_color_manual(values=c("#d2dae2","#D6604D","#4393C3"))+
  geom_hline(yintercept=0,linetype=4)+              
  geom_vline(xintercept=0,linetype=4)+
  theme_bw()+geom_text_repel(data=data[data$significant!="others",],
                             aes(label=symbol),col="black",alpha = 0.8,size=3)+scale_size(trans="exp")
dev.off()
  
  
##dot_plot  17q21.31

q21_31_pro <- read.csv("compare_matrix_RNA_17q21.31_Lossdel_Others(PRO).csv",row.names = 1)
q21_31_rna <- read.csv("compare_matrix_RNA_17q21.31_Lossdel_Others.csv",row.names = 1)
common_genes <- intersect(rownames(q21_31_pro),rownames(q21_31_rna))
data=data.frame(`log2FC_RNA`=q21_31_rna[common_genes,"log2FC"],
                `log2FC_PRO`=q21_31_pro[common_genes,"log2FC"])
data$significant <- "0"
data[data$`log2FC_RNA` > 0 & data$`log2FC_PRO` >0,"significant"] <- "1"
data[data$`log2FC_RNA` < 0 & data$`log2FC_PRO` <0,"significant"] <- "2"
data$symbol <- common_genes
data$significant

library(ggplot2)
library(ggrepel)
library(dplyr)

pdf("Dot_plot_log2FC_RNA_log2FC_PRO_17q21.31.pdf",width = 5,height = 4)
ggplot(data=data, aes(x = log2FC_RNA,y = log2FC_PRO, colour=significant)) +
  geom_point()+
  scale_color_manual(values=c("#d2dae2","#D6604D","#4393C3"))+
  geom_hline(yintercept=0,linetype=4)+              
  geom_vline(xintercept=0,linetype=4)+
  theme_bw()+geom_text_repel(data=data[data$significant!="others",],
                             aes(label=symbol),col="black",alpha = 0.8,size=3)+scale_size(trans="exp")
dev.off()

