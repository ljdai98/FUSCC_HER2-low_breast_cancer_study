table(FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE)
# HER2_0      HER2_low HER2_positive 
# 91           434           182 


meta_gene <- read.csv("lipid_genes_Kegg.csv")[,2]


## DEG hr_POS her2_low vs her2_0
list_a_samplenames<-rownames(FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_0") & 
                                                                         FUSCC_HER2_low_project_Cohort.Info$HR_status %in% c("Positive"),])
list_a_typename<-c("HR_pos_HER2_0")

list_b_samplenames<-rownames(FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_low") & 
                                                                         FUSCC_HER2_low_project_Cohort.Info$HR_status %in% c("Positive"),])
list_b_typename<-c("HR_pos_HER2_low")


exprSet<-round(FUSCC_HER2_low_project_RNA_seq_count,0)
exprSet<-apply(exprSet,2,as.numeric)
rownames(exprSet)<-rownames(FUSCC_HER2_low_project_RNA_seq_count)

group_list<-c(rep(list_a_typename,length(list_a_samplenames)),
              rep(list_b_typename,length(list_b_samplenames)))
names(group_list)<-c(list_a_samplenames,list_b_samplenames)
names(group_list)<-paste(names(group_list),"_RNA_T",sep = "")
group_list<-group_list[names(group_list) %in% colnames(exprSet)]
group_list<-factor(group_list)


exprSet<-exprSet[,names(group_list)]
colData<-data.frame(sample=names(group_list),condition=group_list)


dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = round(exprSet),
  colData = colData,
  design = ~ condition)

dds$condition <- relevel(dds$condition, ref = list_a_typename)

dds2 = DESeq2::DESeq(dds)
res = DESeq2::results(dds2)
res.table = as.data.frame(res@listData)
rownames(res.table) = res@rownames

res.table_csv <- res.table[,c("log2FoldChange","padj")]
colnames(res.table_csv) <- c("log2(Fold change)","Adjusted p")

sigegens <- rownames(res.table[abs(res.table$log2FoldChange)>1 & res.table$padj <0.05 ,])
sigegens <- sigegens[!sigegens %in% c("NA")]


library(stringr)
sigegens <- sigegens[!str_detect(sigegens, "NA.")]
length(sigegens)
# 129


res.tableX <- res.table
res.tableX$meta_gene <- "Others"
res.tableX[meta_gene,"meta_gene"] <- "Lipid metabolism genes"
res.tableX <- res.tableX[order(res.tableX$log2FoldChange,decreasing = T),]
write.csv(file = "HR_pos_HER2_0_low_deg.csv",res.tableX,quote = F)

res.table$meta_gene <- "AOthers"
res.table[meta_gene,"meta_gene"] <- "Lipid metabolism genes"


pdf(paste("mRNA_Volcano_",list_b_typename,"_(red)_vs_",list_a_typename,"_(blue)(xlim).pdf",sep = ""),width = 7,height = 5) #根据要画的火山图更改
library(ggplot2)
library(ggrepel)
library(dplyr)
data<-res.table
logFCLableFilter=1 
creteria_logFC=2
creteria_adj_p=0.05 
data$significant <- as.factor(data$padj<creteria_adj_p & abs(data$log2FoldChange) > creteria_logFC)
data$gene <- rownames(data)
ggplot(data=data, aes(x = log2FoldChange,y = -log10(padj), colour=significant)) +
  geom_point(alpha=0.4, size=3.5,col="#d2dae2",aes(shape=factor(meta_gene)))+
  geom_point(data=subset(data, log2FoldChange > logFCLableFilter & padj < (creteria_adj_p)),alpha=0.4, size=3.5,col="#D6604D",aes(shape=factor(meta_gene)))+
  geom_point(data=subset(data, log2FoldChange < -logFCLableFilter & padj < (creteria_adj_p)),alpha=0.4, size=3.5,col="#4393C3",aes(shape=factor(meta_gene)))+
  geom_vline(xintercept=c(-logFCLableFilter,logFCLableFilter),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(creteria_adj_p),lty=4,col="black",lwd=0.8) +
  labs(x="log2(FoldChange)",
       y="-log10(padj)")+xlim(-6.5,11.5)+ylim(0,16.2)+
  theme_bw()+geom_text_repel(data=subset(data, abs(log2FoldChange) > creteria_logFC & -log10(padj)>-log10(creteria_adj_p)),
                             aes(label=gene),col="black",alpha = 0.8)
dev.off()



## DEG hr_neg her2_low vs her2_0
list_a_samplenames<-rownames(FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_0") & 
                                                                  FUSCC_HER2_low_project_Cohort.Info$HR_status %in% c("Negative"),])
list_a_typename<-c("HR_neg_HER2_0")

list_b_samplenames<-rownames(FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_low") & 
                                                                  FUSCC_HER2_low_project_Cohort.Info$HR_status %in% c("Negative"),])
list_b_typename<-c("HR_neg_HER2_low")


exprSet<-round(FUSCC_HER2_low_project_RNA_seq_count,0)
exprSet<-apply(exprSet,2,as.numeric)
rownames(exprSet)<-rownames(FUSCC_HER2_low_project_RNA_seq_count)

group_list<-c(rep(list_a_typename,length(list_a_samplenames)),
              rep(list_b_typename,length(list_b_samplenames)))
names(group_list)<-c(list_a_samplenames,list_b_samplenames)
names(group_list)<-paste(names(group_list),"_RNA_T",sep = "")
group_list<-group_list[names(group_list) %in% colnames(exprSet)]
group_list<-factor(group_list)


exprSet<-exprSet[,names(group_list)]
colData<-data.frame(sample=names(group_list),condition=group_list)


dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = round(exprSet),
  colData = colData,
  design = ~ condition)

dds$condition <- relevel(dds$condition, ref = list_a_typename)

dds2 = DESeq2::DESeq(dds)
res = DESeq2::results(dds2)
res.table = as.data.frame(res@listData)
rownames(res.table) = res@rownames

res.table_csv <- res.table[,c("log2FoldChange","padj")]
colnames(res.table_csv) <- c("log2(Fold change)","Adjusted p")

sigegens <- rownames(res.table[abs(res.table$log2FoldChange)>1 & res.table$padj <0.05 ,])
sigegens <- sigegens[!sigegens %in% c("NA")]


library(stringr)
sigegens <- sigegens[!str_detect(sigegens, "NA.")]
length(sigegens)
# 281


res.tableX <- res.table
res.tableX$meta_gene <- "Others"
res.tableX[meta_gene,"meta_gene"] <- "Lipid metabolism genes"
res.tableX <- res.tableX[order(res.tableX$log2FoldChange,decreasing = T),]
write.csv(file = "HR_neg_HER2_0_low_deg.csv",res.tableX,quote = F)

res.table$meta_gene <- "AOthers"
res.table[meta_gene,"meta_gene"] <- "Lipid metabolism genes"


pdf(paste("mRNA_Volcano_",list_b_typename,"_(red)_vs_",list_a_typename,"_(blue)(xlim).pdf",sep = ""),width = 7,height = 5) #根据要画的火山图更改
library(ggplot2)
library(ggrepel)
library(dplyr)
data<-res.table
logFCLableFilter=1 
creteria_logFC=2
creteria_adj_p=0.05 
data$significant <- as.factor(data$padj<creteria_adj_p & abs(data$log2FoldChange) > creteria_logFC)
data$gene <- rownames(data)
ggplot(data=data, aes(x = log2FoldChange,y = -log10(padj), colour=significant)) +
  geom_point(alpha=0.4, size=3.5,col="#d2dae2",aes(shape=factor(meta_gene)))+
  geom_point(data=subset(data, log2FoldChange > logFCLableFilter & padj < (creteria_adj_p)),alpha=0.4, size=3.5,col="#D6604D",aes(shape=factor(meta_gene)))+
  geom_point(data=subset(data, log2FoldChange < -logFCLableFilter & padj < (creteria_adj_p)),alpha=0.4, size=3.5,col="#4393C3",aes(shape=factor(meta_gene)))+
  geom_vline(xintercept=c(-logFCLableFilter,logFCLableFilter),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(creteria_adj_p),lty=4,col="black",lwd=0.8) +
  labs(x="log2(FoldChange)",
       y="-log10(padj)")+xlim(-6.5,11.5)+ylim(0,16.2)+
  theme_bw()+geom_text_repel(data=subset(data, abs(log2FoldChange) > creteria_logFC & -log10(padj)>-log10(creteria_adj_p)),
                             aes(label=gene),col="black",alpha = 0.8)
dev.off()







