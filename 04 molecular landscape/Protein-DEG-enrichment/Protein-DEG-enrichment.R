library(ggplot2)
library(ggrepel)
library(dplyr)
library(org.Hs.eg.db)
library(clusterProfiler)

pheno <- FUSCC_HER2_low_project_Cohort.Info

pheno <- pheno[pheno$HR_status=="Negative",]
# pheno <- pheno[pheno$HR_status=="Positive",]

sample_low <- pheno[pheno$HER2_low_status_RE=="HER2_low","PatientCode"]
sample_0 <- pheno[pheno$HER2_low_status_RE=="HER2_0","PatientCode"]

exp <- FUSCC_HER2_low_project_TMT_PRO_normalized
exp <- exp[,substring(colnames(exp),10,10)=="T"] %>% as.matrix()
colnames(exp) <- substring(colnames(exp),1,4)

exp_low <- exp[,colnames(exp)%in%sample_low]
exp_0 <- exp[,colnames(exp)%in%sample_0]

DEG <- matrix(ncol = 5,nrow = nrow(exp_low))
rownames(DEG) <- rownames(exp_low)
colnames(DEG) <- c("median_low","median_0","log2FC_low/0","pvalue","fdr")

for (i in rownames(DEG)) {
  DEG[i,1] <- median(as.numeric(na.omit(exp_low[i,])))
  DEG[i,2] <- median(as.numeric(na.omit(exp_0[i,])))
  DEG[i,3] <- mean(as.numeric(na.omit(exp_low[i,])))-mean(as.numeric(na.omit(exp_0[i,])))
  DEG[i,4] <- wilcox.test(as.numeric(na.omit(exp_low[i,])),as.numeric(na.omit(exp_0[i,])))$p.value
  print(which(rownames(DEG)==i))
}
DEG[,5] <- p.adjust(DEG[,4],method = "fdr") 

data <- as.data.frame(DEG)
data <- data[order(data$`log2FC_low/0`),]
data$rank <- seq(1,nrow(data))

data$significant <- as.factor(abs(data$`log2FC_low/0`) > 0.5)
data$gene <- rownames(data)

data$For_size <- abs(data$median_low-data$median_0)

ggplot(data=data, aes(x = rank,y = `log2FC_low/0`, colour=significant, size=For_size)) +
  geom_point(alpha=0.6,col="#d2dae2")+
  geom_point(data=subset(data, `log2FC_low/0` > 0.5),alpha=0.6,col="#D6604D")+
  geom_point(data=subset(data, `log2FC_low/0` < -0.5),alpha=0.6,col="#4393C3")+
  labs(y="log2 Fold Change(HER2-low/HER2-0)")+ylim(c(-2,2))+
  geom_hline(yintercept=0.5,linetype=4)+
  geom_hline(yintercept=-0.5,linetype=4)+
  theme_bw()+geom_text_repel(data=subset(data, abs(`log2FC_low/0`) > 0.5),
                             aes(label=gene),col="black",alpha = 0.8,size=3)+scale_size(trans="exp")

up <- rownames(DEG[DEG[,3] >= 0.5,])
down <- rownames(DEG[DEG[,3] <= -0.5,])

gene <- up
# gene <- down

gene <- mapIds(org.Hs.eg.db,gene,column = 'ENTREZID', keytype = 'SYMBOL',multiVals = 'filter')

GO <- enrichGO(gene=na.omit(gene),OrgDb=org.Hs.eg.db,ont="ALL",pAdjustMethod = "fdr",
               pvalueCutoff=0.05,qvalueCutoff  = 0.05) %>% as.data.frame()

GO <- GO[order(-GO$Count),][1:10,]
GO$number <- factor(rev(1:nrow(GO)))

ggplot(data=GO, aes(x=number,y=Count))+
  geom_bar(stat="identity", width=0.8) + coord_flip()+ 
  theme_test()+ 
  scale_x_discrete(labels=rev(GO$Description)) +
  theme(axis.text=element_text(size = 10,face = "bold", color="black")) +            
  labs(title = "Enriched GO terms")



