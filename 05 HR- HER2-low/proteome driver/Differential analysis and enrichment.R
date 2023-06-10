pheno <- FUSCC_HER2_low_project_Cohort.Info
pheno <- pheno[pheno$HR_status=="Negative" & pheno$HER2_low_status_RE=="HER2_low" &is.na(pheno$PAM50)==F,]

sample_nonbasal <- pheno[pheno$PAM50!="Basal","PatientCode"]
sample_basal <- pheno[pheno$PAM50=="Basal","PatientCode"]

exp <- FUSCC_HER2_low_project_RNA_seq_log2FPKM %>% as.matrix()
# exp <- FUSCC_HER2_low_project_TMT_PRO_normalized %>% as.matrix()

exp <- exp[,substring(colnames(exp),10,10)=="T"]
colnames(exp) <- substring(colnames(exp),1,4)

exp_nonbasal <- exp[,colnames(exp)%in%sample_nonbasal]
exp_basal <- exp[,colnames(exp)%in%sample_basal]

DEG <- matrix(ncol = 5,nrow = nrow(exp_nonbasal))
rownames(DEG) <- rownames(exp_nonbasal)
colnames(DEG) <- c("median_nonbasal","median_basal","log2FC_nonbasal/basal","pvalue","fdr")

for (i in rownames(DEG)) {
  DEG[i,1] <- median(as.numeric(na.omit(exp_nonbasal[i,])))
  DEG[i,2] <- median(as.numeric(na.omit(exp_basal[i,])))
  DEG[i,3] <- mean(as.numeric(na.omit(exp_nonbasal[i,])))-mean(as.numeric(na.omit(exp_basal[i,])))
  DEG[i,4] <- wilcox.test(as.numeric(na.omit(exp_nonbasal[i,])),as.numeric(na.omit(exp_basal[i,])))$p.value
  print(which(rownames(DEG)==i))
}
DEG[,5] <- p.adjust(DEG[,4],method = "fdr") 

data <- as.data.frame(DEG)

up_RNA <- rownames(data[data$`log2FC_nonbasal/basal` > 1.5 & data$fdr<0.05,])
down_RNA <- rownames(data[data$`log2FC_nonbasal/basal` < -1.5 & data$fdr<0.05,])

# up_pro <- rownames(data[data$`log2FC_nonbasal/basal` > 0.5 & data$fdr<0.2,])
# down_pro <- rownames(data[data$`log2FC_nonbasal/basal` < -0.5 & data$fdr<0.2,])

gene <- up_RNA

gene <- mapIds(org.Hs.eg.db,gene,column = 'ENTREZID', keytype = 'SYMBOL',multiVals = 'filter')

GO <- enrichGO(gene=na.omit(gene),OrgDb=org.Hs.eg.db,ont="ALL",pAdjustMethod = "fdr",
               pvalueCutoff=0.05,qvalueCutoff  = 0.3) %>% as.data.frame()

GO <- GO[order(-GO$Count),]
GO <- GO[1:10,]
GO$number <- factor(rev(1:nrow(GO)))

ggplot(data=GO, aes(x=number,y=Count))+
  geom_bar(stat="identity", width=0.8) + coord_flip()+ 
  theme_test()+ 
  scale_x_discrete(labels=rev(GO$Description)) +
  theme(axis.text=element_text(size = 10,face = "bold", color="black")) +            
  labs(title = "Enriched GO terms")

