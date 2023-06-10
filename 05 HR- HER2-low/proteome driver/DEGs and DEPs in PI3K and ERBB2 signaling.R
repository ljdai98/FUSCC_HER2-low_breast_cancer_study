c2_genesets <- read.gmt("c2.all.v7.5.symbols.gmt")
genesets <- c("REACTOME_PI3K_AKT_ACTIVATION","REACTOME_PI3K_AKT_SIGNALING_IN_CANCER","REACTOME_SIGNALING_BY_ERBB2")
gene <- unique(c2_genesets[c2_genesets$term%in%genesets,"gene"])

pheno <- FUSCC_HER2_low_project_Cohort.Info
pheno <- pheno[pheno$HR_status=="Negative" & pheno$HER2_low_status_RE=="HER2_low" &is.na(pheno$PAM50)==F,]

sample_nonbasal <- pheno[pheno$PAM50!="Basal","PatientCode"]
sample_basal <- pheno[pheno$PAM50=="Basal","PatientCode"]

exp <- FUSCC_HER2_low_project_RNA_seq_log2FPKM %>% as.matrix()
# exp <- FUSCC_HER2_low_project_TMT_PRO_normalized %>% as.matrix()
exp <- exp[,substring(colnames(exp),10,10)=="T"]
colnames(exp) <- substring(colnames(exp),1,4)

exp_nonbasal <- exp[rownames(exp)%in%gene,colnames(exp)%in%sample_nonbasal]
exp_basal <- exp[rownames(exp)%in%gene,colnames(exp)%in%sample_basal]

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

DEG_RNA <- as.data.frame(DEG)
colnames(DEG_RNA)[3] <- "log2FC_nonbasal/basal_RNA"
# DEG_pro <- as.data.frame(DEG)
# colnames(DEG_pro)[3] <- "log2FC_nonbasal/basal_pro"

DEG_RNA$symbol <- rownames(DEG_RNA)
DEG_pro$symbol <- rownames(DEG_pro)

data <- merge(DEG_RNA[,c(3,6)],DEG_pro[,c(3,6)],all.x = F,all.y = F,by = "symbol")

data[data$`log2FC_nonbasal/basal_RNA`>0 & data$`log2FC_nonbasal/basal_pro`>0,"significant"] <- "up"
data[data$`log2FC_nonbasal/basal_RNA`<0 & data$`log2FC_nonbasal/basal_pro`<0,"significant"] <- "down"
data[is.na(data$significant)==T,"significant"] <- "others"

ggplot(data=data, aes(x = `log2FC_nonbasal/basal_RNA`,y = `log2FC_nonbasal/basal_pro`, colour=significant)) +
  geom_point()+
  scale_color_manual(values=c("#4393C3","#d2dae2","#D6604D"))+
  xlim(c(-1.4,1.4))+ylim(c(-1.2,1.2))+
  geom_hline(yintercept=0,linetype=4)+              
  geom_vline(xintercept=0,linetype=4)+
  theme_bw()+geom_text_repel(data=data[data$significant!="others",],
                             aes(label=symbol),col="black",alpha = 0.8,size=3)+scale_size(trans="exp")


