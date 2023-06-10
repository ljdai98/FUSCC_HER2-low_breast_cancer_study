load("KEGG_lipid_metabolism_annotation.RData")
load("ssGSEA_GSVA_curated_pathway.RData")

pheno <- FUSCC_HER2_low_project_Cohort.Info

sample <- pheno[pheno$HR_status=="Negative" & pheno$HER2_low_status_RE=="HER2_low","PatientCode"]

data <- FUSCC_HER2_low_project_WES_Somatic[FUSCC_HER2_low_project_WES_Somatic$Hugo_Symbol=="PIK3CA",]
sample_ex <- sample[sample%in%data$Tumor_Sample_Barcode]
sample_ct <- setdiff(sample,sample_ex)

exp <- FUSCC_HER2_low_project_RNA_seq_log2FPKM
exp <- exp[,substring(colnames(exp),10,10)=="T"] %>% as.matrix()
colnames(exp) <- substring(colnames(exp),1,4)

exp <- exp[,colnames(exp)%in%c(sample_ex,sample_ct)] %>% t() %>% scale() %>% t()

exp_ex <- exp[rownames(exp)%in%lipid_genes,colnames(exp)%in%sample_ex]
exp_ct <- exp[rownames(exp)%in%lipid_genes,colnames(exp)%in%sample_ct]

result <- matrix(ncol = 2,nrow = nrow(exp_ex))
rownames(result) <- rownames(exp_ex)
colnames(result) <- c("median_wt","median_mut")

result[,1] <- apply(exp_ct,1,median,na.rm=T)
result[,2] <- apply(exp_ex,1,median,na.rm=T)

result <- result[is.na(result[,1])==F & is.na(result[,2])==F,]

anno_row <- as.data.frame(names(lipid_genes))
rownames(anno_row) <- lipid_genes

pheatmap(result,cluster_rows = F,cluster_cols = F,show_rownames = F,cellwidth = 10,cellheight = 1,annotation_row = anno_row)

pathway <- toupper(unique(names(lipid_genes)))
pathway <- gsub(" ","_",pathway)
pathway <- gsub("-","_",pathway)
pathway <- paste("KEGG_",pathway,sep = "")
pathway <- c(pathway,"KEGG_FATTY_ACID_METABOLISM")

score <- GSVA_RNA_low[rownames(GSVA_RNA_low)%in%pathway,colnames(GSVA_RNA_low)%in%c(sample_ct,sample_ex)]
score <- t(score) %>% scale() %>% t()

sample_ct <- sample_ct[sample_ct%in%colnames(score)]
sample_ex <- sample_ex[sample_ex%in%colnames(score)]

result <- matrix(ncol = 2,nrow = nrow(score))
rownames(result) <- rownames(score)
colnames(result) <- c("PIK3CA WT","PIK3CA Mut")

result[,1] <- apply(score[,sample_ct],1,median,na.rm=T)
result[,2] <- apply(score[,sample_ex],1,median,na.rm=T)

pheatmap(result,cluster_rows = F,cluster_cols = F)

pheno <- FUSCC_HER2_low_project_Cohort.Info

sample <- pheno[pheno$HR_status=="Negative" & pheno$HER2_low_status_RE=="HER2_low","PatientCode"]

exp <- FUSCC_HER2_low_project_RNA_seq_log2FPKM
exp <- exp[,substring(colnames(exp),10,10)=="T"] %>% as.matrix()
colnames(exp) <- substring(colnames(exp),1,4)
exp <- exp[,colnames(exp)%in%sample] %>% as.matrix()

cor <- matrix(ncol = 6,nrow = length(lipid_genes[lipid_genes%in%rownames(exp)]))
colnames(cor) <- c("rho_FGFR4","p_FGFR4","rho_PTK6","p_PTK6","rho_ERBB4","p_ERBB4")
rownames(cor) <- lipid_genes[lipid_genes%in%rownames(exp)]

for (i in rownames(cor)) {
  cor[i,1] <- cor.test(exp["FGFR4",],exp[i,],method = "spearman",exact = F)$estimate
  cor[i,2] <- cor.test(exp["FGFR4",],exp[i,],method = "spearman",exact = F)$p.value
  cor[i,3] <- cor.test(exp["PTK6",],exp[i,],method = "spearman",exact = F)$estimate
  cor[i,4] <- cor.test(exp["PTK6",],exp[i,],method = "spearman",exact = F)$p.value
  cor[i,5] <- cor.test(exp["ERBB4",],exp[i,],method = "spearman",exact = F)$estimate
  cor[i,6] <- cor.test(exp["ERBB4",],exp[i,],method = "spearman",exact = F)$p.value
}

cor <- as.data.frame(cor)
cor$fdr_FGFR4 <- p.adjust(cor$p_FGFR4,method = "fdr")
cor$fdr_PTK6 <- p.adjust(cor$p_PTK6,method = "fdr")
cor$fdr_ERBB4 <- p.adjust(cor$p_ERBB4,method = "fdr")

anno_row <- as.data.frame(names(lipid_genes))
rownames(anno_row) <- lipid_genes

pheatmap(cor[,c(1,3,5)],cluster_rows = F,cluster_cols = F,show_rownames = F,annotation_row = anno_row
         ,cellwidth = 10,cellheight = 1)

data <- exp[c("FGFR4","PTK6","ERBB4"),] %>% as.matrix()

pathway <- toupper(unique(names(lipid_genes)))
pathway <- gsub(" ","_",pathway)
pathway <- gsub("-","_",pathway)
pathway <- paste("KEGG_",pathway,sep = "")
pathway <- c(pathway,"KEGG_FATTY_ACID_METABOLISM")

score <- GSVA_RNA_all[rownames(GSVA_RNA_all)%in%pathway,colnames(GSVA_RNA_all)%in%sample] 

cor <- matrix(ncol = 6,nrow = nrow(score))
colnames(cor) <- c("rho_FGFR4","p_FGFR4","rho_PTK6","p_PTK6","rho_ERBB4","p_ERBB4")
rownames(cor) <- rownames(score)

for (i in rownames(cor)) {
  cor[i,1] <- cor.test(data["FGFR4",],score[i,],method = "spearman",exact = F)$estimate
  cor[i,2] <- cor.test(data["FGFR4",],score[i,],method = "spearman",exact = F)$p.value
  cor[i,3] <- cor.test(data["PTK6",],score[i,],method = "spearman",exact = F)$estimate
  cor[i,4] <- cor.test(data["PTK6",],score[i,],method = "spearman",exact = F)$p.value
  cor[i,5] <- cor.test(data["ERBB4",],score[i,],method = "spearman",exact = F)$estimate
  cor[i,6] <- cor.test(data["ERBB4",],score[i,],method = "spearman",exact = F)$p.value
}

cor <- as.data.frame(cor)
cor$fdr_FGFR4 <- p.adjust(cor$p_FGFR4,method = "fdr")
cor$fdr_PTK6 <- p.adjust(cor$p_PTK6,method = "fdr")
cor$fdr_ERBB4 <- p.adjust(cor$p_ERBB4,method = "fdr")

pheatmap(cor[,c(1,3,5)],cluster_rows = F,cluster_cols = F)


