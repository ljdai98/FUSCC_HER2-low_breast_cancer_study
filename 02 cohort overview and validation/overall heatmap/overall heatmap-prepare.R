table(FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE)
# HER2_0      HER2_low HER2_positive 
# 91           434           182 

her2_0_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_0"),"PatientCode"]
her2_low_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_low"),"PatientCode"]
her2_pos_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_positive"),"PatientCode"]

group_list<-c(rep("HER2_0",length(her2_0_list)),rep("HER2_low",length(her2_low_list)),rep("HER2_positive",length(her2_pos_list)))
names(group_list)<-c(her2_0_list,her2_low_list,her2_pos_list)
table(group_list)
# group_list
# HER2_0      HER2_low HER2_positive 
# 91           434           182 

##general meta dataframe#####################################
select_var <- c("PatientCode","HR_status","HER2_low_status_RE","HER2_low_status_long_RE","PAM50","Menopause","DMFS_status","Stage_pN")
clinical_anno <- FUSCC_HER2_low_project_Cohort.Info[names(group_list),select_var]

clinical_anno[,c("DMFS_status")] <- as.character(clinical_anno[,c("DMFS_status")])
clinical_anno[clinical_anno$Stage_pN %in% "pN0",c("Stage_pN")] <- "Negative"
clinical_anno[clinical_anno$Stage_pN %in% c("pN1","pN2","pN3"),c("Stage_pN")] <- "Positive"
colnames(clinical_anno)[colnames(clinical_anno) == "Stage_pN"] <- "Lymph_nodes"


clinical_anno[rownames(clinical_anno)[paste(rownames(clinical_anno),"_RNA_T",sep="") %in% colnames(FUSCC_HER2_low_project_RNA_seq_log2FPKM)],"HER2_RNA"] <- 
  as.numeric(FUSCC_HER2_low_project_RNA_seq_log2FPKM["ERBB2",paste(rownames(clinical_anno),"_RNA_T",sep="")[paste(rownames(clinical_anno),"_RNA_T",sep="") %in% colnames(FUSCC_HER2_low_project_RNA_seq_log2FPKM)]])

clinical_anno[rownames(clinical_anno)[paste(rownames(clinical_anno),"_PRO_T",sep="") %in% colnames(FUSCC_HER2_low_project_TMT_PRO_normalized)],"HER2_PRO"] <- 
  as.numeric(FUSCC_HER2_low_project_TMT_PRO_normalized["ERBB2",paste(rownames(clinical_anno),"_PRO_T",sep="")[paste(rownames(clinical_anno),"_PRO_T",sep="") %in% colnames(FUSCC_HER2_low_project_TMT_PRO_normalized)]])


clinical_anno <- clinical_anno[clinical_anno$HER2_low_status_RE %in% "HER2_low",]



her2_low_HR_pos_IHC_1_ori <- clinical_anno[clinical_anno$HR_status %in% c("Positive") & clinical_anno$HER2_low_status_long_RE %in% c("HER2_1"),"PatientCode"]
her2_low_HR_pos_IHC_2_ori <- clinical_anno[clinical_anno$HR_status %in% c("Positive") & clinical_anno$HER2_low_status_long_RE %in% c("HER2_2"),"PatientCode"]
her2_low_HR_neg_IHC_1_ori <- clinical_anno[clinical_anno$HR_status %in% c("Negative") & clinical_anno$HER2_low_status_long_RE %in% c("HER2_1"),"PatientCode"]
her2_low_HR_neg_IHC_2_ori <- clinical_anno[clinical_anno$HR_status %in% c("Negative") & clinical_anno$HER2_low_status_long_RE %in% c("HER2_2"),"PatientCode"]

##mutation matrix#####################################

library(maftools)

vc_cols = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(vc_cols) = c(  'Frame_Shift_Del',  'Missense_Mutation',  'Nonsense_Mutation',  'Multi_Hit',
                     'Frame_Shift_Ins',  'In_Frame_Ins',  'Splice_Site',  'In_Frame_Del')



maf_for_maf=FUSCC_HER2_low_project_WES_Somatic[FUSCC_HER2_low_project_WES_Somatic$Tumor_Sample_Barcode %in% c(her2_low_list),]
genes=read.csv("CAG.csv")[,1]
genes<-genes[genes %in% unique(maf_for_maf$Hugo_Symbol)]
maf_for_maf<-maf_for_maf[maf_for_maf$Hugo_Symbol %in% genes,]

clinicalData_for_maf=clinical_anno[unique(maf_for_maf$Tumor_Sample_Barcode),c("PatientCode","HER2_low_status_long_RE","HR_status")]
colnames(clinicalData_for_maf)[1]<-c("Tumor_Sample_Barcode")


laml = read.maf(maf = maf_for_maf,
                clinicalData = clinicalData_for_maf ,
                verbose = FALSE)

##rename the onco_matrix.txt
#pdf(file=paste("Mutation_waterfall_all_HER2_low.pdf"),width = 8,height = 6)
oncoplot(colors = vc_cols,maf = laml, draw_titv = TRUE,writeMatrix = T,clinicalFeatures = c("HER2_low_status_long_RE","HR_status"),sortByAnnotation = TRUE,
         annotationOrder=c("HER2_1","HER2_2"),top = 40,removeNonMutated=F)#writeMatrix可自动产生一个对应的txt文件
#dev.off()

file.rename("onco_matrix.txt","Mutation_waterfall_all_HER2_low_40.txt")


mut_ori_40 <- read.csv("Mutation_waterfall_all_HER2_low_40.txt",sep = "\t")
#colnames(mut_ori) <- paste(colnames(mut_ori),"_PRO_T",sep = "")##为了匹配到空白模板上


#DEG(HR+)

maf_for_maf=FUSCC_HER2_low_project_WES_Somatic[FUSCC_HER2_low_project_WES_Somatic$Tumor_Sample_Barcode %in% her2_low_HR_pos_IHC_1_ori,]
genes=read.csv("CAG.csv")[,1]
genes<-genes[genes %in% unique(maf_for_maf$Hugo_Symbol)]
maf_for_maf<-maf_for_maf[maf_for_maf$Hugo_Symbol %in% genes,]

clinicalData_for_maf=clinical_anno[unique(maf_for_maf$Tumor_Sample_Barcode),c("PatientCode","HER2_low_status_long_RE","HR_status")]
colnames(clinicalData_for_maf)[1]<-c("Tumor_Sample_Barcode")


laml1 = read.maf(maf = maf_for_maf,
                 clinicalData = clinicalData_for_maf ,
                 verbose = FALSE)


maf_for_maf=FUSCC_HER2_low_project_WES_Somatic[FUSCC_HER2_low_project_WES_Somatic$Tumor_Sample_Barcode %in% her2_low_HR_pos_IHC_2_ori,]
genes=read.csv("CAG.csv")[,1]
genes<-genes[genes %in% unique(maf_for_maf$Hugo_Symbol)]
maf_for_maf<-maf_for_maf[maf_for_maf$Hugo_Symbol %in% genes,]

clinicalData_for_maf=clinical_anno[unique(maf_for_maf$Tumor_Sample_Barcode),c("PatientCode","HER2_low_status_long_RE","HR_status")]
colnames(clinicalData_for_maf)[1]<-c("Tumor_Sample_Barcode")


laml2 = read.maf(maf = maf_for_maf,
                 clinicalData = clinicalData_for_maf ,
                 verbose = FALSE)


pt.vs.rt <- mafCompare(m1 = laml1, m2 = laml2, m1Name = 'HR_pos_IHC_1_ori', m2Name = 'HR_pos_IHC_2_ori', minMut = 5)
print(pt.vs.rt)
# $results
# Hugo_Symbol HR_pos_IHC_1_ori HR_pos_IHC_2_ori       pval        or      ci.up     ci.low   adjPval
# 1:        AKT1               18                6 0.01896619 3.0535288   9.675793 1.11991759 0.3314790
# 2:         NF1                6                0 0.03013445       Inf        Inf 1.12971295 0.3314790
# $SampleSummary
# Cohort SampleSize
# 1: HR_pos_IHC_1_ori        159
# 2: HR_pos_IHC_2_ori        150

sig_dif_names_hr_pos <- data.frame(pt.vs.rt$results)
sig_dif_names_hr_pos <- as.vector(sig_dif_names_hr_pos[sig_dif_names_hr_pos$pval < 0.05,"Hugo_Symbol"])
sig_dif_names_hr_pos
# [1] "AKT1" "NF1"

#DEG(HR-)
maf_for_maf=FUSCC_HER2_low_project_WES_Somatic[FUSCC_HER2_low_project_WES_Somatic$Tumor_Sample_Barcode %in% her2_low_HR_neg_IHC_1_ori,]
genes=read.csv("CAG.csv")[,1]
genes<-genes[genes %in% unique(maf_for_maf$Hugo_Symbol)]
maf_for_maf<-maf_for_maf[maf_for_maf$Hugo_Symbol %in% genes,]

clinicalData_for_maf=clinical_anno[unique(maf_for_maf$Tumor_Sample_Barcode),c("PatientCode","HER2_low_status_long_RE","HR_status")]
colnames(clinicalData_for_maf)[1]<-c("Tumor_Sample_Barcode")


laml1 = read.maf(maf = maf_for_maf,
                 clinicalData = clinicalData_for_maf ,
                 verbose = FALSE)


maf_for_maf=FUSCC_HER2_low_project_WES_Somatic[FUSCC_HER2_low_project_WES_Somatic$Tumor_Sample_Barcode %in% her2_low_HR_neg_IHC_2_ori,]
genes=read.csv("CAG.csv")[,1]
genes<-genes[genes %in% unique(maf_for_maf$Hugo_Symbol)]
maf_for_maf<-maf_for_maf[maf_for_maf$Hugo_Symbol %in% genes,]

clinicalData_for_maf=clinical_anno[unique(maf_for_maf$Tumor_Sample_Barcode),c("PatientCode","HER2_low_status_long_RE","HR_status")]
colnames(clinicalData_for_maf)[1]<-c("Tumor_Sample_Barcode")


laml2 = read.maf(maf = maf_for_maf,
                 clinicalData = clinicalData_for_maf ,
                 verbose = FALSE)


pt.vs.rt <- mafCompare(m1 = laml1, m2 = laml2, m1Name = 'HR_neg_IHC_1_ori', m2Name = 'HR_neg_IHC_2_ori', minMut = 5)
print(pt.vs.rt)
# $results
# Hugo_Symbol HR_neg_IHC_1_ori HR_neg_IHC_2_ori       pval        or     ci.up     ci.low    adjPval
# 1:      PIK3CA                5                7 0.03301101 0.2237408  1.020515 0.04475334 0.06602203
# 2:        TP53               31               11 0.18983801 2.3734372 10.473302 0.53375209 0.18983801
# 
# $SampleSummary
# Cohort SampleSize
# 1: HR_neg_IHC_1_ori         38
# 2: HR_neg_IHC_2_ori         17
sig_dif_names_hr_neg <- data.frame(pt.vs.rt$results)
sig_dif_names_hr_neg <- as.vector(sig_dif_names_hr_neg[sig_dif_names_hr_neg$pval < 0.05,"Hugo_Symbol"])
sig_dif_names_hr_neg
# [1] "PIK3CA"

sig_dif_names <- c(sig_dif_names_hr_pos,sig_dif_names_hr_neg)

temp_gene_order <- rownames(mut_ori_40)[rownames(mut_ori_40) %in% c(rownames(mut_ori_40)[1:10],sig_dif_names)]
gene_order <- temp_gene_order


######################for sample order

maf_for_maf=FUSCC_HER2_low_project_WES_Somatic[FUSCC_HER2_low_project_WES_Somatic$Tumor_Sample_Barcode %in% her2_low_HR_pos_IHC_1_ori,]
genes=read.csv("CAG.csv")[,1]
genes<-genes[genes %in% unique(maf_for_maf$Hugo_Symbol)]
maf_for_maf<-maf_for_maf[maf_for_maf$Hugo_Symbol %in% genes,]

clinicalData_for_maf=clinical_anno[unique(maf_for_maf$Tumor_Sample_Barcode),c("PatientCode","HER2_low_status_long_RE","HR_status")]
colnames(clinicalData_for_maf)[1]<-c("Tumor_Sample_Barcode")


laml = read.maf(maf = maf_for_maf,
                clinicalData = clinicalData_for_maf ,
                verbose = FALSE)

##rename the onco_matrix.txt
oncoplot(colors = vc_cols,maf = laml, draw_titv = TRUE,writeMatrix = T,clinicalFeatures = c("HER2_low_status_long_RE","HR_status"),sortByAnnotation = TRUE,
         annotationOrder=c("0","HER2_1","HER2_2","3"),genes  = gene_order,keepGeneOrder = T,removeNonMutated=F)#writeMatrix可自动产生一个对应的txt文件


file.rename("onco_matrix.txt","Mutation_waterfall_her2_low_HR_pos_IHC_1_ori.txt")

mut_ori <- read.csv("Mutation_waterfall_her2_low_HR_pos_IHC_1_ori.txt",sep = "\t")
her2_low_HR_pos_IHC_1 <- c(colnames(mut_ori),her2_low_HR_pos_IHC_1_ori)
her2_low_HR_pos_IHC_1 <- her2_low_HR_pos_IHC_1[!duplicated(her2_low_HR_pos_IHC_1)]

#################

maf_for_maf=FUSCC_HER2_low_project_WES_Somatic[FUSCC_HER2_low_project_WES_Somatic$Tumor_Sample_Barcode %in% her2_low_HR_pos_IHC_2_ori,]
genes=read.csv("CAG.csv")[,1]
genes<-genes[genes %in% unique(maf_for_maf$Hugo_Symbol)]
maf_for_maf<-maf_for_maf[maf_for_maf$Hugo_Symbol %in% genes,]

clinicalData_for_maf=clinical_anno[unique(maf_for_maf$Tumor_Sample_Barcode),c("PatientCode","HER2_low_status_long_RE","HR_status")]
colnames(clinicalData_for_maf)[1]<-c("Tumor_Sample_Barcode")


laml = read.maf(maf = maf_for_maf,
                clinicalData = clinicalData_for_maf ,
                verbose = FALSE)

##rename the onco_matrix.txt
oncoplot(colors = vc_cols,maf = laml, draw_titv = TRUE,writeMatrix = T,clinicalFeatures = c("HER2_low_status_long_RE","HR_status"),sortByAnnotation = TRUE,
         annotationOrder=c("0","HER2_1","HER2_2","3"),genes  = gene_order,keepGeneOrder = T,removeNonMutated=F)#writeMatrix可自动产生一个对应的txt文件


file.rename("onco_matrix.txt","Mutation_waterfall_her2_low_HR_pos_IHC_2_ori.txt")

mut_ori <- read.csv("Mutation_waterfall_her2_low_HR_pos_IHC_2_ori.txt",sep = "\t")
her2_low_HR_pos_IHC_2 <- c(colnames(mut_ori),her2_low_HR_pos_IHC_2_ori)
her2_low_HR_pos_IHC_2 <- her2_low_HR_pos_IHC_2[!duplicated(her2_low_HR_pos_IHC_2)]

#################

maf_for_maf=FUSCC_HER2_low_project_WES_Somatic[FUSCC_HER2_low_project_WES_Somatic$Tumor_Sample_Barcode %in% her2_low_HR_neg_IHC_1_ori,]
genes=read.csv("CAG.csv")[,1]
genes<-genes[genes %in% unique(maf_for_maf$Hugo_Symbol)]
maf_for_maf<-maf_for_maf[maf_for_maf$Hugo_Symbol %in% genes,]

clinicalData_for_maf=clinical_anno[unique(maf_for_maf$Tumor_Sample_Barcode),c("PatientCode","HER2_low_status_long_RE","HR_status")]
colnames(clinicalData_for_maf)[1]<-c("Tumor_Sample_Barcode")


laml = read.maf(maf = maf_for_maf,
                clinicalData = clinicalData_for_maf ,
                verbose = FALSE)

##rename the onco_matrix.txt
oncoplot(colors = vc_cols,maf = laml, draw_titv = TRUE,writeMatrix = T,clinicalFeatures = c("HER2_low_status_long_RE","HR_status"),sortByAnnotation = TRUE,
         annotationOrder=c("0","HER2_1","HER2_2","3"),genes  = gene_order,keepGeneOrder = T,removeNonMutated=F)#writeMatrix可自动产生一个对应的txt文件


file.rename("onco_matrix.txt","Mutation_waterfall_her2_low_HR_neg_IHC_1_ori.txt")

mut_ori <- read.csv("Mutation_waterfall_her2_low_HR_neg_IHC_1_ori.txt",sep = "\t")
her2_low_HR_neg_IHC_1 <- c(colnames(mut_ori),her2_low_HR_neg_IHC_1_ori)
her2_low_HR_neg_IHC_1 <- her2_low_HR_neg_IHC_1[!duplicated(her2_low_HR_neg_IHC_1)]
#################

maf_for_maf=FUSCC_HER2_low_project_WES_Somatic[FUSCC_HER2_low_project_WES_Somatic$Tumor_Sample_Barcode %in% her2_low_HR_neg_IHC_2_ori,]
genes=read.csv("CAG.csv")[,1]
genes<-genes[genes %in% unique(maf_for_maf$Hugo_Symbol)]
maf_for_maf<-maf_for_maf[maf_for_maf$Hugo_Symbol %in% genes,]

clinicalData_for_maf=clinical_anno[unique(maf_for_maf$Tumor_Sample_Barcode),c("PatientCode","HER2_low_status_long_RE","HR_status")]
colnames(clinicalData_for_maf)[1]<-c("Tumor_Sample_Barcode")


laml = read.maf(maf = maf_for_maf,
                clinicalData = clinicalData_for_maf ,
                verbose = FALSE)

##rename the onco_matrix.txt
oncoplot(colors = vc_cols,maf = laml, draw_titv = TRUE,writeMatrix = T,clinicalFeatures = c("HER2_low_status_long_RE","HR_status"),sortByAnnotation = TRUE,
         annotationOrder=c("0","HER2_1","HER2_2","3"),genes  = gene_order,keepGeneOrder = T,removeNonMutated=F)#writeMatrix可自动产生一个对应的txt文件


file.rename("onco_matrix.txt","Mutation_waterfall_her2_low_HR_neg_IHC_2_ori.txt")

mut_ori <- read.csv("Mutation_waterfall_her2_low_HR_neg_IHC_2_ori.txt",sep = "\t")
her2_low_HR_neg_IHC_2 <- c(colnames(mut_ori),her2_low_HR_neg_IHC_2_ori)
her2_low_HR_neg_IHC_2 <- her2_low_HR_neg_IHC_2[!duplicated(her2_low_HR_neg_IHC_2)]


clinical_anno <- clinical_anno[c(her2_low_HR_pos_IHC_1,her2_low_HR_pos_IHC_2,her2_low_HR_neg_IHC_1,her2_low_HR_neg_IHC_2),]


##blank template
mut_matrix <- matrix("NA" ,ncol = length(rownames(clinical_anno)), nrow = length(gene_order))
colnames(mut_matrix) <- rownames(clinical_anno)
mut_matrix <- as.data.frame(mut_matrix)


mut_ori <- read.csv("Mutation_waterfall_all_HER2_low_40.txt",sep = "\t")[gene_order,]

rownames(mut_matrix) <- rownames(mut_ori)
mut_matrix[,colnames(mut_ori)] <- mut_ori

## anno WT 
mut_matrix[mut_matrix == c("")] <- c("WT")
mut_matrix[mut_matrix == c("0")] <- c("WT")
mut_matrix[mut_matrix == 0] <- c("WT")
mut_matrix[mut_matrix =="NA"] <- NA
#table(mut_matrix)

#### for separate
total_order_cut_col <- factor(paste(clinical_anno$HR_status,clinical_anno$HER2_low_status_long_RE),levels = c("Positive HER2_1","Positive HER2_2","Negative HER2_1","Negative HER2_2"))


##CNV matrix#####################################
## frequent CNV driving genes
cna_list <- unique(c(c("FGFR1","MYC","CCND1","CCNE1","ERBB2","UPD","EGFR","CDK4","IGF1R","AURKA","CDK6"),c("RB1","PTEN","CDKN2A","MLL3","MAP2K4","TP53","ERBB2","PIK3CA")))
cna_list <- cna_list[cna_list %in% rownames(FUSCC_HER2_low_project_GISTICgene.thre)]

sig_amp_genes <- cna_list

show_gene_number <- 5

gene_amp_matrix <- matrix(NA ,ncol = length(rownames(clinical_anno)), nrow = show_gene_number)
colnames(mut_matrix) <- rownames(clinical_anno)

check_list <- her2_low_list
check_list <- check_list[check_list %in% colnames(FUSCC_HER2_low_project_GISTICgene.thre)]
check_list_all <- rownames(clinical_anno)
check_list_all <- check_list_all[check_list_all %in% colnames(FUSCC_HER2_low_project_GISTICgene.thre)]


##sort
TEMP_1 <- rep(NA,length(sig_amp_genes))

for (i in 1:length(TEMP_1)){
  TEMP_1[i] <- table(FUSCC_HER2_low_project_GISTICgene.thre[sig_amp_genes[i],check_list] %in% c(2,"2"))[2]
}
names(TEMP_1) <- sig_amp_genes
sig_amp_genes_count <- (sort(TEMP_1,decreasing = T))

##top genes
top_sig_amp_genes <- names(sig_amp_genes_count)

##record
rownames(gene_amp_matrix) <- top_sig_amp_genes[1:show_gene_number]
gene_amp_matrix <- as.data.frame(gene_amp_matrix)
colnames(gene_amp_matrix) <- rownames(clinical_anno)
gene_amp_matrix[,check_list_all] <- FUSCC_HER2_low_project_GISTICgene.thre[top_sig_amp_genes[1:show_gene_number],check_list_all]




sig_del_genes <- cna_list

##blank template
gene_del_matrix <- matrix(NA ,ncol = length(rownames(clinical_anno)), nrow = show_gene_number)
colnames(mut_matrix) <- rownames(clinical_anno)

check_list <- her2_low_list
check_list <- check_list[check_list %in% colnames(FUSCC_HER2_low_project_GISTICgene.thre)]
check_list_all <- rownames(clinical_anno)
check_list_all <- check_list_all[check_list_all %in% colnames(FUSCC_HER2_low_project_GISTICgene.thre)]


##sort
TEMP_1 <- rep(NA,length(sig_del_genes))

for (i in 1:length(TEMP_1)){
  TEMP_1[i] <- table(FUSCC_HER2_low_project_GISTICgene.thre[sig_del_genes[i],check_list] %in% c(-2,"-2"))[2]
}
names(TEMP_1) <- sig_del_genes
sig_del_genes_count <- (sort(TEMP_1,decreasing = T))

##top genes
top_sig_del_genes <- names(sig_del_genes_count)

##record
rownames(gene_del_matrix) <- top_sig_del_genes[1:show_gene_number]
gene_del_matrix <- as.data.frame(gene_del_matrix)
colnames(gene_del_matrix) <- rownames(clinical_anno)
gene_del_matrix[,check_list_all] <- FUSCC_HER2_low_project_GISTICgene.thre[top_sig_del_genes[1:show_gene_number],check_list_all]



################################################
## RNA matrix


RNA_gene_list <- rownames(FUSCC_HER2_low_project_RNA_seq_log2FPKM)
RNA_group_list <- paste(clinical_anno$HR_status,clinical_anno$HER2_low_status_long_RE)
names(RNA_group_list) <- paste(rownames(clinical_anno),"RNA_T",sep = "_")
RNA_group_list <- RNA_group_list[names(RNA_group_list) %in% colnames(FUSCC_HER2_low_project_RNA_seq_log2FPKM)]

RNA_group_list_1 <- RNA_group_list[RNA_group_list %in% c("Positive HER2_1","Positive HER2_2")]
RNA_group_list_2 <- RNA_group_list[RNA_group_list %in% c("Negative HER2_1","Negative HER2_2")]

RNA_gene_list_record1 <- rep(NA,length(RNA_gene_list))
names(RNA_gene_list_record1) <- RNA_gene_list

RNA_gene_list_record2 <- rep(NA,length(RNA_gene_list))
names(RNA_gene_list_record2) <- RNA_gene_list

for (i in RNA_gene_list){
  temp_for_test <- data.frame(Expr=as.numeric(FUSCC_HER2_low_project_RNA_seq_log2FPKM[i,names(RNA_group_list_1)]),Cluster=RNA_group_list_1) 
  fit<-kruskal.test(Expr~Cluster,data = temp_for_test)
  RNA_gene_list_record1[i] <- fit[["p.value"]]
}



for (i in RNA_gene_list){
  temp_for_test <- data.frame(Expr=as.numeric(FUSCC_HER2_low_project_RNA_seq_log2FPKM[i,names(RNA_group_list_2)]),Cluster=RNA_group_list_2) 
  fit<-kruskal.test(Expr~Cluster,data = temp_for_test)
  RNA_gene_list_record2[i] <- fit[["p.value"]]
}


##direction
RNA_gene_list_record1_dir <-rep(NA,length(RNA_gene_list_record1)) 
names(RNA_gene_list_record1_dir) <- names(RNA_gene_list_record1)
check_temp <- (apply(FUSCC_HER2_low_project_RNA_seq_log2FPKM[,names(RNA_group_list_1)[RNA_group_list_1 %in% "Positive HER2_1"]], 1, median,na.rm=T) > 
                 apply(FUSCC_HER2_low_project_RNA_seq_log2FPKM[,names(RNA_group_list_1)[RNA_group_list_1 %in% "Positive HER2_2"]], 1, median,na.rm=T) )
table(check_temp)
RNA_gene_list_record1_dir[names(check_temp)[check_temp]] <- "Up"
RNA_gene_list_record1_dir[names(check_temp)[!check_temp]] <- "Down"


RNA_gene_list_record2_dir <-rep(NA,length(RNA_gene_list_record2)) 
names(RNA_gene_list_record2_dir) <- names(RNA_gene_list_record2)
check_temp <- (apply(FUSCC_HER2_low_project_RNA_seq_log2FPKM[,names(RNA_group_list_2)[RNA_group_list_2 %in% "Negative HER2_1"]], 1, median,na.rm=T) > 
                 apply(FUSCC_HER2_low_project_RNA_seq_log2FPKM[,names(RNA_group_list_2)[RNA_group_list_2 %in% "Negative HER2_2"]], 1, median,na.rm=T) )
table(check_temp)
RNA_gene_list_record2_dir[names(check_temp)[check_temp]] <- "Up"
RNA_gene_list_record2_dir[names(check_temp)[!check_temp]] <- "Down"

##filter P
RNA_gene_list_record1_FDR <- p.adjust(RNA_gene_list_record1,method = c("fdr"))

RNA_gene_list_record1_FDR_sort <- sort(RNA_gene_list_record1_FDR,decreasing = F)
RNA_gene_list_record1_FDR_sort_sig <- RNA_gene_list_record1_FDR_sort

RNA_gene_list_record1_FDR_order <- c(names(RNA_gene_list_record1_FDR_sort_sig)[names(RNA_gene_list_record1_FDR_sort_sig) %in% names(RNA_gene_list_record1_dir)[RNA_gene_list_record1_dir %in% "Up"]][1:50],
                                     names(RNA_gene_list_record1_FDR_sort_sig)[names(RNA_gene_list_record1_FDR_sort_sig) %in% names(RNA_gene_list_record1_dir)[RNA_gene_list_record1_dir %in% "Down"]][1:50])
RNA_gene_list_record1_FDR_order <- RNA_gene_list_record1_FDR_order[!is.na(RNA_gene_list_record1_FDR_order)]


RNA_gene_list_record2_FDR <- p.adjust(RNA_gene_list_record2,method = c("fdr"))

RNA_gene_list_record2_FDR_sort <- sort(RNA_gene_list_record2_FDR,decreasing = F)
RNA_gene_list_record2_FDR_sort_sig <- RNA_gene_list_record2_FDR_sort

RNA_gene_list_record2_FDR_order <- c(names(RNA_gene_list_record2_FDR_sort_sig)[names(RNA_gene_list_record2_FDR_sort_sig) %in% names(RNA_gene_list_record2_dir)[RNA_gene_list_record2_dir %in% "Up"]][1:50],
                                     names(RNA_gene_list_record2_FDR_sort_sig)[names(RNA_gene_list_record2_FDR_sort_sig) %in% names(RNA_gene_list_record2_dir)[RNA_gene_list_record2_dir %in% "Down"]][1:50])
RNA_gene_list_record2_FDR_order <- RNA_gene_list_record2_FDR_order[!is.na(RNA_gene_list_record2_FDR_order)]


intersect(RNA_gene_list_record1_FDR_order,RNA_gene_list_record2_FDR_order)

RNA_gene_list_record_FDR_order <- c(RNA_gene_list_record1_FDR_order,RNA_gene_list_record2_FDR_order)
RNA_gene_list_record_FDR_order[duplicated(RNA_gene_list_record_FDR_order)] <- paste0(RNA_gene_list_record_FDR_order[duplicated(RNA_gene_list_record_FDR_order)],"_dup")

RNA_row_split <- c(rep("A",50),rep("B",50),rep("C",50),rep("D",50))
RNA_row_split <- RNA_row_split[!duplicated(RNA_gene_list_record_FDR_order)]
table(RNA_row_split)

RNA_matrix_total <- matrix(NA,ncol = length(rownames(clinical_anno)),nrow = length(RNA_gene_list_record_FDR_order))
RNA_matrix_total <- data.frame(RNA_matrix_total)
rownames(RNA_matrix_total) <- RNA_gene_list_record_FDR_order
colnames(RNA_matrix_total) <- rownames(clinical_anno)


RNA_matrix_total[RNA_gene_list_record_FDR_order,colnames(RNA_matrix_total)[paste(colnames(RNA_matrix_total),"_RNA_T",sep = "") %in% colnames(FUSCC_HER2_low_project_RNA_seq_log2FPKM)]] <- FUSCC_HER2_low_project_RNA_seq_log2FPKM[as.character(strsplit(RNA_gene_list_record_FDR_order,"_dup")),paste(colnames(RNA_matrix_total),"_RNA_T",sep = "")[paste(colnames(RNA_matrix_total),"_RNA_T",sep = "") %in% colnames(FUSCC_HER2_low_project_RNA_seq_log2FPKM)]]


gene_anno_RNA <- rbind(data.frame(gene=RNA_gene_list_record_FDR_order[1:100], FDR =RNA_gene_list_record1_FDR_sort[as.character(strsplit(RNA_gene_list_record_FDR_order[1:100],"_dup"))]),
                       data.frame(gene=RNA_gene_list_record_FDR_order[101:200], FDR =RNA_gene_list_record2_FDR_sort[as.character(strsplit(RNA_gene_list_record_FDR_order[101:200],"_dup"))]))
gene_anno_RNA <- gene_anno_RNA[!duplicated(gene_anno_RNA$gene),]



check_temp <- gene_anno_RNA[,"FDR"]<0.1
gene_anno_RNA[check_temp,"FDR"] <- "<0.1"
gene_anno_RNA[!check_temp,"FDR"] <- ">=0.1"




################################################
## PRO matrix
PRO_gene_list <- rownames(FUSCC_HER2_low_project_TMT_PRO_normalized)
PRO_group_list <- paste(clinical_anno$HR_status,clinical_anno$HER2_low_status_long_RE)
names(PRO_group_list) <- paste(rownames(clinical_anno),"PRO_T",sep = "_")
PRO_group_list <- PRO_group_list[names(PRO_group_list) %in% colnames(FUSCC_HER2_low_project_TMT_PRO_normalized)]

PRO_group_list_1 <- PRO_group_list[PRO_group_list %in% c("Positive HER2_1","Positive HER2_2")]
PRO_group_list_2 <- PRO_group_list[PRO_group_list %in% c("Negative HER2_1","Negative HER2_2")]

PRO_gene_list_record1 <- rep(NA,length(PRO_gene_list))
names(PRO_gene_list_record1) <- PRO_gene_list

PRO_gene_list_record2 <- rep(NA,length(PRO_gene_list))
names(PRO_gene_list_record2) <- PRO_gene_list

for (i in PRO_gene_list){
  temp_for_test <- data.frame(Expr=as.numeric(FUSCC_HER2_low_project_TMT_PRO_normalized[i,names(PRO_group_list_1)]),Cluster=PRO_group_list_1) 
  fit<-kruskal.test(Expr~Cluster,data = temp_for_test)
  PRO_gene_list_record1[i] <- fit[["p.value"]]
}

for (i in PRO_gene_list){
  temp_for_test <- data.frame(Expr=as.numeric(FUSCC_HER2_low_project_TMT_PRO_normalized[i,names(PRO_group_list_2)]),Cluster=PRO_group_list_2) 
  fit<-kruskal.test(Expr~Cluster,data = temp_for_test)
  PRO_gene_list_record2[i] <- fit[["p.value"]]
}



##direction
PRO_gene_list_record1_dir <-rep(NA,length(PRO_gene_list_record1)) 
names(PRO_gene_list_record1_dir) <- names(PRO_gene_list_record1)
check_temp <- (apply(FUSCC_HER2_low_project_TMT_PRO_normalized[,names(PRO_group_list_1)[PRO_group_list_1 %in% "Positive HER2_1"]], 1, median,na.rm=T) > 
                 apply(FUSCC_HER2_low_project_TMT_PRO_normalized[,names(PRO_group_list_1)[PRO_group_list_1 %in% "Positive HER2_2"]], 1, median,na.rm=T) )
table(check_temp)
PRO_gene_list_record1_dir[names(check_temp)[check_temp]] <- "Up"
PRO_gene_list_record1_dir[names(check_temp)[!check_temp]] <- "Down"


PRO_gene_list_record2_dir <-rep(NA,length(PRO_gene_list_record2)) 
names(PRO_gene_list_record2_dir) <- names(PRO_gene_list_record2)
check_temp <- (apply(FUSCC_HER2_low_project_TMT_PRO_normalized[,names(PRO_group_list_2)[PRO_group_list_2 %in% "Negative HER2_1"]], 1, median,na.rm=T) > 
                 apply(FUSCC_HER2_low_project_TMT_PRO_normalized[,names(PRO_group_list_2)[PRO_group_list_2 %in% "Negative HER2_2"]], 1, median,na.rm=T) )
table(check_temp)
PRO_gene_list_record2_dir[names(check_temp)[check_temp]] <- "Up"
PRO_gene_list_record2_dir[names(check_temp)[!check_temp]] <- "Down"

##filter P
PRO_gene_list_record1_FDR <- p.adjust(PRO_gene_list_record1,method = c("fdr"))

PRO_gene_list_record1_FDR_sort <- sort(PRO_gene_list_record1_FDR,decreasing = F)
PRO_gene_list_record1_FDR_sort_sig <- PRO_gene_list_record1_FDR_sort

PRO_gene_list_record1_FDR_order <- c(names(PRO_gene_list_record1_FDR_sort_sig)[names(PRO_gene_list_record1_FDR_sort_sig) %in% names(PRO_gene_list_record2_dir)[PRO_gene_list_record2_dir %in% "Up"]][1:50],
                                     names(PRO_gene_list_record1_FDR_sort_sig)[names(PRO_gene_list_record1_FDR_sort_sig) %in% names(PRO_gene_list_record2_dir)[PRO_gene_list_record2_dir %in% "Down"]][1:50])
PRO_gene_list_record1_FDR_order <- PRO_gene_list_record1_FDR_order[!is.na(PRO_gene_list_record1_FDR_order)]


PRO_gene_list_record2_FDR <- p.adjust(PRO_gene_list_record2,method = c("fdr"))

PRO_gene_list_record2_FDR_sort <- sort(PRO_gene_list_record2_FDR,decreasing = F)
PRO_gene_list_record2_FDR_sort_sig <- PRO_gene_list_record2_FDR_sort

PRO_gene_list_record2_FDR_order <- c(names(PRO_gene_list_record2_FDR_sort_sig)[names(PRO_gene_list_record2_FDR_sort_sig) %in% names(PRO_gene_list_record2_dir)[PRO_gene_list_record2_dir %in% "Up"]][1:50],
                                     names(PRO_gene_list_record2_FDR_sort_sig)[names(PRO_gene_list_record2_FDR_sort_sig) %in% names(PRO_gene_list_record2_dir)[PRO_gene_list_record2_dir %in% "Down"]][1:50])
PRO_gene_list_record2_FDR_order <- PRO_gene_list_record2_FDR_order[!is.na(PRO_gene_list_record2_FDR_order)]


intersect(PRO_gene_list_record1_FDR_order,PRO_gene_list_record2_FDR_order)
# [1] "SMC1A" "FASN"  "ACACA"

PRO_gene_list_record_FDR_order <- c(PRO_gene_list_record1_FDR_order,PRO_gene_list_record2_FDR_order)
PRO_gene_list_record_FDR_order[duplicated(PRO_gene_list_record_FDR_order)] <- paste0(PRO_gene_list_record_FDR_order[duplicated(PRO_gene_list_record_FDR_order)],"_dup")

PRO_row_split <- c(rep("A",50),rep("B",50),rep("C",50),rep("D",50))
PRO_row_split <- PRO_row_split[!duplicated(PRO_gene_list_record_FDR_order)]
table(PRO_row_split)

PRO_matrix_total <- matrix(NA,ncol = length(rownames(clinical_anno)),nrow = length(PRO_gene_list_record_FDR_order))
PRO_matrix_total <- data.frame(PRO_matrix_total)
rownames(PRO_matrix_total) <- PRO_gene_list_record_FDR_order
colnames(PRO_matrix_total) <- rownames(clinical_anno)


PRO_matrix_total[PRO_gene_list_record_FDR_order,colnames(PRO_matrix_total)[paste(colnames(PRO_matrix_total),"_PRO_T",sep = "") %in% colnames(FUSCC_HER2_low_project_TMT_PRO_normalized)]] <- FUSCC_HER2_low_project_TMT_PRO_normalized[as.character(strsplit(PRO_gene_list_record_FDR_order,"_dup")),paste(colnames(PRO_matrix_total),"_PRO_T",sep = "")[paste(colnames(PRO_matrix_total),"_PRO_T",sep = "") %in% colnames(FUSCC_HER2_low_project_TMT_PRO_normalized)]]

gene_anno_PRO <- rbind(data.frame(gene=PRO_gene_list_record_FDR_order[1:100], FDR =PRO_gene_list_record1_FDR_sort[as.character(strsplit(PRO_gene_list_record_FDR_order[1:100],"_dup"))]),
                       data.frame(gene=PRO_gene_list_record_FDR_order[101:200], FDR =PRO_gene_list_record2_FDR_sort[as.character(strsplit(PRO_gene_list_record_FDR_order[101:200],"_dup"))]))
gene_anno_PRO <- gene_anno_PRO[!duplicated(gene_anno_PRO$gene),]


check_temp <- gene_anno_PRO[,"FDR"]<0.1
gene_anno_PRO[check_temp,"FDR"] <- "<0.1"
gene_anno_PRO[!check_temp,"FDR"] <- ">=0.1"


################################################
## POL matrix


POL_gene_list <- rownames(FUSCC_HER2_low_project_metabolite_polar)
POL_group_list <- paste(clinical_anno$HR_status,clinical_anno$HER2_low_status_long_RE)
names(POL_group_list) <- paste(rownames(clinical_anno),"pol_T",sep = "_")
POL_group_list <- POL_group_list[names(POL_group_list) %in% colnames(FUSCC_HER2_low_project_metabolite_polar)]

POL_group_list_1 <- POL_group_list[POL_group_list %in% c("Positive HER2_1","Positive HER2_2")]
POL_group_list_2 <- POL_group_list[POL_group_list %in% c("Negative HER2_1","Negative HER2_2")]

POL_gene_list_record1 <- rep(NA,length(POL_gene_list))
names(POL_gene_list_record1) <- POL_gene_list

POL_gene_list_record2 <- rep(NA,length(POL_gene_list))
names(POL_gene_list_record2) <- POL_gene_list


for (i in POL_gene_list){
  temp_for_test <- data.frame(Expr=as.numeric(FUSCC_HER2_low_project_metabolite_polar[i,names(POL_group_list_1)]),Cluster=POL_group_list_1) 
  fit<-kruskal.test(Expr~Cluster,data = temp_for_test)
  POL_gene_list_record1[i] <- fit[["p.value"]]
}



for (i in POL_gene_list){
  temp_for_test <- data.frame(Expr=as.numeric(FUSCC_HER2_low_project_metabolite_polar[i,names(POL_group_list_2)]),Cluster=POL_group_list_2) 
  fit<-kruskal.test(Expr~Cluster,data = temp_for_test)
  POL_gene_list_record2[i] <- fit[["p.value"]]
}

##direction
POL_gene_list_record1_dir <-rep(NA,length(POL_gene_list_record1)) 
names(POL_gene_list_record1_dir) <- names(POL_gene_list_record1)
check_temp <- (apply(FUSCC_HER2_low_project_metabolite_polar[,names(POL_group_list_1)[POL_group_list_1 %in% "Positive HER2_1"]], 1, median,na.rm=T) > 
                 apply(FUSCC_HER2_low_project_metabolite_polar[,names(POL_group_list_1)[POL_group_list_1 %in% "Positive HER2_2"]], 1, median,na.rm=T) )
table(check_temp)
POL_gene_list_record1_dir[names(check_temp)[check_temp]] <- "Up"
POL_gene_list_record1_dir[names(check_temp)[!check_temp]] <- "Down"


POL_gene_list_record2_dir <-rep(NA,length(POL_gene_list_record2)) 
names(POL_gene_list_record2_dir) <- names(POL_gene_list_record2)
check_temp <- (apply(FUSCC_HER2_low_project_metabolite_polar[,names(POL_group_list_2)[POL_group_list_2 %in% "Negative HER2_1"]], 1, median,na.rm=T) > 
                 apply(FUSCC_HER2_low_project_metabolite_polar[,names(POL_group_list_2)[POL_group_list_2 %in% "Negative HER2_2"]], 1, median,na.rm=T) )
table(check_temp)
POL_gene_list_record2_dir[names(check_temp)[check_temp]] <- "Up"
POL_gene_list_record2_dir[names(check_temp)[!check_temp]] <- "Down"

##filter P
POL_gene_list_record1_FDR <- p.adjust(POL_gene_list_record1,method = c("fdr"))

POL_gene_list_record1_FDR_sort <- sort(POL_gene_list_record1_FDR,decreasing = F)
POL_gene_list_record1_FDR_sort_sig <- POL_gene_list_record1_FDR_sort

POL_gene_list_record1_FDR_order <- c(names(POL_gene_list_record1_FDR_sort_sig)[names(POL_gene_list_record1_FDR_sort_sig) %in% names(POL_gene_list_record2_dir)[POL_gene_list_record2_dir %in% "Up"]][1:25],
                                     names(POL_gene_list_record1_FDR_sort_sig)[names(POL_gene_list_record1_FDR_sort_sig) %in% names(POL_gene_list_record2_dir)[POL_gene_list_record2_dir %in% "Down"]][1:25])
POL_gene_list_record1_FDR_order <- POL_gene_list_record1_FDR_order[!is.na(POL_gene_list_record1_FDR_order)]


POL_gene_list_record2_FDR <- p.adjust(POL_gene_list_record2,method = c("fdr"))

POL_gene_list_record2_FDR_sort <- sort(POL_gene_list_record2_FDR,decreasing = F)
POL_gene_list_record2_FDR_sort_sig <- POL_gene_list_record2_FDR_sort

POL_gene_list_record2_FDR_order <- c(names(POL_gene_list_record2_FDR_sort_sig)[names(POL_gene_list_record2_FDR_sort_sig) %in% names(POL_gene_list_record2_dir)[POL_gene_list_record2_dir %in% "Up"]][1:25],
                                     names(POL_gene_list_record2_FDR_sort_sig)[names(POL_gene_list_record2_FDR_sort_sig) %in% names(POL_gene_list_record2_dir)[POL_gene_list_record2_dir %in% "Down"]][1:25])
POL_gene_list_record2_FDR_order <- POL_gene_list_record2_FDR_order[!is.na(POL_gene_list_record2_FDR_order)]


intersect(POL_gene_list_record1_FDR_order,POL_gene_list_record2_FDR_order)
# [1] "M101T93_NEG"  "M133T340_POS" "M226T56_POS"  "M219T57_POS"  "M253T50_POS"  "M341T388_NEG" "M147T278_NEG" "M276T437_POS" "M191T65_POS"

POL_gene_list_record_FDR_order <- c(POL_gene_list_record1_FDR_order,POL_gene_list_record2_FDR_order)
POL_gene_list_record_FDR_order[duplicated(POL_gene_list_record_FDR_order)] <- paste(POL_gene_list_record_FDR_order[duplicated(POL_gene_list_record_FDR_order)],"_dup")


POL_row_split <- c(rep("A",25),rep("B",25),rep("C",25),rep("D",25))
POL_row_split <- POL_row_split[!duplicated(POL_gene_list_record_FDR_order)]
table(POL_row_split)

POL_matrix_total <- matrix(NA,ncol = length(rownames(clinical_anno)),nrow = length(POL_gene_list_record_FDR_order))
POL_matrix_total <- data.frame(POL_matrix_total)
rownames(POL_matrix_total) <- POL_gene_list_record_FDR_order
colnames(POL_matrix_total) <- rownames(clinical_anno)


POL_matrix_total[POL_gene_list_record_FDR_order,colnames(POL_matrix_total)[paste(colnames(POL_matrix_total),"_pol_T",sep = "") %in% colnames(FUSCC_HER2_low_project_metabolite_polar)]] <- FUSCC_HER2_low_project_metabolite_polar[trimws(as.character(strsplit(POL_gene_list_record_FDR_order,"_dup"))),paste(colnames(POL_matrix_total),"_pol_T",sep = "")[paste(colnames(POL_matrix_total),"_pol_T",sep = "") %in% colnames(FUSCC_HER2_low_project_metabolite_polar)]]

gene_anno_POL <- rbind(data.frame(gene=POL_gene_list_record_FDR_order[1:50], FDR =POL_gene_list_record1_FDR_sort[trimws(as.character(strsplit(POL_gene_list_record_FDR_order[1:50],"_dup")))]),
                       data.frame(gene=POL_gene_list_record_FDR_order[51:100], FDR =POL_gene_list_record2_FDR_sort[trimws(as.character(strsplit(POL_gene_list_record_FDR_order[51:100],"_dup")))]))
gene_anno_POL <- gene_anno_POL[!duplicated(gene_anno_POL$gene),]

check_temp <- gene_anno_POL[,"FDR"]<0.1
gene_anno_POL[check_temp,"FDR"] <- "<0.1"
gene_anno_POL[!check_temp,"FDR"] <- ">=0.1"

################################################
## LIP matrix


LIP_gene_list <- rownames(FUSCC_HER2_low_project_metabolite_lipid)
LIP_group_list <- paste(clinical_anno$HR_status,clinical_anno$HER2_low_status_long_RE)
names(LIP_group_list) <- paste(rownames(clinical_anno),"pol_T",sep = "_")
LIP_group_list <- LIP_group_list[names(LIP_group_list) %in% colnames(FUSCC_HER2_low_project_metabolite_lipid)]

LIP_group_list_1 <- LIP_group_list[LIP_group_list %in% c("Positive HER2_1","Positive HER2_2")]
LIP_group_list_2 <- LIP_group_list[LIP_group_list %in% c("Negative HER2_1","Negative HER2_2")]

LIP_gene_list_record1 <- rep(NA,length(LIP_gene_list))
names(LIP_gene_list_record1) <- LIP_gene_list

LIP_gene_list_record2 <- rep(NA,length(LIP_gene_list))
names(LIP_gene_list_record2) <- LIP_gene_list


for (i in LIP_gene_list){
  temp_for_test <- data.frame(Expr=as.numeric(FUSCC_HER2_low_project_metabolite_lipid[i,names(LIP_group_list_1)]),Cluster=LIP_group_list_1) 
  fit<-kruskal.test(Expr~Cluster,data = temp_for_test)
  LIP_gene_list_record1[i] <- fit[["p.value"]]
}



for (i in LIP_gene_list){
  temp_for_test <- data.frame(Expr=as.numeric(FUSCC_HER2_low_project_metabolite_lipid[i,names(LIP_group_list_2)]),Cluster=LIP_group_list_2) 
  fit<-kruskal.test(Expr~Cluster,data = temp_for_test)
  LIP_gene_list_record2[i] <- fit[["p.value"]]
}

##direction
LIP_gene_list_record1_dir <-rep(NA,length(LIP_gene_list_record1)) 
names(LIP_gene_list_record1_dir) <- names(LIP_gene_list_record1)
check_temp <- (apply(FUSCC_HER2_low_project_metabolite_lipid[,names(LIP_group_list_1)[LIP_group_list_1 %in% "Positive HER2_1"]], 1, median,na.rm=T) > 
                 apply(FUSCC_HER2_low_project_metabolite_lipid[,names(LIP_group_list_1)[LIP_group_list_1 %in% "Positive HER2_2"]], 1, median,na.rm=T) )
table(check_temp)
LIP_gene_list_record1_dir[names(check_temp)[check_temp]] <- "Up"
LIP_gene_list_record1_dir[names(check_temp)[!check_temp]] <- "Down"


LIP_gene_list_record2_dir <-rep(NA,length(LIP_gene_list_record2)) 
names(LIP_gene_list_record2_dir) <- names(LIP_gene_list_record2)
check_temp <- (apply(FUSCC_HER2_low_project_metabolite_lipid[,names(LIP_group_list_2)[LIP_group_list_2 %in% "Negative HER2_1"]], 1, median,na.rm=T) > 
                 apply(FUSCC_HER2_low_project_metabolite_lipid[,names(LIP_group_list_2)[LIP_group_list_2 %in% "Negative HER2_2"]], 1, median,na.rm=T) )
table(check_temp)
LIP_gene_list_record2_dir[names(check_temp)[check_temp]] <- "Up"
LIP_gene_list_record2_dir[names(check_temp)[!check_temp]] <- "Down"

##filter P
LIP_gene_list_record1_FDR <- p.adjust(LIP_gene_list_record1,method = c("fdr"))
LIP_gene_list_record1_FDR_sort <- sort(LIP_gene_list_record1_FDR,decreasing = F)
LIP_gene_list_record1_FDR_sort_sig <- LIP_gene_list_record1_FDR_sort

LIP_gene_list_record1_FDR_order <- c(names(LIP_gene_list_record1_FDR_sort_sig)[names(LIP_gene_list_record1_FDR_sort_sig) %in% names(LIP_gene_list_record2_dir)[LIP_gene_list_record2_dir %in% "Up"]][1:25],
                                     names(LIP_gene_list_record1_FDR_sort_sig)[names(LIP_gene_list_record1_FDR_sort_sig) %in% names(LIP_gene_list_record2_dir)[LIP_gene_list_record2_dir %in% "Down"]][1:25])
LIP_gene_list_record1_FDR_order <- LIP_gene_list_record1_FDR_order[!is.na(LIP_gene_list_record1_FDR_order)]


LIP_gene_list_record2_FDR <- p.adjust(LIP_gene_list_record2,method = c("fdr"))

LIP_gene_list_record2_FDR_sort <- sort(LIP_gene_list_record2_FDR,decreasing = F)
LIP_gene_list_record2_FDR_sort_sig <- LIP_gene_list_record2_FDR_sort

LIP_gene_list_record2_FDR_order <- c(names(LIP_gene_list_record2_FDR_sort_sig)[names(LIP_gene_list_record2_FDR_sort_sig) %in% names(LIP_gene_list_record2_dir)[LIP_gene_list_record2_dir %in% "Up"]][1:25],
                                     names(LIP_gene_list_record2_FDR_sort_sig)[names(LIP_gene_list_record2_FDR_sort_sig) %in% names(LIP_gene_list_record2_dir)[LIP_gene_list_record2_dir %in% "Down"]][1:25])
LIP_gene_list_record2_FDR_order <- LIP_gene_list_record2_FDR_order[!is.na(LIP_gene_list_record2_FDR_order)]


intersect(LIP_gene_list_record1_FDR_order,LIP_gene_list_record2_FDR_order)
# [1] "M831T602_POS"  "M803T573_POS"  "M1472T437_NEG" "M879T533_NEG"  "M657T565_POS"

LIP_gene_list_record_FDR_order <- c(LIP_gene_list_record1_FDR_order,LIP_gene_list_record2_FDR_order)
LIP_gene_list_record_FDR_order[duplicated(LIP_gene_list_record_FDR_order)] <- paste(LIP_gene_list_record_FDR_order[duplicated(LIP_gene_list_record_FDR_order)],"_dup")

LIP_row_split <- c(rep("A",25),rep("B",25),rep("C",25),rep("D",25))
LIP_row_split <- LIP_row_split[!duplicated(LIP_gene_list_record_FDR_order)]
table(LIP_row_split)

LIP_matrix_total <- matrix(NA,ncol = length(rownames(clinical_anno)),nrow = length(LIP_gene_list_record_FDR_order))
LIP_matrix_total <- data.frame(LIP_matrix_total)
rownames(LIP_matrix_total) <- LIP_gene_list_record_FDR_order
colnames(LIP_matrix_total) <- rownames(clinical_anno)


LIP_matrix_total[LIP_gene_list_record_FDR_order,colnames(LIP_matrix_total)[paste(colnames(LIP_matrix_total),"_pol_T",sep = "") %in% colnames(FUSCC_HER2_low_project_metabolite_lipid)]] <- FUSCC_HER2_low_project_metabolite_lipid[trimws(as.character(strsplit(LIP_gene_list_record_FDR_order,"_dup"))),paste(colnames(LIP_matrix_total),"_pol_T",sep = "")[paste(colnames(LIP_matrix_total),"_pol_T",sep = "") %in% colnames(FUSCC_HER2_low_project_metabolite_lipid)]]
gene_anno_LIP <- rbind(data.frame(gene=LIP_gene_list_record_FDR_order[1:50], FDR =LIP_gene_list_record1_FDR_sort[trimws(as.character(strsplit(LIP_gene_list_record_FDR_order[1:50],"_dup")))]),
                       data.frame(gene=LIP_gene_list_record_FDR_order[51:100], FDR =LIP_gene_list_record2_FDR_sort[trimws(as.character(strsplit(LIP_gene_list_record_FDR_order[51:100],"_dup")))]))
gene_anno_LIP <- gene_anno_LIP[!duplicated(gene_anno_LIP$gene),]

check_temp <- gene_anno_LIP[,"FDR"]<0.1
gene_anno_LIP[check_temp,"FDR"] <- "<0.1"
gene_anno_LIP[!check_temp,"FDR"] <- ">=0.1"



