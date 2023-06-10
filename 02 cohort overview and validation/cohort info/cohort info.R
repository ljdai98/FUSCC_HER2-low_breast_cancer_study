table(FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE)
# HER2_0      HER2_low HER2_positive 
# 91           434           182 

##omics number 
her2_low_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_low"),"PatientCode"]

## WES -----------------
table(substr(unique(FUSCC_HER2_low_project_WES_Somatic$Tumor_Sample_Barcode),1,4) %in% her2_low_list)
# FALSE  TRUE 
# 231   374

table(substr((FUSCC_HER2_low_project_WES_Somatic$Tumor_Sample_Barcode),1,4) %in% her2_low_list)
# FALSE  TRUE 
# 37304 49122 

## GISTIC copy number -----------------
table(substr(colnames(FUSCC_HER2_low_project_GISTICgene.alldata),1,4) %in% her2_low_list)
# FALSE  TRUE 
# 245   379  

dim(FUSCC_HER2_low_project_GISTICgene.alldata)
# [1] 27100   624

dim(FUSCC_HER2_low_project_GISTICpeaks.val)
# [1]  76 624

## RNA_seq -----------------
table(substr(colnames(FUSCC_HER2_low_project_RNA_seq_log2FPKM),1,4) %in% her2_low_list)
# FALSE  TRUE 
# 269   421

dim(FUSCC_HER2_low_project_RNA_seq_log2FPKM)
# [1] 19892   690

## TMT Proteome -----------------
table(substr(colnames(FUSCC_HER2_low_project_TMT_PRO_normalized),1,4) %in% her2_low_list)
# FALSE  TRUE 
# 98   156

dim(FUSCC_HER2_low_project_TMT_PRO_normalized)
# [1] 7952  254

## Metabolomics --------------
table(substr(colnames(FUSCC_HER2_low_project_metabolite_polar),1,4) %in% her2_low_list)
# FALSE  TRUE 
# 169   254 

dim(FUSCC_HER2_low_project_metabolite_polar)
# [1] 669 423

dim(FUSCC_HER2_low_project_metabolite_lipid)
# [1] 1312  423
