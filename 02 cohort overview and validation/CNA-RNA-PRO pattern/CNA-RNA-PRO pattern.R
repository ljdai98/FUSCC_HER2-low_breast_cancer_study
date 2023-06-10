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

##group_list
clinical_anno<- FUSCC_HER2_low_project_Cohort.Info
her2_low_HR_pos_IHC_1_ori <- clinical_anno[clinical_anno$HR_status %in% c("Positive") & clinical_anno$HER2_low_status_long_RE %in% c("HER2_1"),"PatientCode"]
her2_low_HR_pos_IHC_2_ori <- clinical_anno[clinical_anno$HR_status %in% c("Positive") & clinical_anno$HER2_low_status_long_RE %in% c("HER2_2"),"PatientCode"]
her2_low_HR_neg_IHC_1_ori <- clinical_anno[clinical_anno$HR_status %in% c("Negative") & clinical_anno$HER2_low_status_long_RE %in% c("HER2_1"),"PatientCode"]
her2_low_HR_neg_IHC_2_ori <- clinical_anno[clinical_anno$HR_status %in% c("Negative") & clinical_anno$HER2_low_status_long_RE %in% c("HER2_2"),"PatientCode"]

group_list_combo <- c(rep("HR_pos_IHC_1",length(her2_low_HR_pos_IHC_1_ori)),
                      rep("HR_pos_IHC_2",length(her2_low_HR_pos_IHC_2_ori)),
                      rep("HR_neg_IHC_1",length(her2_low_HR_neg_IHC_1_ori)),
                      rep("HR_neg_IHC_2",length(her2_low_HR_neg_IHC_2_ori)))
names(group_list_combo) <-  c(her2_low_HR_pos_IHC_1_ori,her2_low_HR_pos_IHC_2_ori,her2_low_HR_neg_IHC_1_ori,her2_low_HR_neg_IHC_2_ori)
group_list_combo <- factor(group_list_combo,levels = c("HR_pos_IHC_1","HR_pos_IHC_2","HR_neg_IHC_1","HR_neg_IHC_2"))

#####################
## CNA
#####################
cna_list <- unique(c(c("FGFR1","MYC","CCND1","CCNE1","ERBB2","UPD","EGFR","CDK4","IGF1R","AURKA","CDK6"),c("RB1","PTEN","CDKN2A","MLL3","MAP2K4","TP53","ERBB2","PIK3CA")))

#########################
##dataframes
########################
cna_list_with_sig <- cna_list
group_list_combo2 <- group_list_combo

## included sample and genes
table(group_list_combo2)
# HR_pos_IHC_1 HR_pos_IHC_2 HR_neg_IHC_1 HR_neg_IHC_2 
# 179          182           48           25

group_list_combo2 <- group_list_combo2[names(group_list_combo2) %in% colnames(FUSCC_HER2_low_project_GISTICgene.thre)]
group_list_combo2 <- group_list_combo2[paste0(names(group_list_combo2),"_RNA_T") %in% colnames(FUSCC_HER2_low_project_RNA_seq_log2FPKM)]
group_list_combo2 <- group_list_combo2[paste0(names(group_list_combo2),"_PRO_T") %in% colnames(FUSCC_HER2_low_project_TMT_PRO_normalized)]
table(group_list_combo2)
# HR_pos_IHC_1 HR_pos_IHC_2 HR_neg_IHC_1 HR_neg_IHC_2 
# 60           51           21            8

length(cna_list_with_sig)
# [1] 18
cna_list_with_sig2 <- cna_list_with_sig
cna_list_with_sig2 <- cna_list_with_sig2[cna_list_with_sig2 %in% rownames(FUSCC_HER2_low_project_RNA_seq_log2FPKM)]
# cna_list_with_sig2 <- cna_list_with_sig2[cna_list_with_sig2 %in% rownames(CBCGA.Extended_PRO_expanded_NA_not_filtered)]
cna_list_with_sig2 
# [1] "FGFR1"  "MYC"    "CCND1"  "CCNE1"  "ERBB2"  "EGFR"   "CDK4"   "IGF1R"  "AURKA"  "CDK6"   "RB1"    "PTEN"   "CDKN2A" "MAP2K4" "TP53"  
# [16] "PIK3CA"
length(cna_list_with_sig2)
# [1] 19

##reorder
cna_list_with_sig3 <- c("ERBB2","EGFR","FGFR1","IGF1R","CDK4","CDK6","CCND1","CCNE1","CDKN2A","AURKA","RB1","MYC","PIK3CA","MAP2K4","PTEN","TP53")


##CNA_data.frame
CNA_data.frame <- matrix(NA,ncol=(1+length(group_list_combo2)),nrow=length(cna_list_with_sig3))
CNA_data.frame <- data.frame(CNA_data.frame)
rownames(CNA_data.frame) <- cna_list_with_sig3
colnames(CNA_data.frame) <- c("Index",names(group_list_combo2))
CNA_data.frame[rownames(CNA_data.frame)[rownames(CNA_data.frame) %in% rownames(FUSCC_HER2_low_project_GISTICgene.thre)],names(group_list_combo2)] <- FUSCC_HER2_low_project_GISTICgene.thre[rownames(CNA_data.frame)[rownames(CNA_data.frame) %in% rownames(FUSCC_HER2_low_project_GISTICgene.thre)],names(group_list_combo2)]
rownames(CNA_data.frame) <- paste0("CNA_",rownames(CNA_data.frame))
CNA_data.frame$Index <- 1:length(cna_list_with_sig3)
CNA_data.frame$Index <- CNA_data.frame$Index *3-2



##RNA_data.frame
RNA_data.frame <- matrix(NA,ncol=(1+length(group_list_combo2)),nrow=length(cna_list_with_sig3))
RNA_data.frame <- data.frame(RNA_data.frame)
rownames(RNA_data.frame) <- cna_list_with_sig3
colnames(RNA_data.frame) <- c("Index",names(group_list_combo2))

FUSCC_HER2_low_project_RNA_seq_log2FPKM_scale <- t( scale( t( FUSCC_HER2_low_project_RNA_seq_log2FPKM) ) )
FUSCC_HER2_low_project_RNA_seq_log2FPKM_scale[FUSCC_HER2_low_project_RNA_seq_log2FPKM_scale > 2] = 2
FUSCC_HER2_low_project_RNA_seq_log2FPKM_scale[FUSCC_HER2_low_project_RNA_seq_log2FPKM_scale < -2] = -2
RNA_data.frame[rownames(RNA_data.frame)[rownames(RNA_data.frame) %in% rownames(FUSCC_HER2_low_project_RNA_seq_log2FPKM_scale)],names(group_list_combo2)] <- FUSCC_HER2_low_project_RNA_seq_log2FPKM_scale[rownames(RNA_data.frame)[rownames(RNA_data.frame) %in% rownames(FUSCC_HER2_low_project_RNA_seq_log2FPKM_scale)],paste0(names(group_list_combo2),"_RNA_T")]

rownames(RNA_data.frame) <- paste0("RNA_",rownames(RNA_data.frame))


RNA_data.frame$Index <- 1:length(cna_list_with_sig3)
RNA_data.frame$Index <- RNA_data.frame$Index *3-1

##PRO_data.frame
PRO_data.frame <- matrix(NA,ncol=(1+length(group_list_combo2)),nrow=length(cna_list_with_sig3))
PRO_data.frame <- data.frame(PRO_data.frame)
rownames(PRO_data.frame) <- cna_list_with_sig3
colnames(PRO_data.frame) <- c("Index",names(group_list_combo2))

FUSCC_HER2_low_project_TMT_PRO_normalized_scale <- t( scale( t( FUSCC_HER2_low_project_TMT_PRO_normalized) ) )
FUSCC_HER2_low_project_TMT_PRO_normalized_scale[FUSCC_HER2_low_project_TMT_PRO_normalized_scale > 2] = 2
FUSCC_HER2_low_project_TMT_PRO_normalized_scale[FUSCC_HER2_low_project_TMT_PRO_normalized_scale < -2] = -2
PRO_data.frame[rownames(PRO_data.frame)[rownames(PRO_data.frame) %in% rownames(FUSCC_HER2_low_project_TMT_PRO_normalized_scale)],names(group_list_combo2)] <- FUSCC_HER2_low_project_TMT_PRO_normalized_scale[rownames(PRO_data.frame)[rownames(PRO_data.frame) %in% rownames(FUSCC_HER2_low_project_TMT_PRO_normalized_scale)],paste0(names(group_list_combo2),"_PRO_T")]

rownames(PRO_data.frame) <- paste0("PRO_",rownames(PRO_data.frame))


PRO_data.frame$Index <- 1:length(cna_list_with_sig3)
PRO_data.frame$Index <- PRO_data.frame$Index *3

##Commo_data.frame
Commo_data.frame <- rbind(CNA_data.frame,RNA_data.frame,PRO_data.frame)
Commo_data.frame <- Commo_data.frame[order(Commo_data.frame$Index,decreasing = F),]
Commo_data.frame$Index <- NULL

library(ComplexHeatmap)
##anno col#################################
df = FUSCC_HER2_low_project_Cohort.Info[names(group_list_combo2),c("HR_status","HER2_low_status_long_RE","PAM50")]
colnames(df) <- c("HR status","HER2 IHC score","PAM50")

df[df=="Her2"] <- "HER2"

ha = HeatmapAnnotation(df = df,
                       col  = list(`HR status` = c("Positive" = "#4d8ee4", "Negative" = "#e4574d"),
                                   `HER2 IHC score` = c("HER2_1" = "#7AB801","HER2_2" = "#F2AF01"),
                                   `PAM50`=c("LumA" = "#000286", "LumB" = "#87CCFA","HER2" = "#FF67B8", "Basal" = "#FE0302","Normal"="#36E200")),
                       na_col = "#EAEAEA")
ha
draw(ha)

# pdf("overall heatmap CNA-RNA-PRO.pdf",width =20 ,height = 10)
# ht_mut = Heatmap(as.matrix(Commo_data.frame),cluster_rows = F,cluster_columns = F,
#                  border_gp = gpar(col = "Grey", lwd = 1.5),
#                  rect_gp = gpar(col = 'white', lwd = 1),
#                  row_title=cna_list_with_sig3,row_title_rot =0,row_title_gp = gpar(fontsize = 8),
#                  row_labels=rep(c("CNA","RNA","Protein"),length(cna_list_with_sig3)),
#                  # col = colors,
#                  show_column_names = F,height = 2,column_split =group_list_combo2,column_gap  = unit(1.5, "mm"),column_title = NULL,
#                  row_split  =rep(1:length(cna_list_with_sig3),each=3),  row_gap = unit(1.5, "mm") ,
#                  top_annotation=ha,
#                  show_row_names=T,na_col = "#F1F1F1",row_names_gp = gpar(fontsize = 8))
# draw(ht_mut)
# dev.off()

colors <- circlize::colorRamp2(seq(-2, 2, length=3), c("#3A53A3","#E7E6E6", "#EC1D3B"))

pdf("overall heatmap CNA-RNA-PRO order.pdf",width =20 ,height = 7.5)
ht_mut = Heatmap(as.matrix(Commo_data.frame),cluster_rows = F,cluster_columns = T,show_column_dend = F,column_dend_reorder=T,
                 border_gp = gpar(col = "Grey", lwd = 1.5),
                 rect_gp = gpar(col = 'white', lwd = 1),
                 row_title=cna_list_with_sig3,row_title_rot =0,row_title_gp = gpar(fontsize = 8),
                 row_labels=rep(c("CNA","RNA","Protein"),length(cna_list_with_sig3)),
                 col = colors,
                 show_column_names = F,height = 2,column_split =group_list_combo2,column_gap  = unit(1.5, "mm"),column_title = NULL,
                 row_split  =rep(1:length(cna_list_with_sig3),each=3),  row_gap = unit(1.5, "mm") ,
                 top_annotation=ha,
                 show_row_names=T,na_col = "#F1F1F1",row_names_gp = gpar(fontsize = 8))
draw(ht_mut)
dev.off()


