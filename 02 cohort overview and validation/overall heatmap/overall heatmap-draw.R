
library(ComplexHeatmap)
library(circlize)


##############################################################
##mutation anno
##############################################################

##anno row#################################
clinical_anno$HR_status <- factor(clinical_anno$HR_status,levels = c("Positive","Negative"))
total_order <- rownames(clinical_anno)

mut_anno_matrix <- matrix(NA,nrow = 11 ,ncol = length(rownames(mut_matrix[,total_order])))
rownames(mut_anno_matrix) <- c(
  'Frame_Shift_Del',  'Missense_Mutation',  'Nonsense_Mutation',  'Multi_Hit',  'Frame_Shift_Ins',  'In_Frame_Ins',  'Splice_Site',
  'In_Frame_Del',  "WT",  "NA",  "Percentage"
)
ttt <- c(  'Frame_Shift_Del',  'Missense_Mutation',  'Nonsense_Mutation',  'Multi_Hit',  'Frame_Shift_Ins',
           'In_Frame_Ins',  'Splice_Site',  'In_Frame_Del',  "WT",  "NA")

##
colnames(mut_anno_matrix) <- (rownames(mut_matrix[,total_order]))
mut_anno_matrix <- data.frame(mut_anno_matrix)
for (i in ((rownames(mut_matrix[,total_order])))){
  for (j in ttt[1:9]){
    mut_anno_matrix[j,i] <- length(mut_matrix[i,][mut_matrix[i,] %in% j]) 
  }
  tttt2 <- as.numeric(na.omit(mut_anno_matrix[,i]))
  mut_anno_matrix[11,i] <- (round(sum(tttt2[1:8])/sum(tttt2[1:9]),2)*100)
}

##gene name with percentage
mut_matrix_per <- mut_matrix
rownames(mut_matrix_per)[rownames(mut_matrix_per) %in% sig_dif_names_hr_pos] <- paste(rownames(mut_matrix_per)[rownames(mut_matrix_per) %in% sig_dif_names_hr_pos],"*",sep = " ")
rownames(mut_matrix_per)[rownames(mut_matrix_per) %in% sig_dif_names_hr_neg] <- paste(rownames(mut_matrix_per)[rownames(mut_matrix_per) %in% sig_dif_names_hr_neg],"#",sep = " ")

rownames(mut_matrix_per) <- paste(mut_anno_matrix["Percentage",],"%"," ",rownames(mut_matrix_per),sep = "")



##anno col#################################
df = clinical_anno[total_order,c("HR_status","HER2_low_status_long_RE","PAM50","Menopause","DMFS_status","Lymph_nodes")]
colnames(df) <- c("HR status","HER2 IHC score","PAM50","Menopause","DMFS status","Lymph nodes")

df[df=="Her2"] <- "HER2"

ha = HeatmapAnnotation(df = df,
                       col  = list(`HR status` = c("Positive" = "#4d8ee4", "Negative" = "#e4574d"),
                                   `HER2 IHC score` = c("HER2_1" = "#7AB801","HER2_2" = "#F2AF01"),
                                   `Lymph nodes` = c("Positive" = "#BC102B", "Negative" = "#EAEAEA"),
                                   `PAM50`=c("LumA" = "#000286", "LumB" = "#87CCFA","HER2" = "#FF67B8", "Basal" = "#FE0302","Normal"="#36E200"),
                                   `Menopause`=c("Yes" = "#293241", "No" = "#98C1D9","Male" = "#EAEAEA"),
                                   `DMFS status`=c("1" = "#293241", "0" = "#CCC8DE")
                       ),
                       na_col = "#EAEAEA")
ha
draw(ha)

colors <- structure(c(c("#606e8a","#0C55A4","#7db2c1","#018a9c","#fe586f","#ff9601","#8d2b07","#015638","#EEEBCE"),"#EAEAEA","#EAEAEA"), names=c(  'Frame_Shift_Del',  'Missense_Mutation',  'Nonsense_Mutation',  'Multi_Hit',
                                                                                                                                                   'Frame_Shift_Ins',  'In_Frame_Ins',  'Splice_Site',  'In_Frame_Del',"Translation_Start_Site",  "WT",  "NA"))

#draw(ht_mut)
mut_anno_matrix_per <- mut_anno_matrix
colnames(mut_anno_matrix_per) <- rownames(mut_matrix_per)

ha2 = rowAnnotation(` ` = row_anno_barplot(data.frame(t((mut_anno_matrix_per[1:8,]))),gp = gpar(fill = colors ,col=NA), width = unit(2, "cm"),border = F,axis = F)    )


mut_row_split <- rep(2,length(gene_order))
mut_row_split[1:10] <- 1

  ht_mut = Heatmap(as.matrix(mut_matrix_per[,total_order]), name = "Variant classification",row_title_gp = gpar(fontsize = 10) ,row_title = "Somatic mutation",cluster_rows = F,cluster_columns = F,
                  col = colors,show_column_names = F,height = 2,column_split =total_order_cut_col,column_gap  = unit(1.5, "mm"),column_title = NULL,
                  row_split    =mut_row_split,  row_gap = unit(1, "mm") ,
                  top_annotation=ha,
                  show_row_names=T,right_annotation = ha2,na_col = "#F1F1F1",row_names_gp = gpar(fontsize = 8))

draw(ht_mut)


##############################################################
##CNV
##############################################################

##amp#############################################

#col
colors <- structure(c("#F41D1D","#F98181","#EAEAEA","#6B7BA4","#0A2A7C"),names=c("Amplification","Gain","Neutral","Loss","Deletion"))
colors <- structure(c("#BC102B","#DF989E","#EAEAEA","#AEC4D6","#5385AC"),names=c("Amplification","Gain","Neutral","Loss","Deletion"))
##sub
gene_amp_matrix_2 <- gene_amp_matrix
gene_amp_matrix_2[gene_amp_matrix_2 == 2] <- "Amplification"
gene_amp_matrix_2[gene_amp_matrix_2 == 1] <- "Gain"
gene_amp_matrix_2[gene_amp_matrix_2 == 0] <- "Neutral"
gene_amp_matrix_2[gene_amp_matrix_2 == -1] <- "Loss"
gene_amp_matrix_2[gene_amp_matrix_2 == -2] <- "Deletion"

amp_anno_matrix <- matrix(NA,nrow = 6 ,ncol = length(rownames(gene_amp_matrix_2[,total_order])))
rownames(amp_anno_matrix) <- c("Amplification","Gain","Neutral","Loss","Deletion","Percentage")
ttt <- c("Amplification","Gain","Neutral","Loss","Deletion","Percentage")

colnames(amp_anno_matrix) <- (rownames(gene_amp_matrix_2[,total_order]))
amp_anno_matrix <- data.frame(amp_anno_matrix)
for (i in ((rownames(gene_amp_matrix_2[,total_order])))){##样本
  for (j in ttt[1:5]){##种类
    amp_anno_matrix[j,i] <- length(gene_amp_matrix_2[i,][gene_amp_matrix_2[i,] %in% j]) 
  }
  tttt2 <- as.numeric(na.omit(amp_anno_matrix[,i]))
  amp_anno_matrix[6,i] <- (round(sum(tttt2[1:2])/sum(tttt2[1:5]),2)*100)
}

amp_anno_matrix_per <- amp_anno_matrix
colnames(amp_anno_matrix_per) <- paste(amp_anno_matrix_per["Percentage",],"%"," ",colnames(amp_anno_matrix_per),sep = "")

rownames(gene_amp_matrix_2) <- colnames(amp_anno_matrix_per)

##del#############################################


gene_del_matrix_2 <- gene_del_matrix
gene_del_matrix_2[gene_del_matrix_2 == 2] <- "Amplification"
gene_del_matrix_2[gene_del_matrix_2 == 1] <- "Gain"
gene_del_matrix_2[gene_del_matrix_2 == 0] <- "Neutral"
gene_del_matrix_2[gene_del_matrix_2 == -1] <- "Loss"
gene_del_matrix_2[gene_del_matrix_2 == -2] <- "Deletion"

del_anno_matrix <- matrix(NA,nrow = 6 ,ncol = length(rownames(gene_del_matrix_2[,total_order])))
rownames(del_anno_matrix) <- c("Amplification","Gain","Neutral","Loss","Deletion","Percentage")
ttt <- c("Amplification","Gain","Neutral","Loss","Deletion","Percentage")

colnames(del_anno_matrix) <- (rownames(gene_del_matrix_2[,total_order]))
del_anno_matrix <- data.frame(del_anno_matrix)
for (i in ((rownames(gene_del_matrix_2[,total_order])))){##样本
  for (j in ttt[1:5]){##种类
    del_anno_matrix[j,i] <- length(gene_del_matrix_2[i,][gene_del_matrix_2[i,] %in% j]) 
  }
  tttt2 <- as.numeric(na.omit(del_anno_matrix[,i]))
  del_anno_matrix[6,i] <- (round(sum(tttt2[4:5])/sum(tttt2[1:5]),2)*100)
}

del_anno_matrix_per <- del_anno_matrix
colnames(del_anno_matrix_per) <- paste(del_anno_matrix_per["Percentage",],"%"," ",colnames(del_anno_matrix_per),sep = "")

rownames(gene_del_matrix_2) <- colnames(del_anno_matrix_per)


###row anno###################################

amp_anno_clean <- amp_anno_matrix_per[1:5,]
amp_anno_clean[3:5,] <- 0


rownames(amp_anno_clean)[4:5] <- c("Deletion","Loss")


del_anno_clean <- del_anno_matrix_per[1:5,]
del_anno_clean[1:3,] <- 0
del_anno_clean[3,] <- del_anno_clean[4,] 
del_anno_clean[4,] <- del_anno_clean[5,] 
del_anno_clean[5,] <- del_anno_clean[3,] 
del_anno_clean[1:3,] <- 0
rownames(del_anno_clean)[4:5] <- c("Deletion","Loss")

ha4 = rowAnnotation( ` ` = row_anno_barplot(data.frame(t(cbind(amp_anno_clean,del_anno_clean))),gp = gpar(fill = structure(c("#BC102B","#DF989E","#EAEAEA","#5385AC","#AEC4D6"),names=c("Amplification","Gain","Neutral","Deletion","Loss")) ,col=NA), width = unit(2, "cm"),border = F,axis = F)    )

###draw together###################################

ht_CNV_gene=Heatmap(as.matrix(rbind(gene_amp_matrix_2,gene_del_matrix_2)[,total_order]), name = "CNA", 
                    row_title_gp = gpar(fontsize = 10) ,row_title = "CNA",cluster_rows = F,cluster_columns = F,col = colors,show_column_names = F,na_col = "#F1F1F1",
                    height = 1.5,column_split =total_order_cut_col,column_gap  = unit(1.5, "mm"),right_annotation = ha4,column_title = NULL,row_names_gp = gpar(fontsize = 8)) 
draw(ht_CNV_gene)


##############################################################
##RNA 
##############################################################

choose_matrix <- RNA_matrix_total
choose_matrix = t( scale( t( choose_matrix ) ) )

choose_matrix[choose_matrix > 2] = 2
choose_matrix[choose_matrix < -2] = -2

RNA_anno_matrix<- rowAnnotation(` ` = gene_anno_RNA[,"FDR"],
                                col  = list(` ` = c("<0.1" = "#E42417", ">=0.1" = "#EAEAEA")))

label_genes <- gene_anno_RNA$gene[gene_anno_RNA$gene %in%  c(read.csv("COSMIC CGC tier12.csv")[,1],paste0(read.csv("COSMIC CGC tier12.csv")[,1],"_dup"))]
label_genes <- c(label_genes,paste0(label_genes,"_dup"))

label_genes <- label_genes[gene_anno_RNA[label_genes,"FDR"] %in% c("<0.1")]

RNA_anno_matrix <- rowAnnotation(`FDR` = gene_anno_RNA[,"FDR"],
                                 col  = list(`FDR` = c("<0.1" = "#E42417", ">=0.1" = "#EAEAEA")),
                                 foo = anno_mark(at = which(rownames(choose_matrix) %in% label_genes),labels = label_genes, labels_gp  = gpar(fontsize = 8)),
                                 show_annotation_name = F)


library(circlize)
ht_RNA=Heatmap(as.matrix(choose_matrix[,total_order]), name = "Normalized RNA log2(FPKM+1)",  row_title_gp = gpar(fontsize = 10) ,row_title = "RNA",cluster_rows = F,cluster_columns = F,
               show_row_names = F,col = colorRamp2(seq(min(choose_matrix,na.rm = T), max(choose_matrix,na.rm = T), length=3), c("#009FFF","White", "#ec2F4B")),show_column_names = F,height = 1.2,na_col = "#F1F1F1",
               column_split =total_order_cut_col,column_gap  = unit(1.5, "mm"),
               row_split    =RNA_row_split,  row_gap = unit(0.5, "mm") ,column_title = NULL
               ,right_annotation = RNA_anno_matrix)
draw(ht_RNA)


##Sep_plot_RNA
ha_short = HeatmapAnnotation(df = df[total_order[! colSums(is.na(choose_matrix)) == nrow(choose_matrix)],c("HR status","HER2 IHC score","PAM50")],
                             col  = list(`HR status` = c("Positive" = "#4d8ee4", "Negative" = "#e4574d"),
                                         `HER2 IHC score` = c("HER2_1" = "#7AB801","HER2_2" = "#F2AF01"),
                                         `PAM50`=c("LumA" = "#000286", "LumB" = "#87CCFA","HER2" = "#FF67B8", "Basal" = "#FE0302","Normal"="#36E200"))
                             ,na_col = "#EAEAEA",simple_anno_size = unit(0.2, "cm"),annotation_name_gp = gpar(fontsize = 5))

pdf("Sep_plot_RNA.pdf",width = 12,height = 6)
Heatmap(as.matrix(choose_matrix[,total_order[! colSums(is.na(choose_matrix)) == nrow(choose_matrix)]]), name = "Normalized RNA log2(FPKM+1)",  row_title_gp = gpar(fontsize = 10) ,row_title = "RNA",cluster_rows = F,cluster_columns = F,
        show_row_names = T,col = colorRamp2(seq(min(choose_matrix,na.rm = T), max(choose_matrix,na.rm = T), length=3), c("#009FFF","White", "#ec2F4B")),show_column_names = F,height = 1.2,na_col = "#F1F1F1",
        column_split =total_order_cut_col[! colSums(is.na(choose_matrix)) == nrow(choose_matrix)],column_gap  = unit(1.5, "mm"),
        row_split    =RNA_row_split,  row_gap = unit(0.5, "mm") ,column_title = NULL
        ,right_annotation = RNA_anno_matrix,top_annotation = ha_short,row_names_gp = gpar(fontsize = 2),row_names_side="left")
dev.off()

{total_list <- paste(FUSCC_HER2_low_project_Cohort.Info$HR_status,FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_long_RE)
  total_list <- gsub("Positive ","HR-pos ",total_list)
  total_list <- gsub("Negative ","HR-neg ",total_list)
  total_list <- gsub("HER2_1","HER2 1+",total_list)
  total_list <- gsub("HER2_2","HER2 2+",total_list)
  names(total_list) <- rownames(FUSCC_HER2_low_project_Cohort.Info)
  table(total_list)
  
  RNA_matrix_total_sup <- cbind(RNA_matrix_total,Compare=rep(c("HR-pos HER2 1+ vs HER2 2+ (1+ upregualted)","HR-pos HER2 1+ vs HER2 2+ (2+ upregualted)","HR-neg HER2 1+ vs HER2 2+ (1+ upregualted)","HR-neg HER2 1+ vs HER2 2+ (2+ upregualted)"),each=nrow(RNA_matrix_total)/4),FDR=gene_anno_RNA$FDR)
  RNA_matrix_total_sup <- RNA_matrix_total_sup[,! colSums(is.na(RNA_matrix_total_sup)) == nrow(RNA_matrix_total_sup)]
  RNA_matrix_total_sup <- cbind(Gene_symbol=rownames(RNA_matrix_total_sup),RNA_matrix_total_sup)
  RNA_matrix_total_sup$Gene_symbol <- gsub("_dup","",RNA_matrix_total_sup$Gene_symbol)
  
  RNA_matrix_total_sup <- rbind(Group=c("Gene_symbol",total_list[colnames(RNA_matrix_total_sup)[2:(ncol(RNA_matrix_total_sup)-2)]],"Compare","FDR"),RNA_matrix_total_sup)
  RNA_matrix_total_sup$FDR <- gsub(">=","≥",RNA_matrix_total_sup$FDR)
  xlsx::write.xlsx(file = "RNA_matrix_total_sup.xlsx",RNA_matrix_total_sup,row.names = F)}

##PRO
##############################################################

choose_matrix <- PRO_matrix_total
choose_matrix = t( scale( t( choose_matrix ) ) )

choose_matrix[choose_matrix > 2] = 2
choose_matrix[choose_matrix < -2] = -2


PRO_anno_matrix<- rowAnnotation(`FDR` = gene_anno_PRO[,"FDR"],
                                col  = list(`FDR` = c("<0.1" = "#E42417", ">=0.1" = "#EAEAEA")))

label_genes <- gene_anno_PRO$gene[gene_anno_PRO$gene %in%  c(read.csv("COSMIC CGC tier12.csv")[,1],paste0(read.csv("COSMIC CGC tier12.csv")[,1],"_dup"))]


label_genes <- label_genes[gene_anno_PRO[label_genes,"FDR"] %in% c("<0.1")]

PRO_anno_matrix <- rowAnnotation(  
  `FDR` = gene_anno_PRO[,"FDR"],
  col  = list(`FDR` = c("<0.1" = "#E42417", ">=0.1" = "#EAEAEA")),
  foo = anno_mark(at = which(rownames(choose_matrix) %in% label_genes),labels = label_genes, labels_gp  = gpar(fontsize = 8)),
  show_annotation_name = F)


library(circlize)
ht_PRO=Heatmap(as.matrix(choose_matrix[,total_order]), name = "Normalized protein abundance",  row_title_gp = gpar(fontsize = 10) ,row_title = "Protein",
               cluster_rows = F,cluster_columns = F,show_row_names = F,col = colorRamp2(seq(min(choose_matrix,na.rm = T), max(choose_matrix,na.rm = T), length=3),c("#1D2671","White", "#C33764")),show_column_names = F,height = 1.2,na_col = "#F1F1F1",
               column_split =total_order_cut_col,column_gap  = unit(1.5, "mm"),
               row_split    =PRO_row_split,  row_gap = unit(0.5, "mm") ,column_title = NULL
               ,right_annotation = PRO_anno_matrix)
draw(ht_PRO)

##Sep_plot_PRO
ha_short = HeatmapAnnotation(df = df[total_order[! colSums(is.na(choose_matrix)) == nrow(choose_matrix)],c("HR status","HER2 IHC score","PAM50")],
                             col  = list(`HR status` = c("Positive" = "#4d8ee4", "Negative" = "#e4574d"),
                                         `HER2 IHC score` = c("HER2_1" = "#7AB801","HER2_2" = "#F2AF01"),
                                         `PAM50`=c("LumA" = "#000286", "LumB" = "#87CCFA","HER2" = "#FF67B8", "Basal" = "#FE0302","Normal"="#36E200"))
                             ,na_col = "#EAEAEA",simple_anno_size = unit(0.2, "cm"),annotation_name_gp = gpar(fontsize = 5))

pdf("Sep_plot_PRO.pdf",width = 12,height = 6)
Heatmap(as.matrix(choose_matrix[,total_order[! colSums(is.na(choose_matrix)) == nrow(choose_matrix)]]), name = "Normalized protein abundance",  row_title_gp = gpar(fontsize = 10) ,row_title = "Protein",cluster_rows = F,cluster_columns = F,
        show_row_names = T,col = colorRamp2(seq(min(choose_matrix,na.rm = T), max(choose_matrix,na.rm = T), length=3), c("#1D2671","White", "#C33764")),show_column_names = F,height = 1.2,na_col = "#F1F1F1",
        column_split =total_order_cut_col[! colSums(is.na(choose_matrix)) == nrow(choose_matrix)],column_gap  = unit(1.5, "mm"),
        row_split    =PRO_row_split,  row_gap = unit(0.5, "mm") ,column_title = NULL
        ,right_annotation = PRO_anno_matrix,top_annotation = ha_short,row_names_gp = gpar(fontsize = 2),row_names_side="left")
dev.off()

{total_list <- paste(FUSCC_HER2_low_project_Cohort.Info$HR_status,FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_long_RE)
  total_list <- gsub("Positive ","HR-pos ",total_list)
  total_list <- gsub("Negative ","HR-neg ",total_list)
  total_list <- gsub("HER2_1","HER2 1+",total_list)
  total_list <- gsub("HER2_2","HER2 2+",total_list)
  names(total_list) <- rownames(FUSCC_HER2_low_project_Cohort.Info)
  table(total_list)
  
  PRO_matrix_total_sup <- cbind(PRO_matrix_total,Compare=rep(c("HR-pos HER2 1+ vs HER2 2+ (1+ upregualted)","HR-pos HER2 1+ vs HER2 2+ (2+ upregualted)","HR-neg HER2 1+ vs HER2 2+ (1+ upregualted)","HR-neg HER2 1+ vs HER2 2+ (2+ upregualted)"),each=nrow(PRO_matrix_total)/4),FDR=gene_anno_PRO$FDR)
  PRO_matrix_total_sup <- PRO_matrix_total_sup[,! colSums(is.na(PRO_matrix_total_sup)) == nrow(PRO_matrix_total_sup)]
  PRO_matrix_total_sup <- cbind(Gene_symbol=rownames(PRO_matrix_total_sup),PRO_matrix_total_sup)
  PRO_matrix_total_sup$Gene_symbol <- gsub("_dup","",PRO_matrix_total_sup$Gene_symbol)
  
  PRO_matrix_total_sup <- rbind(Group=c("Gene_symbol",total_list[colnames(PRO_matrix_total_sup)[2:(ncol(PRO_matrix_total_sup)-2)]],"Compare","FDR"),PRO_matrix_total_sup)
  PRO_matrix_total_sup$FDR <- gsub(">=","≥",PRO_matrix_total_sup$FDR)
  xlsx::write.xlsx(file = "PRO_matrix_total_sup.xlsx",PRO_matrix_total_sup,row.names = F)}



##POL
##############################################################

choose_matrix <- POL_matrix_total
choose_matrix = t( scale( t( choose_matrix ) ) )

choose_matrix[choose_matrix > 2] = 2
choose_matrix[choose_matrix < -2] = -2


POL_anno_matrix<- rowAnnotation(`FDR` = gene_anno_POL[,"FDR"],
                                col  = list(`FDR` = c("<0.1" = "#E42417", ">=0.1" = "#EAEAEA")))

label_genes <- label_genes[gene_anno_POL[label_genes,"FDR"] %in% c("<0.1")]

POL_anno_matrix <- rowAnnotation(  
  `FDR` = gene_anno_POL[,"FDR"],
  col  = list(`FDR` = c("<0.1" = "#E42417", ">=0.1" = "#EAEAEA")),
  show_annotation_name = F)

library(circlize)

ht_POL=Heatmap(as.matrix(choose_matrix[,total_order]), name = "Normalized metabolite abundance",  row_title_gp = gpar(fontsize = 10) ,row_title = "     Polar \n  metabolites",
               cluster_rows = F,cluster_columns = F,show_row_names = F,col = colorRamp2(seq(min(choose_matrix,na.rm = T), max(choose_matrix,na.rm = T), length=3),c("#00A69C","White", "#d90057")),show_column_names = F,height = 0.6,na_col = "#F1F1F1",
               column_split =total_order_cut_col,column_gap  = unit(1.5, "mm"),
               row_split    =POL_row_split,  row_gap = unit(0.5, "mm") ,column_title = NULL
               ,right_annotation = POL_anno_matrix)
draw(ht_POL)

##Sep_plot_POL
ha_short = HeatmapAnnotation(df = df[total_order[! colSums(is.na(choose_matrix)) == nrow(choose_matrix)],c("HR status","HER2 IHC score","PAM50")],
                             col  = list(`HR status` = c("Positive" = "#4d8ee4", "Negative" = "#e4574d"),
                                         `HER2 IHC score` = c("HER2_1" = "#7AB801","HER2_2" = "#F2AF01"),
                                         `PAM50`=c("LumA" = "#000286", "LumB" = "#87CCFA","HER2" = "#FF67B8", "Basal" = "#FE0302","Normal"="#36E200"))
                             ,na_col = "#EAEAEA",simple_anno_size = unit(0.2, "cm"),annotation_name_gp = gpar(fontsize = 5))

pdf("Sep_plot_POL.pdf",width = 12,height = 3)
Heatmap(as.matrix(choose_matrix[,total_order[! colSums(is.na(choose_matrix)) == nrow(choose_matrix)]]), name = "Normalized metabolite abundance",  row_title_gp = gpar(fontsize = 10) ,row_title = "     Polar \n  metabolites",cluster_rows = F,cluster_columns = F,
        show_row_names = T,col = colorRamp2(seq(min(choose_matrix,na.rm = T), max(choose_matrix,na.rm = T), length=3), c("#00A69C","White", "#d90057")),show_column_names = F,height = 1.2,na_col = "#F1F1F1",
        column_split =total_order_cut_col[! colSums(is.na(choose_matrix)) == nrow(choose_matrix)],column_gap  = unit(1.5, "mm"),
        row_split    =POL_row_split,  row_gap = unit(0.5, "mm") ,column_title = NULL
        ,right_annotation = POL_anno_matrix,top_annotation = ha_short,row_names_gp = gpar(fontsize = 2),row_names_side="left")
dev.off()

{total_list <- paste(FUSCC_HER2_low_project_Cohort.Info$HR_status,FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_long_RE)
  total_list <- gsub("Positive ","HR-pos ",total_list)
  total_list <- gsub("Negative ","HR-neg ",total_list)
  total_list <- gsub("HER2_1","HER2 1+",total_list)
  total_list <- gsub("HER2_2","HER2 2+",total_list)
  names(total_list) <- rownames(FUSCC_HER2_low_project_Cohort.Info)
  table(total_list)
  
  POL_matrix_total_sup <- cbind(POL_matrix_total,Compare=rep(c("HR-pos HER2 1+ vs HER2 2+ (1+ upregualted)","HR-pos HER2 1+ vs HER2 2+ (2+ upregualted)","HR-neg HER2 1+ vs HER2 2+ (1+ upregualted)","HR-neg HER2 1+ vs HER2 2+ (2+ upregualted)"),each=nrow(POL_matrix_total)/4),FDR=gene_anno_POL$FDR)
  POL_matrix_total_sup <- POL_matrix_total_sup[,! colSums(is.na(POL_matrix_total_sup)) == nrow(POL_matrix_total_sup)]
  POL_matrix_total_sup <- cbind(Peak=rownames(POL_matrix_total_sup),POL_matrix_total_sup)
  POL_matrix_total_sup$Peak <- gsub("_dup","",POL_matrix_total_sup$Peak)
  POL_matrix_total_sup[,colnames(FUSCC_HER2_low_project_metabolite_polar_anno)] <- FUSCC_HER2_low_project_metabolite_polar_anno[POL_matrix_total_sup$Peak,]
  
  
  POL_matrix_total_sup <- rbind(Group=c("Peak",total_list[colnames(POL_matrix_total_sup)[2:(ncol(POL_matrix_total_sup)-2-ncol(FUSCC_HER2_low_project_metabolite_polar_anno))]],"Compare","FDR",colnames(FUSCC_HER2_low_project_metabolite_polar_anno)),POL_matrix_total_sup)
  POL_matrix_total_sup$FDR <- gsub(">=","≥",POL_matrix_total_sup$FDR)
  xlsx::write.xlsx(file = "POL_matrix_total_sup.xlsx",POL_matrix_total_sup,row.names = F)}


##LIP
##############################################################

choose_matrix <- LIP_matrix_total
choose_matrix = t( scale( t( choose_matrix ) ) )

choose_matrix[choose_matrix > 2] = 2
choose_matrix[choose_matrix < -2] = -2


LIP_anno_matrix<- rowAnnotation(`FDR` = gene_anno_LIP[,"FDR"],
                                col  = list(`FDR` = c("<0.1" = "#E42417", ">=0.1" = "#EAEAEA")))

label_genes <- label_genes[gene_anno_LIP[label_genes,"FDR"] %in% c("<0.1")]

LIP_anno_matrix <- rowAnnotation(  
  `FDR` = gene_anno_LIP[,"FDR"],
  col  = list(`FDR` = c("<0.1" = "#E42417", ">=0.1" = "#EAEAEA")),
  show_annotation_name = F)


library(circlize)
ht_LIP=Heatmap(as.matrix(choose_matrix[,total_order]), name = "Normalized metabolite abundance",  row_title_gp = gpar(fontsize = 10) ,row_title = "Lipids",
               cluster_rows = F,cluster_columns = F,show_row_names = F,col = colorRamp2(seq(min(choose_matrix,na.rm = T), max(choose_matrix,na.rm = T), length=3),c("#00A69C","White", "#d90057")),show_column_names = F,height = 0.6,na_col = "#F1F1F1",
               column_split =total_order_cut_col,column_gap  = unit(1.5, "mm"),
               row_split    =LIP_row_split,  row_gap = unit(0.5, "mm") ,column_title = NULL
               ,right_annotation = LIP_anno_matrix)
draw(ht_LIP)

##Sep_plot_LIP
ha_short = HeatmapAnnotation(df = df[total_order[! colSums(is.na(choose_matrix)) == nrow(choose_matrix)],c("HR status","HER2 IHC score","PAM50")],
                             col  = list(`HR status` = c("Positive" = "#4d8ee4", "Negative" = "#e4574d"),
                                         `HER2 IHC score` = c("HER2_1" = "#7AB801","HER2_2" = "#F2AF01"),
                                         `PAM50`=c("LumA" = "#000286", "LumB" = "#87CCFA","HER2" = "#FF67B8", "Basal" = "#FE0302","Normal"="#36E200"))
                             ,na_col = "#EAEAEA",simple_anno_size = unit(0.2, "cm"),annotation_name_gp = gpar(fontsize = 5))

pdf("Sep_plot_LIP.pdf",width = 12,height = 3)
Heatmap(as.matrix(choose_matrix[,total_order[! colSums(is.na(choose_matrix)) == nrow(choose_matrix)]]), name = "Normalized metabolite abundance",  row_title_gp = gpar(fontsize = 10) ,row_title = "Lipids",cluster_rows = F,cluster_columns = F,
        show_row_names = T,col = colorRamp2(seq(min(choose_matrix,na.rm = T), max(choose_matrix,na.rm = T), length=3), c("#00A69C","White", "#d90057")),show_column_names = F,height = 1.2,na_col = "#F1F1F1",
        column_split =total_order_cut_col[! colSums(is.na(choose_matrix)) == nrow(choose_matrix)],column_gap  = unit(1.5, "mm"),
        row_split    =LIP_row_split,  row_gap = unit(0.5, "mm") ,column_title = NULL
        ,right_annotation = LIP_anno_matrix,top_annotation = ha_short,row_names_gp = gpar(fontsize = 2),row_names_side="left")
dev.off()

{total_list <- paste(FUSCC_HER2_low_project_Cohort.Info$HR_status,FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_long_RE)
  total_list <- gsub("Positive ","HR-pos ",total_list)
  total_list <- gsub("Negative ","HR-neg ",total_list)
  total_list <- gsub("HER2_1","HER2 1+",total_list)
  total_list <- gsub("HER2_2","HER2 2+",total_list)
  names(total_list) <- rownames(FUSCC_HER2_low_project_Cohort.Info)
  table(total_list)
  
  LIP_matrix_total_sup <- cbind(LIP_matrix_total,Compare=rep(c("HR-pos HER2 1+ vs HER2 2+ (1+ upregualted)","HR-pos HER2 1+ vs HER2 2+ (2+ upregualted)","HR-neg HER2 1+ vs HER2 2+ (1+ upregualted)","HR-neg HER2 1+ vs HER2 2+ (2+ upregualted)"),each=nrow(LIP_matrix_total)/4),FDR=gene_anno_LIP$FDR)
  LIP_matrix_total_sup <- LIP_matrix_total_sup[,! colSums(is.na(LIP_matrix_total_sup)) == nrow(LIP_matrix_total_sup)]
  LIP_matrix_total_sup <- cbind(Peak=rownames(LIP_matrix_total_sup),LIP_matrix_total_sup)
  LIP_matrix_total_sup$Peak <- gsub("_dup","",LIP_matrix_total_sup$Peak)
  LIP_matrix_total_sup[,colnames(FUSCC_HER2_low_project_metabolite_lipid_anno)] <- FUSCC_HER2_low_project_metabolite_lipid_anno[LIP_matrix_total_sup$Peak,]
  
  
  LIP_matrix_total_sup <- rbind(Group=c("Peak",total_list[colnames(LIP_matrix_total_sup)[2:(ncol(LIP_matrix_total_sup)-2-ncol(FUSCC_HER2_low_project_metabolite_lipid_anno))]],"Compare","FDR",colnames(FUSCC_HER2_low_project_metabolite_lipid_anno)),LIP_matrix_total_sup)
  LIP_matrix_total_sup$FDR <- gsub(">=","≥",LIP_matrix_total_sup$FDR)
  xlsx::write.xlsx(file = "LIP_matrix_total_sup.xlsx",LIP_matrix_total_sup,row.names = F)}



##############################################################
##final output
##############################################################

pdf(file="Merge_ALL_HER2_low_4_SUBTYPE_temp.pdf",width = 12,height = 8)
ht_list = ht_mut %v% ht_CNV_gene %v% ht_RNA %v% ht_PRO %v% ht_POL %v% ht_LIP
draw(ht_list)
dev.off()
