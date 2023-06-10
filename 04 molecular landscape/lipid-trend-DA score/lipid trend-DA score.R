load("metabolome_trend.Rdata")
load("KEGG_annotation.Rdata")

plot_meta <- c(rownames(trend_HR_negative_lipid[trend_HR_negative_lipid$Rho>=0.3 & trend_HR_negative_lipid$Rho_FDR<0.5,]),rownames(trend_HR_positive_lipid[trend_HR_positive_lipid$Rho>=0.3 &trend_HR_positive_lipid$Rho_FDR<0.5,]))
plot <- cbind(trend_HR_positive_lipid[plot_meta,1:3],trend_HR_negative_lipid[plot_meta,1:3]) %>% t() %>% scale() %>% t()

anno <- FUSCC_HER2_low_project_metabolite_lipid_anno[plot_meta,]
anno <- arrange(anno,`Lipid.super.class`,`Putative_metabolite_name`)
anno <- anno[,1:2] %>% as.matrix()
plot <- plot[rownames(anno),] %>% as.matrix()

rownames(anno) <- anno[,1]
rownames(plot) <- anno[,1]

pheatmap(plot,cluster_rows = F,cluster_cols = F,border_color = NA)

lipid_anno <- FUSCC_HER2_low_project_metabolite_lipid_anno[,c(4,7)]
lipid_anno$peak_name <- rownames(lipid_anno)
lipid_anno <- lipid_anno[,c(3:1)]
colnames(lipid_anno) <- c("peak_name","KEGG_ID","Lipid_abbreviation")

polar_anno <- FUSCC_HER2_low_project_metabolite_polar_anno[,c(1,4)]
polar_anno$Lipid_abbreviation <- ""
colnames(polar_anno) <- c("peak_name","KEGG_ID","Lipid_abbreviation")

peak_anno <- rbind(lipid_anno,polar_anno)

data <- rbind(trend_HR_negative_lipid,trend_HR_negative_polar)
# data <- rbind(trend_HR_positive_lipid,trend_HR_positive_polar)

Up_peak_name <- rownames(data[data$level1<data$level2 & data$level2<data$level3 & data$Rho>0.3,])
Down_peak_name <- rownames(data[data$level1>data$level2 & data$level2>data$level3 & data$Rho< -0.3,])

DA_score_matrix <- matrix(ncol=13,nrow=nrow(KEGG_annotation))
colnames(DA_score_matrix) <- c("Path_Name","Path_Class","Path_KEGG_Num","Anno_KEGG_Num","Anno_Lip_Num","Anno_Num",
                               "Up_KEGG_Num","Down_KEGG_Num","Up_Lip_Num","Down_Lip_Num","Up_Num","Down_Num",
                               "DA_Score")
rownames(DA_score_matrix) <- KEGG_annotation$Pathway
DA_score_matrix[,1:2] <- as.matrix(KEGG_annotation[,2:3])
DA_score_matrix <- as.data.frame(DA_score_matrix,stringsAsFactors=F)

for (i in 1:nrow(DA_score_matrix)){
  
  a <- strsplit(KEGG_annotation$Metabolites_1[i],";")[[1]]
  DA_score_matrix[i,"Path_KEGG_Num"] <- as.numeric(length(a))
  DA_score_matrix[i,"Anno_KEGG_Num"] <- as.numeric(length(rownames(peak_anno)[peak_anno$KEGG_ID%in%a]))
  
  b <- strsplit(KEGG_annotation$Metabolites_2[i],";")[[1]]
  DA_score_matrix[i,"Anno_Lip_Num"] <- as.numeric(length(rownames(peak_anno)[peak_anno$Lipid_abbreviation%in%b]))
  DA_score_matrix[i,"Anno_Num"] <- as.numeric(length(rownames(peak_anno)[peak_anno$KEGG_ID%in%a | peak_anno$Lipid_abbreviation%in%b]))
  
  DA_score_matrix[i,"Up_KEGG_Num"] <- as.numeric(length(rownames(peak_anno)[peak_anno$KEGG_ID%in%a & rownames(peak_anno)%in%Up_peak_name]))
  DA_score_matrix[i,"Down_KEGG_Num"] <- as.numeric(length(rownames(peak_anno)[peak_anno$KEGG_ID%in%a & rownames(peak_anno)%in%Down_peak_name]))
  DA_score_matrix[i,"Up_Lip_Num"] <- as.numeric(length(rownames(peak_anno)[is.na(peak_anno$Lipid_abbreviation)==F & peak_anno$Lipid_abbreviation%in%b & rownames(peak_anno)%in%Up_peak_name]))
  DA_score_matrix[i,"Down_Lip_Num"] <- as.numeric(length(rownames(peak_anno)[is.na(peak_anno$Lipid_abbreviation)==F & peak_anno$Lipid_abbreviation%in%b & rownames(peak_anno)%in%Down_peak_name]))
  DA_score_matrix[i,"Up_Num"] <- as.numeric(length(rownames(peak_anno)[rownames(peak_anno)%in%Up_peak_name & (peak_anno$KEGG_ID%in%a | (is.na(peak_anno$Lipid_abbreviation)==F & peak_anno$Lipid_abbreviation%in%b))]))
  DA_score_matrix[i,"Down_Num"] <- as.numeric(length(rownames(peak_anno)[rownames(peak_anno)%in%Down_peak_name & (peak_anno$KEGG_ID%in%a | (is.na(peak_anno$Lipid_abbreviation)==F & peak_anno$Lipid_abbreviation%in%b))]))
  if (as.numeric(DA_score_matrix[i,"Anno_Num"])<=3) {DA_score_matrix[i,"DA_Score"] <- NA} 
  else DA_score_matrix[i,"DA_Score"] <- (as.numeric(DA_score_matrix[i,"Up_Num"])-as.numeric(DA_score_matrix[i,"Down_Num"]))/as.numeric(DA_score_matrix[i,"Anno_Num"])
  
}


for(i in 3:ncol(DA_score_matrix)){
  DA_score_matrix[,i] <- as.numeric(DA_score_matrix[,i])
}

DA_score_matrix <- DA_score_matrix[is.na(DA_score_matrix$DA_Score)==F,]
DA_score_matrix$Pathway_order <- c(1:nrow(DA_score_matrix))
DA_score_matrix$Path_Name_factor <- factor(DA_score_matrix$Path_Name,levels=DA_score_matrix$Path_Name)

ggplot(data=DA_score_matrix,aes(x=DA_Score,y=Path_Name_factor))+
  geom_point(aes(size=log2(DA_score_matrix$Anno_Num),color=DA_score_matrix$Path_Class)) +
  geom_segment(aes(color=DA_score_matrix$Path_Class),x=rep(0,nrow(DA_score_matrix)),xend=DA_score_matrix$DA_Score,y=DA_score_matrix$Pathway_order,yend=DA_score_matrix$Pathway_order)+
  theme_light(base_size = 12, base_family = "") +
  scale_x_continuous(limits= c(-0.3,0.3),breaks = c(-0.3,0,0.3))+
  scale_color_manual(values=brewer.pal(11,"Paired"))


