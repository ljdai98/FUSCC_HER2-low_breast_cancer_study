load("proteome_trend.Rdata")
load("KEGG_annotation.Rdata")

library(pheatmap)

gene <- KEGG_annotation_symbol[KEGG_annotation_symbol$Sub_Catagory=="Lipid metabolism",2]

trend_HR_negative_pro <- trend_HR_negative_pro[rownames(trend_HR_negative_pro)%in%gene,]
trend_HR_positive_pro <- trend_HR_positive_pro[rownames(trend_HR_positive_pro)%in%gene,]

plot_gene <- c(rownames(trend_HR_negative_pro[trend_HR_negative_pro$Rho>=0.3 & trend_HR_negative_pro$Rho_FDR<0.5,]),rownames(trend_HR_positive_pro[trend_HR_positive_pro$Rho>=0.3 & trend_HR_positive_pro$Rho_FDR<0.5,]))

plot <- cbind(trend_HR_positive_pro[plot_gene,1:3],trend_HR_negative_pro[plot_gene,1:3]) %>% t() %>% scale() %>% t()

anno <- KEGG_annotation_symbol[KEGG_annotation_symbol$Symbol%in%plot_gene &KEGG_annotation_symbol$Sub_Catagory=="Lipid metabolism",]
anno <- anno[!duplicated(anno$Symbol),]
rownames(anno) <- anno$Symbol
anno <- anno[,-2]

pheatmap(plot,cluster_rows = F,cluster_cols = F,border_color = NA,annotation_row = anno)



