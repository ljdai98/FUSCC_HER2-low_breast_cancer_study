table(FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE)
# HER2_0      HER2_low HER2_positive 
# 91           434           182 

load("../GSVA/RNA_C2.Rdata")
esmicro <- data.frame(esmicro)
esmicro <- esmicro[substr(rownames(esmicro),1,8) %in% "REACTOME",]

##HR+ 0 1 2
list_a_samplenames<-rownames(FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_long_RE %in% c("HER2_0") & 
                                                                  FUSCC_HER2_low_project_Cohort.Info$HR_status %in% c("Positive") ,])
list_a_typename<-c("HR_pos_HER2_0")

list_b_samplenames<-rownames(FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_long_RE %in% c("HER2_1") & 
                                                                  FUSCC_HER2_low_project_Cohort.Info$HR_status %in% c("Positive")  ,])
list_b_typename<-c("HR_pos_HER2_1")

list_c_samplenames<-rownames(FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_long_RE %in% c("HER2_2") & 
                                                                  FUSCC_HER2_low_project_Cohort.Info$HR_status %in% c("Positive")  ,])
list_c_typename<-c("HR_pos_HER2_2")

group_list <- c(rep(list_a_typename,length(list_a_samplenames)),
                rep(list_b_typename,length(list_b_samplenames)),
                rep(list_c_typename,length(list_c_samplenames)))

names(group_list) <- c(list_a_samplenames,list_b_samplenames,list_c_samplenames)
table(group_list)
# group_list
# HR_pos_HER2_0     HR_pos_HER2_1 HR_pos_HER2_2_low 
# 63               179               182

RNA_group_list <- group_list
names(RNA_group_list) <- paste(names(RNA_group_list),"_RNA_T",sep = "")

RNA_group_list <- RNA_group_list[names(RNA_group_list) %in% colnames(esmicro)]


RNA_group_list_num <- RNA_group_list
RNA_group_list_num[RNA_group_list_num == "HR_pos_HER2_0"] <- 1
RNA_group_list_num[RNA_group_list_num == "HR_pos_HER2_1"] <- 2
RNA_group_list_num[RNA_group_list_num == "HR_pos_HER2_2"] <- 3
table(RNA_group_list_num)


esmicro_median_centered <- t(scale(t(esmicro)))

choose_gene <- rownames(esmicro_median_centered)

temp_matrix <- matrix(NA,8,nrow = length(choose_gene))
temp_matrix <- data.frame(temp_matrix)
rownames(temp_matrix) <-choose_gene
colnames(temp_matrix) <- c("P","level1","level2","level3","FDR","Rho","Rho_P","Rho_FDR")


library(ggplot2)
library(ggpubr)

for (i in rownames(temp_matrix)){
  temp_for_test <- data.frame(Expr=as.numeric(esmicro_median_centered[i,names(RNA_group_list)]),Cluster=RNA_group_list)
  if (all(is.nan(temp_for_test$Expr))){
    next
  }
  fit<-kruskal.test(Expr~Cluster,data = temp_for_test)
  temp_matrix[i,"P"]<- fit[["p.value"]]
  temp_matrix[i,2:4] <- aggregate(as.numeric(esmicro_median_centered[i,names(RNA_group_list)]),by = list(RNA_group_list),FUN= median)[,2]
  temp_matrix[i,"Rho"]<- cor(as.numeric(esmicro_median_centered[i,names(RNA_group_list)]),as.numeric(RNA_group_list_num),method = "spearman") ;
  temp_matrix[i,"Rho_P"]<- cor.test(as.numeric(esmicro_median_centered[i,names(RNA_group_list)]),as.numeric(RNA_group_list_num),method = "spearman")[["p.value"]]
}
temp_matrix$FDR <- p.adjust(temp_matrix$P,method = "fdr")
temp_matrix$Rho_FDR <- p.adjust(temp_matrix$Rho_P,method = "fdr")


write.csv(temp_matrix,file = "mRNA_c2_with_0_1_2_low(HR+)(REACTOME).csv")
temp_matrix <- read.csv(file = "mRNA_c2_with_0_1_2_low(HR+)(REACTOME).csv",header = T,row.names = 1)

temp_matrix <- na.omit(temp_matrix)


sig_genes_up <- rownames(temp_matrix[ temp_matrix$Rho > 0.2 & temp_matrix$Rho_FDR < 0.25,])
sig_genes_up <- na.omit(sig_genes_up)
length(sig_genes_up)
# [1] 2

sig_genes_down <- rownames(temp_matrix[ temp_matrix$Rho < -0.2 & temp_matrix$Rho_FDR < 0.25,])
sig_genes_down <- na.omit(sig_genes_down)
length(sig_genes_down)
# [1] 11

sig_genes_up <- sig_genes_up[order(temp_matrix[sig_genes_up,"Rho_FDR"],decreasing = F)]
sig_genes_down <- sig_genes_down[order(temp_matrix[sig_genes_down,"Rho_FDR"],decreasing = F)]

sig_genes_up_matrx <- temp_matrix[sig_genes_up,]
sig_genes_down_matrx <- temp_matrix[sig_genes_down,]

write.csv(file = "mRNA_c2_with_0_1_2_low(HR+)sig_genes_up(REACTOME).csv",sig_genes_up_matrx,quote = F)
write.csv(file = "mRNA_c2_with_0_1_2_low(HR+)sig_genes_down(REACTOME).csv",sig_genes_down_matrx,quote = F)

##increasing
plot_matrix <- temp_matrix[sig_genes_up,2:4]
plot_matrix_2 <- data.frame(Expr=c(plot_matrix[,1],plot_matrix[,2],plot_matrix[,3]),
                            Cluster = c(rep(colnames(plot_matrix)[1],each = (nrow(plot_matrix))),rep(colnames(plot_matrix)[2],each = (nrow(plot_matrix))),rep(colnames(plot_matrix)[3],each = (nrow(plot_matrix)))),
                            Gene = c(rownames(plot_matrix),rownames(plot_matrix),rownames(plot_matrix)))
plot_matrix_2 <- na.omit(plot_matrix_2)

library(stringr)
plot_matrix_2$iffat <- "No"
plot_matrix_2[str_detect(plot_matrix_2$Gene, "FATTY_ACID|BETA_OXIDATION|LIPID|STEROIDS|LDL_|BILE_ACID|GLYCOSPHINGOLIPID"),"iffat"]<- "Yes"
table(plot_matrix_2$iffat)/3
# No 
# 2
library(ggplot2)


pdf("Trend_mRNA_c2_with_0_1_2_low(HR+)_up(REACTOME).pdf",width = 2,height = 4)
ggplot(plot_matrix_2, aes(x = Cluster, y = Expr)) +
  geom_line(data=plot_matrix_2[plot_matrix_2$iffat %in%  "No",],aes(group = Gene), color = '#D3D3D3', lwd = .5)+
  geom_line(data=plot_matrix_2[plot_matrix_2$iffat %in%  "Yes",],aes(group = Gene), color = 'Red', lwd = .5)+
  geom_boxplot( show.legend = FALSE, width = 0.6,outlier.size = 0) +  
  ylim(c(-2.1,2.1))+
  theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black'), 
        panel.background = element_blank(), 
        axis.text = element_text(size = 9, color = 'black'), axis.title = element_text(size = 10, color = 'black')) 
dev.off()


##decreasing
plot_matrix <- temp_matrix[sig_genes_down,2:4]
plot_matrix_2 <- data.frame(Expr=c(plot_matrix[,1],plot_matrix[,2],plot_matrix[,3]),
                            Cluster = c(rep(colnames(plot_matrix)[1],each = (nrow(plot_matrix))),rep(colnames(plot_matrix)[2],each = (nrow(plot_matrix))),rep(colnames(plot_matrix)[3],each = (nrow(plot_matrix)))),
                            Gene = c(rownames(plot_matrix),rownames(plot_matrix),rownames(plot_matrix)))
plot_matrix_2 <- na.omit(plot_matrix_2)

library(stringr)
plot_matrix_2$iffat <- "No"
plot_matrix_2[str_detect(plot_matrix_2$Gene, "FATTY_ACID|BETA_OXIDATION|LIPID|STEROIDS|LDL_|BILE_ACID|GLYCOSPHINGOLIPID"),"iffat"]<- "Yes"
table(plot_matrix_2$iffat)/3
# No 
# 11
library(ggplot2)

pdf("Trend_mRNA_c2_with_0_1_2_low(HR+)_down(REACTOME).pdf",width = 2,height = 4)
ggplot(plot_matrix_2, aes(x = Cluster, y = Expr)) +
  geom_line(data=plot_matrix_2[plot_matrix_2$iffat %in%  "No",],aes(group = Gene), color = '#D3D3D3', lwd = .5)+
  geom_line(data=plot_matrix_2[plot_matrix_2$iffat %in%  "Yes",],aes(group = Gene), color = 'Blue', lwd = .5)+
  geom_boxplot(show.legend = FALSE, width = 0.6,outlier.size = 0) +  
  ylim(c(-2.1,2.1))+
  theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black'), 
        panel.background = element_blank(), 
        axis.text = element_text(size = 9, color = 'black'), axis.title = element_text(size = 10, color = 'black')) 
dev.off()


go_up_and_down <- rbind(cbind(temp_matrix[substr(rownames(temp_matrix),1,8) %in% "REACTOME",],Category="Not significant"))
go_up_and_down$Description <-  rownames(go_up_and_down)
go_up_and_down[sig_genes_up,"Category"] <- "Increasing"
go_up_and_down[sig_genes_down,"Category"] <- "Decreasing"
go_up_and_down$Category <- as.factor(go_up_and_down$Category)
go_up_and_down$iffat <- "No"
go_up_and_down[sig_genes_up[str_detect(sig_genes_up, "FATTY_ACID|BETA_OXIDATION|LIPID|STEROIDS|LDL_|BILE_ACID|GLYCOSPHINGOLIPID")],"iffat"]<- "Yes"

write.csv(file = "Trend_mRNA_c2_with_0_1_2_low(HR+).csv",go_up_and_down,quote = F)

library(ggpubr)
library(ggrepel)

pdf("Trend_mRNA_c2_with_0_1_2_low(HR+)(vol)(labled)(with fat anno).pdf",height = 6,width =6)
ggplot(data=go_up_and_down, aes(x = Rho,y = -log10(Rho_P), colour=iffat)) +
  geom_point(alpha=0.8, size=3.5,col="Grey")+
  geom_point(data=go_up_and_down[go_up_and_down$iffat %in% "Yes",],alpha=0.4, size=3.5,col="#D6604D")+
  labs(x="Spearman Rho",
       y="-log10 (adjusted P)")+
  theme_bw()+geom_text_repel(data=go_up_and_down[go_up_and_down$iffat %in% "Yes",],
                             aes(label=Description),col="black",alpha = 0.8,size = 2)+xlim(-0.5,0.5)+ylim(0,6)
dev.off()


##HR- 0 1 2
list_a_samplenames<-rownames(FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_long_RE %in% c("HER2_0") & 
                                                                  FUSCC_HER2_low_project_Cohort.Info$HR_status %in% c("Negative") ,])
list_a_typename<-c("HR_neg_HER2_0")

list_b_samplenames<-rownames(FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_long_RE %in% c("HER2_1") & 
                                                                  FUSCC_HER2_low_project_Cohort.Info$HR_status %in% c("Negative")  ,])
list_b_typename<-c("HR_neg_HER2_1")

list_c_samplenames<-rownames(FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_long_RE %in% c("HER2_2") & 
                                                                  FUSCC_HER2_low_project_Cohort.Info$HR_status %in% c("Negative")  ,])
list_c_typename<-c("HR_neg_HER2_2")

group_list <- c(rep(list_a_typename,length(list_a_samplenames)),
                rep(list_b_typename,length(list_b_samplenames)),
                rep(list_c_typename,length(list_c_samplenames)))

names(group_list) <- c(list_a_samplenames,list_b_samplenames,list_c_samplenames)
table(group_list)
# HR_neg_HER2_0 HR_neg_HER2_1 HR_neg_HER2_2 
# 28            48            25

RNA_group_list <- group_list
names(RNA_group_list) <- paste(names(RNA_group_list),"_RNA_T",sep = "")

RNA_group_list <- RNA_group_list[names(RNA_group_list) %in% colnames(esmicro)]


RNA_group_list_num <- RNA_group_list
RNA_group_list_num[RNA_group_list_num == "HR_neg_HER2_0"] <- 1
RNA_group_list_num[RNA_group_list_num == "HR_neg_HER2_1"] <- 2
RNA_group_list_num[RNA_group_list_num == "HR_neg_HER2_2"] <- 3
table(RNA_group_list_num)


esmicro_median_centered <- t(scale(t(esmicro)))

choose_gene <- rownames(esmicro_median_centered)

temp_matrix <- matrix(NA,8,nrow = length(choose_gene))
temp_matrix <- data.frame(temp_matrix)
rownames(temp_matrix) <-choose_gene
colnames(temp_matrix) <- c("P","level1","level2","level3","FDR","Rho","Rho_P","Rho_FDR")


library(ggplot2)
library(ggpubr)

for (i in rownames(temp_matrix)){
  temp_for_test <- data.frame(Expr=as.numeric(esmicro_median_centered[i,names(RNA_group_list)]),Cluster=RNA_group_list)
  if (all(is.nan(temp_for_test$Expr))){
    next
  }
  fit<-kruskal.test(Expr~Cluster,data = temp_for_test)
  temp_matrix[i,"P"]<- fit[["p.value"]]
  temp_matrix[i,2:4] <- aggregate(as.numeric(esmicro_median_centered[i,names(RNA_group_list)]),by = list(RNA_group_list),FUN= median)[,2]
  temp_matrix[i,"Rho"]<- cor(as.numeric(esmicro_median_centered[i,names(RNA_group_list)]),as.numeric(RNA_group_list_num),method = "spearman") ;
  temp_matrix[i,"Rho_P"]<- cor.test(as.numeric(esmicro_median_centered[i,names(RNA_group_list)]),as.numeric(RNA_group_list_num),method = "spearman")[["p.value"]]
}
temp_matrix$FDR <- p.adjust(temp_matrix$P,method = "fdr")
temp_matrix$Rho_FDR <- p.adjust(temp_matrix$Rho_P,method = "fdr")


write.csv(temp_matrix,file = "mRNA_c2_with_0_1_2_low(HR-)(REACTOME).csv")
temp_matrix <- read.csv(file = "mRNA_c2_with_0_1_2_low(HR-)(REACTOME).csv",header = T,row.names = 1)

temp_matrix <- na.omit(temp_matrix)


sig_genes_up <- rownames(temp_matrix[ temp_matrix$Rho > 0.2 & temp_matrix$Rho_FDR < 0.25,])
sig_genes_up <- na.omit(sig_genes_up)
length(sig_genes_up)
# [1] 161

sig_genes_down <- rownames(temp_matrix[ temp_matrix$Rho < -0.2 & temp_matrix$Rho_FDR < 0.25,])
sig_genes_down <- na.omit(sig_genes_down)
length(sig_genes_down)
# [1] 285

sig_genes_up <- sig_genes_up[order(temp_matrix[sig_genes_up,"Rho_FDR"],decreasing = F)]
sig_genes_down <- sig_genes_down[order(temp_matrix[sig_genes_down,"Rho_FDR"],decreasing = F)]

sig_genes_up_matrx <- temp_matrix[sig_genes_up,]
sig_genes_down_matrx <- temp_matrix[sig_genes_down,]

write.csv(file = "mRNA_c2_with_0_1_2_low(HR-)sig_genes_up(REACTOME).csv",sig_genes_up_matrx,quote = F)
write.csv(file = "mRNA_c2_with_0_1_2_low(HR-)sig_genes_down(REACTOME).csv",sig_genes_down_matrx,quote = F)

##increasing
plot_matrix <- temp_matrix[sig_genes_up,2:4]
plot_matrix_2 <- data.frame(Expr=c(plot_matrix[,1],plot_matrix[,2],plot_matrix[,3]),
                            Cluster = c(rep(colnames(plot_matrix)[1],each = (nrow(plot_matrix))),rep(colnames(plot_matrix)[2],each = (nrow(plot_matrix))),rep(colnames(plot_matrix)[3],each = (nrow(plot_matrix)))),
                            Gene = c(rownames(plot_matrix),rownames(plot_matrix),rownames(plot_matrix)))
plot_matrix_2 <- na.omit(plot_matrix_2)

library(stringr)
plot_matrix_2$iffat <- "No"
plot_matrix_2[str_detect(plot_matrix_2$Gene, "FATTY_ACID|BETA_OXIDATION|LIPID|STEROIDS|LDL_|BILE_ACID|GLYCOSPHINGOLIPID"),"iffat"]<- "Yes"
table(plot_matrix_2$iffat)/3
# No Yes 
# 142  19

library(ggplot2)


pdf("Trend_mRNA_c2_with_0_1_2_low(HR-)_up(REACTOME).pdf",width = 2,height = 4)
ggplot(plot_matrix_2, aes(x = Cluster, y = Expr)) +
  geom_line(data=plot_matrix_2[plot_matrix_2$iffat %in%  "No",],aes(group = Gene), color = '#D3D3D3', lwd = .5)+
  geom_line(data=plot_matrix_2[plot_matrix_2$iffat %in%  "Yes",],aes(group = Gene), color = 'Red', lwd = .5)+
  geom_boxplot( show.legend = FALSE, width = 0.6,outlier.size = 0) +  
  ylim(c(-2.1,2.1))+
  theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black'), 
        panel.background = element_blank(), 
        axis.text = element_text(size = 9, color = 'black'), axis.title = element_text(size = 10, color = 'black')) 
dev.off()


##decreasing
plot_matrix <- temp_matrix[sig_genes_down,2:4]
plot_matrix_2 <- data.frame(Expr=c(plot_matrix[,1],plot_matrix[,2],plot_matrix[,3]),
                            Cluster = c(rep(colnames(plot_matrix)[1],each = (nrow(plot_matrix))),rep(colnames(plot_matrix)[2],each = (nrow(plot_matrix))),rep(colnames(plot_matrix)[3],each = (nrow(plot_matrix)))),
                            Gene = c(rownames(plot_matrix),rownames(plot_matrix),rownames(plot_matrix)))
plot_matrix_2 <- na.omit(plot_matrix_2)

library(stringr)
plot_matrix_2$iffat <- "No"
plot_matrix_2[str_detect(plot_matrix_2$Gene, "FATTY_ACID|BETA_OXIDATION|LIPID|STEROIDS|LDL_|BILE_ACID|GLYCOSPHINGOLIPID"),"iffat"]<- "Yes"
table(plot_matrix_2$iffat)/3
# No 
# 285 
library(ggplot2)

pdf("Trend_mRNA_c2_with_0_1_2_low(HR-)_down(REACTOME).pdf",width = 2,height = 4)
ggplot(plot_matrix_2, aes(x = Cluster, y = Expr)) +
  geom_line(data=plot_matrix_2[plot_matrix_2$iffat %in%  "No",],aes(group = Gene), color = '#D3D3D3', lwd = .5)+
  geom_line(data=plot_matrix_2[plot_matrix_2$iffat %in%  "Yes",],aes(group = Gene), color = 'Blue', lwd = .5)+
  geom_boxplot(show.legend = FALSE, width = 0.6,outlier.size = 0) +  
  ylim(c(-2.1,2.1))+
  theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black'), 
        panel.background = element_blank(), 
        axis.text = element_text(size = 9, color = 'black'), axis.title = element_text(size = 10, color = 'black')) 
dev.off()


go_up_and_down <- rbind(cbind(temp_matrix[substr(rownames(temp_matrix),1,8) %in% "REACTOME",],Category="Not significant"))
go_up_and_down$Description <-  rownames(go_up_and_down)
go_up_and_down[sig_genes_up,"Category"] <- "Increasing"
go_up_and_down[sig_genes_down,"Category"] <- "Decreasing"
go_up_and_down$Category <- as.factor(go_up_and_down$Category)
go_up_and_down$iffat <- "No"
go_up_and_down[sig_genes_up[str_detect(sig_genes_up, "FATTY_ACID|BETA_OXIDATION|LIPID|STEROIDS|LDL_|BILE_ACID|GLYCOSPHINGOLIPID")],"iffat"]<- "Yes"

write.csv(file = "Trend_mRNA_c2_with_0_1_2_low(HR-).csv",go_up_and_down,quote = F)

library(ggpubr)
library(ggrepel)

pdf("Trend_mRNA_c2_with_0_1_2_low(HR-)(vol)(labled)(with fat anno).pdf",height = 6,width =6)
ggplot(data=go_up_and_down, aes(x = Rho,y = -log10(Rho_P), colour=iffat)) +
  geom_point(alpha=0.8, size=3.5,col="Grey")+
  geom_point(data=go_up_and_down[go_up_and_down$iffat %in% "Yes",],alpha=0.4, size=3.5,col="#D6604D")+
  labs(x="Spearman Rho",
       y="-log10 (adjusted P)")+
  theme_bw()+geom_text_repel(data=go_up_and_down[go_up_and_down$iffat %in% "Yes",],
                             aes(label=Description),col="black",alpha = 0.8,size = 2)+xlim(-0.5,0.5)+ylim(0,6)
dev.off()
