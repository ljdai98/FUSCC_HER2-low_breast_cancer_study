table(FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE)
# HER2_0      HER2_low HER2_positive 
# 91           434           182 


##############################################################
HR_neg_HER2_0_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$Clinical_Subtype %in% c("HR-HER2+","TNBC") & FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_0"),"PatientCode"]
HR_neg_HER2_low_list_basal<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$Clinical_Subtype %in% c("HR-HER2+","TNBC") & FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_low") & FUSCC_HER2_low_project_Cohort.Info$PAM50 %in% c("Basal"),"PatientCode"]
HR_neg_HER2_low_list_nonbasal<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$Clinical_Subtype %in% c("HR-HER2+","TNBC") & FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_low") & FUSCC_HER2_low_project_Cohort.Info$PAM50 %in% c("Her2","LumA","LumB","Normal"),"PatientCode"]
HR_neg_HER2_low_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$Clinical_Subtype %in% c("HR-HER2+","TNBC") & FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_low"),"PatientCode"]

HR_neg_HER2_pos_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$Clinical_Subtype %in% c("HR-HER2+","TNBC") & FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_positive"),"PatientCode"]

group_list<-c(rep("HER2_0",length(HR_neg_HER2_0_list)),
              rep("HER2_low_basal",length(HR_neg_HER2_low_list_basal)),
              rep("HER2_low_nonbasal",length(HR_neg_HER2_low_list_nonbasal)),
              rep("HER2_pos",length(HR_neg_HER2_pos_list)))

names(group_list)<-c(HR_neg_HER2_0_list,HR_neg_HER2_low_list_basal,HR_neg_HER2_low_list_nonbasal,HR_neg_HER2_pos_list)
group_list<-factor(group_list,levels=c("HER2_0","HER2_low_basal","HER2_low_nonbasal","HER2_pos"))
table(group_list)



##3D PCA RNA
library(pca3d)
library(ggfortify)
plot_names<-paste(names(group_list),"_RNA_T",sep="")
plot_names<-plot_names[plot_names %in% colnames(FUSCC_HER2_low_project_RNA_seq_log2FPKM)]
plot_list<-group_list[substr(plot_names,1,4)]

exp<-FUSCC_HER2_low_project_RNA_seq_log2FPKM[,plot_names]
dat<-cbind(as.character(plot_list),t(exp))
colnames(dat)[1]<-c("label")
dat<-t(na.omit(t(dat)))

dim(dat)

out_pca <- prcomp(apply(dat[,-1],2,as.numeric))
dat[,1]<-as.character(dat[,1])


for_review <- data.frame(group_list=dat[,1],out_pca[["x"]][,c("PC1","PC2","PC3")])

color_pca <- c("#DC978D","#4DAF4A","#984EA3","#83B5D4")
names(color_pca) <- c("HER2_0","HER2_low_basal","HER2_low_nonbasal","HER2_pos")

for_review$color <- color_pca[for_review$group_list]

library(scatterplot3d )

pdf("3D_PCA_RNA.pdf",width = 5,height = 5)
s3d <- scatterplot3d(x = for_review$PC1,
                     y = for_review$PC2,
                     z = for_review$PC3,
                     xlab = "PC1", ylab = "PC2", zlab = "PC3",
                     pch = 21,
                     bg=for_review$color,
                     angle = 120,
                     scale.y = 0.7,
                     grid=T,
                     cex.symbols = 1.5,
                     col.axis = "#444444",col.grid = "#CCCCCC",col.lab="black")
dev.off()

##normalized inter-center distance
PCA_result <- data.frame(out_pca$x)
rownames(PCA_result) <- names(dat[,1])

PCA_result$group <- dat[,1]

coordinate <- aggregate(PCA_result[,c('PC1','PC2','PC3')], by=list(PCA_result$group),mean)
rownames(coordinate) <- coordinate$Group.1
coordinate <- coordinate[,-1]

dist_result <- dist(coordinate,p=2)
dist_result
#                     HER2_0 HER2_low_basal HER2_low_nonbasal
# HER2_low_basal    10.60286                                 
# HER2_low_nonbasal 50.42615       52.90672                  
# HER2_pos          50.56113       49.11357          20.24019

dist_result <- round(dist_result/dist_result[1],1)
dist_result
#                     HER2_0 HER2_low_basal HER2_low_nonbasal
# HER2_low_basal      1.00                                 
# HER2_low_nonbasal   4.76           4.99                  
# HER2_pos            4.77           4.63              1.91


pdf("3D_PCA_RNA_anno.pdf",width = 3,height = 3)

lengh_seq <- seq(1,max(dist_result),0.1)
length_col <- length(lengh_seq)
cols<-c('#F2F2F2','#DC2418')
pal<-colorRampPalette(cols)(length_col)

image(x=1:length_col,y=1,z=as.matrix(1:length_col),col=pal)

names(pal) <- rev(lengh_seq)
pal <- pal[as.character(c(round(dist_result,1)))]

plot(2, 2, col = "white", xlab = "", ylab = "")
segments(1.5,2.5,1.5,1.5,lwd=5,col=pal[1])
segments(1.5,2.5,2.5,1.5,lwd=5,col=pal[2])
segments(1.5,2.5,2.5,2.5,lwd=5,col=pal[3])
segments(1.5,1.5,2.5,1.5,lwd=5,col=pal[4])
segments(1.5,1.5,2.5,2.5,lwd=5,col=pal[5])
segments(2.5,1.5,2.5,2.5,lwd=5,col=pal[6])

dev.off()
