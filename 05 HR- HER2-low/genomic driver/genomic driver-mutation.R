table(FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE)
# HER2_0      HER2_low HER2_positive 
# 91           434           182 

HR_neg_HER2_0_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HR_status %in% c("Negative")& FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_0"),"PatientCode"]
HR_neg_HER2_low_list_basal<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HR_status %in% c("Negative") & FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_low") & FUSCC_HER2_low_project_Cohort.Info$PAM50 %in% c("Basal"),"PatientCode"]
HR_neg_HER2_low_list_nonbasal<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HR_status %in% c("Negative")& FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_low") & FUSCC_HER2_low_project_Cohort.Info$PAM50 %in% c("Her2","LumA","LumB","Normal"),"PatientCode"]
HR_neg_HER2_low_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HR_status %in% c("Negative")& FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_low"),"PatientCode"]
HR_neg_HER2_pos_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HR_status %in% c("Negative")& FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_positive"),"PatientCode"]


group_list<-c(rep("HER2_0",length(HR_neg_HER2_0_list)),
              rep("HER2_low_basal",length(HR_neg_HER2_low_list_basal)),
              rep("HER2_low_nonbasal",length(HR_neg_HER2_low_list_nonbasal)),
              rep("HER2_pos",length(HR_neg_HER2_pos_list)))

names(group_list)<-c(HR_neg_HER2_0_list,HR_neg_HER2_low_list_basal,HR_neg_HER2_low_list_nonbasal,HR_neg_HER2_pos_list)
group_list<-factor(group_list,levels=c("HER2_0","HER2_low_basal","HER2_low_nonbasal","HER2_pos"))
table(group_list)



##waterfall
library(maftools)
table(group_list)

maf_for_maf=FUSCC_HER2_low_project_WES_Somatic[FUSCC_HER2_low_project_WES_Somatic$Tumor_Sample_Barcode %in% names(group_list),]
genes=read.csv("CAG.csv")[,1]
genes<-genes[genes %in% unique(maf_for_maf$Hugo_Symbol)]
maf_for_maf<-maf_for_maf[maf_for_maf$Hugo_Symbol %in% genes,]


clinicalData_for_maf=FUSCC_HER2_low_project_Cohort.Info[unique(maf_for_maf$Tumor_Sample_Barcode),c("PatientCode","PAM50","HER2_low_status_RE")]
colnames(clinicalData_for_maf)[1]<-c("Tumor_Sample_Barcode")
clinicalData_for_maf$cluster<-group_list[clinicalData_for_maf[,1]]


laml = read.maf(maf = maf_for_maf,
                clinicalData = clinicalData_for_maf ,
                verbose = FALSE)

vc_cols = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(vc_cols) = c(  'Frame_Shift_Del',  'Missense_Mutation',  'Nonsense_Mutation',  'Multi_Hit',
                     'Frame_Shift_Ins',  'In_Frame_Ins',  'Splice_Site',  'In_Frame_Del')

##waterfall
pdf(file=paste("Mutation_waterfall_cancer_associated_genes_HR_neg_4_SUBTYPE.pdf"),width = 8,height = 6)
oncoplot(colors = vc_cols,maf = laml, draw_titv = TRUE,keepGeneOrder=T,writeMatrix = T,annotationOrder=c("HER2_0","HER2_low_basal","HER2_low_nonbasal","HER2_pos"),clinicalFeatures = c('cluster'),sortByAnnotation = TRUE,top = 20,removeNonMutated=F)
dev.off()


gene_fall_list<-rownames(read.csv("onco_matrix.txt",sep = "\t"))


oncoplot(colors = vc_cols,maf = laml, genes = c(gene_fall_list[1:5],"ERBB2"),writeMatrix = T, draw_titv = TRUE,keepGeneOrder=T,removeNonMutated=F,annotationOrder=c("HER2_0","HER2_low_basal","HER2_low_nonbasal","HER2_pos"),clinicalFeatures = c('cluster'),sortByAnnotation = TRUE)
gene_fall_list<-rownames(read.csv("onco_matrix.txt",sep = "\t"))

##histplot
histplot_dataframe <- read.csv("onco_matrix.txt",sep = "\t")
histplot_dataframe[histplot_dataframe == ""] <- "WT"
histplot_dataframe[histplot_dataframe == 0] <- "WT"
histplot_dataframe[!histplot_dataframe == "WT"] <- "MT"
histplot_dataframe <- as.data.frame(t(histplot_dataframe))
histplot_dataframe$cluster <- as.character(group_list[rownames(histplot_dataframe)])


col <- c("#4472C4","#D1D1D1")
names(col) <- c("MT","WT")

table_temp <- table(histplot_dataframe[,c("cluster","TP53")])
table_temp <- as.data.frame.matrix(prop.table(table_temp,1))
barplot(as.matrix(t(table_temp)),col = col,border = F,horiz = T)

##histplot
histplot_dataframe <- read.csv("onco_matrix.txt",sep = "\t")
histplot_dataframe[histplot_dataframe == ""] <- "WT"
histplot_dataframe[histplot_dataframe == 0] <- "WT"
histplot_dataframe[!histplot_dataframe == "WT"] <- "MT"
histplot_dataframe <- as.data.frame(t(histplot_dataframe))
histplot_dataframe$cluster <- as.character(group_list[rownames(histplot_dataframe)])

histplot_dataframe_m<-reshape2::melt(histplot_dataframe, id.vars = c("cluster"),
                                     measure.vars = colnames(histplot_dataframe)[-ncol(histplot_dataframe)],
                                     variable.name = c('Gene'),
                                     value.name = 'Status')

col <- c("#4472C4","#D1D1D1")
names(col) <- c("MT","WT")

table_temp <- table(histplot_dataframe_m[,c("Gene","Status","cluster")])

pdf("Matation_histplot.pdf",width = 3,height = 5)
table_temp2 <- as.data.frame.matrix(prop.table(table_temp[6:1,,"HER2_0"],1))
barplot(as.matrix(t(table_temp2)),col = col,border = F,horiz = T,space = 0.8,cex.names=0.5)

table_temp2 <- as.data.frame.matrix(prop.table(table_temp[6:1,,"HER2_low_basal"],1))
barplot(as.matrix(t(table_temp2)),col = col,border = F,horiz = T,space = 0.8,cex.names=0.5)

table_temp2 <- as.data.frame.matrix(prop.table(table_temp[6:1,,"HER2_low_nonbasal"],1))
barplot(as.matrix(t(table_temp2)),col = col,border = F,horiz = T,space = 0.8,cex.names=0.5)

table_temp2 <- as.data.frame.matrix(prop.table(table_temp[6:1,,"HER2_pos"],1))
barplot(as.matrix(t(table_temp2)),col = col,border = F,horiz = T,space = 0.8,cex.names=0.5)

dev.off()

##waterfall  seperate
##1
maf_for_maf=FUSCC_HER2_low_project_WES_Somatic[FUSCC_HER2_low_project_WES_Somatic$Tumor_Sample_Barcode %in% HR_neg_HER2_0_list,]

clinicalData_for_maf=FUSCC_HER2_low_project_Cohort.Info[unique(maf_for_maf$Tumor_Sample_Barcode),c("PatientCode","PAM50","HER2_low_status_RE")]
colnames(clinicalData_for_maf)[1]<-c("Tumor_Sample_Barcode")
laml = read.maf(maf = maf_for_maf,
                clinicalData = clinicalData_for_maf ,
                verbose = FALSE)

pdf(file=paste("Mutation_waterfall_cancer_associated_genes_HR_neg_4_SUBTYPE_HR_neg_HER2_0_list.pdf"),width = 8,height = 6)
oncoplot(colors = vc_cols,maf = laml, genes = c(gene_fall_list[1:5],"ERBB2"), draw_titv = TRUE,keepGeneOrder=T,top = 20,removeNonMutated=F)
dev.off()

##2
maf_for_maf=FUSCC_HER2_low_project_WES_Somatic[FUSCC_HER2_low_project_WES_Somatic$Tumor_Sample_Barcode %in% HR_neg_HER2_low_list_basal,]

clinicalData_for_maf=FUSCC_HER2_low_project_Cohort.Info[unique(maf_for_maf$Tumor_Sample_Barcode),c("PatientCode","PAM50","HER2_low_status_RE")]
colnames(clinicalData_for_maf)[1]<-c("Tumor_Sample_Barcode")



laml = read.maf(maf = maf_for_maf,
                clinicalData = clinicalData_for_maf ,
                verbose = FALSE)

pdf(file=paste("Mutation_waterfall_cancer_associated_genes_HR_neg_4_SUBTYPE_HR_neg_HER2_low_list_basal.pdf"),width = 8,height = 6)
oncoplot(colors = vc_cols,maf = laml, genes = c(gene_fall_list[1:5],"ERBB2"), draw_titv = TRUE,keepGeneOrder=T,top = 20,removeNonMutated=F)
dev.off()

##3
maf_for_maf=FUSCC_HER2_low_project_WES_Somatic[FUSCC_HER2_low_project_WES_Somatic$Tumor_Sample_Barcode %in% HR_neg_HER2_low_list_nonbasal,]

clinicalData_for_maf=FUSCC_HER2_low_project_Cohort.Info[unique(maf_for_maf$Tumor_Sample_Barcode),c("PatientCode","PAM50","HER2_low_status_RE")]
colnames(clinicalData_for_maf)[1]<-c("Tumor_Sample_Barcode")



laml = read.maf(maf = maf_for_maf,
                clinicalData = clinicalData_for_maf ,
                verbose = FALSE)

pdf(file=paste("Mutation_waterfall_cancer_associated_genes_HR_neg_4_SUBTYPE_HR_neg_HER2_low_list_nonbasal.pdf"),width = 8,height = 6)
oncoplot(colors = vc_cols,maf = laml, genes = c(gene_fall_list[1:5],"ERBB2"), draw_titv = TRUE,keepGeneOrder=T,top = 20,removeNonMutated=F)
dev.off()

##4
maf_for_maf=FUSCC_HER2_low_project_WES_Somatic[FUSCC_HER2_low_project_WES_Somatic$Tumor_Sample_Barcode %in% HR_neg_HER2_pos_list,]

clinicalData_for_maf=FUSCC_HER2_low_project_Cohort.Info[unique(maf_for_maf$Tumor_Sample_Barcode),c("PatientCode","PAM50","HER2_low_status_RE")]
colnames(clinicalData_for_maf)[1]<-c("Tumor_Sample_Barcode")



laml = read.maf(maf = maf_for_maf,
                clinicalData = clinicalData_for_maf ,
                verbose = FALSE)

pdf(file=paste("Mutation_waterfall_cancer_associated_genes_HR_neg_4_SUBTYPE_HR_neg_HER2_pos_listl.pdf"),width = 8,height = 6)
oncoplot(colors = vc_cols,maf = laml, genes = c(gene_fall_list[1:5],"ERBB2"), draw_titv = TRUE,keepGeneOrder=T,top = 20,removeNonMutated=F)
dev.off()

# mutation comparation

library(maftools)

maf_for_maf1=FUSCC_HER2_low_project_WES_Somatic[FUSCC_HER2_low_project_WES_Somatic$Tumor_Sample_Barcode %in% HR_neg_HER2_0_list,]

clinicalData_for_maf=FUSCC_HER2_low_project_Cohort.Info[unique(maf_for_maf1$Tumor_Sample_Barcode),c("PatientCode","PAM50","HER2_low_status_RE")]
colnames(clinicalData_for_maf)[1]<-c("Tumor_Sample_Barcode")



laml1 = read.maf(maf = maf_for_maf1,
                 clinicalData = clinicalData_for_maf ,
                 verbose = FALSE)


maf_for_maf2=FUSCC_HER2_low_project_WES_Somatic[FUSCC_HER2_low_project_WES_Somatic$Tumor_Sample_Barcode %in% HR_neg_HER2_low_list_nonbasal,]

clinicalData_for_maf=FUSCC_HER2_low_project_Cohort.Info[unique(maf_for_maf2$Tumor_Sample_Barcode),c("PatientCode","PAM50","HER2_low_status_RE")]
colnames(clinicalData_for_maf)[1]<-c("Tumor_Sample_Barcode")



laml2 = read.maf(maf = maf_for_maf2,
                 clinicalData = clinicalData_for_maf ,
                 verbose = FALSE)


pt.vs.rt <- mafCompare(m1 = laml1, m2 = laml2, m1Name = 'HR_neg_HER2_0_list', m2Name = 'HR_neg_HER2_low_list_nonbasal', minMut = 5)
print(pt.vs.rt)
# $results
#    Hugo_Symbol HR_neg_HER2_0_list HR_neg_HER2_low_list_nonbasal        pval         or      ci.up      ci.low     adjPval
# 1:      PIK3CA                  2                             8 0.001706689 0.06952037  0.4871188 0.005539091 0.003413379
# 2:        TP53                 19                             7 0.050497285 5.13408757 40.9227959 0.828707979 0.050497285
# 
# $SampleSummary
# Cohort SampleSize
# 1:            HR_neg_HER2_0_list         22
# 2: HR_neg_HER2_low_list_nonbasal         13
forestPlot(mafCompareRes = pt.vs.rt, pVal = 0.05, color = c('royalblue', 'maroon'), geneFontSize = 0.8)


# mutation comparation
library(maftools)

maf_for_maf1=FUSCC_HER2_low_project_WES_Somatic[FUSCC_HER2_low_project_WES_Somatic$Tumor_Sample_Barcode %in% HR_neg_HER2_low_list_basal,]

clinicalData_for_maf=FUSCC_HER2_low_project_Cohort.Info[unique(maf_for_maf1$Tumor_Sample_Barcode),c("PatientCode","PAM50","HER2_low_status_RE")]
colnames(clinicalData_for_maf)[1]<-c("Tumor_Sample_Barcode")



laml1 = read.maf(maf = maf_for_maf1,
                 clinicalData = clinicalData_for_maf ,
                 verbose = FALSE)


maf_for_maf2=FUSCC_HER2_low_project_WES_Somatic[FUSCC_HER2_low_project_WES_Somatic$Tumor_Sample_Barcode %in% HR_neg_HER2_low_list_nonbasal,]

clinicalData_for_maf=FUSCC_HER2_low_project_Cohort.Info[unique(maf_for_maf2$Tumor_Sample_Barcode),c("PatientCode","PAM50","HER2_low_status_RE")]
colnames(clinicalData_for_maf)[1]<-c("Tumor_Sample_Barcode")



laml2 = read.maf(maf = maf_for_maf2,
                 clinicalData = clinicalData_for_maf ,
                 verbose = FALSE)


pt.vs.rt <- mafCompare(m1 = laml1, m2 = laml2, m1Name = 'HR_neg_HER2_low_list_basal', m2Name = 'HR_neg_HER2_low_list_nonbasal', minMut = 5)
print(pt.vs.rt)
# $results
# Hugo_Symbol HR_neg_HER2_low_list_basal HR_neg_HER2_low_list_nonbasal         pval        or      ci.up     ci.low     adjPval
# 1:      PIK3CA                          4                             8 0.0007475214 0.0817753  0.4324626 0.01246389 0.004485129
# 2:        TP53                         32                             7 0.0233395462 5.2535053 29.3474686 1.02118655 0.070018639
# 3:    ARHGAP30                          5                             0 0.3087253865       Inf        Inf 0.32236404 0.370470464
# 4:     CTTNBP2                          5                             0 0.3087253865       Inf        Inf 0.32236404 0.370470464
# 5:      MYO15A                          5                             0 0.3087253865       Inf        Inf 0.32236404 0.370470464
# 6:       CDH23                          6                             2 1.0000000000 1.0632033 12.2970702 0.15698583 1.000000000
# 
# $SampleSummary
# Cohort SampleSize
# 1:    HR_neg_HER2_low_list_basal         37
# 2: HR_neg_HER2_low_list_nonbasal         13

forestPlot(mafCompareRes = pt.vs.rt, pVal = 0.05, color = c('royalblue', 'maroon'), geneFontSize = 0.8)


# mutation comparation
library(maftools)

maf_for_maf1=FUSCC_HER2_low_project_WES_Somatic[FUSCC_HER2_low_project_WES_Somatic$Tumor_Sample_Barcode %in% HR_neg_HER2_low_list_nonbasal,]

clinicalData_for_maf=FUSCC_HER2_low_project_Cohort.Info[unique(maf_for_maf1$Tumor_Sample_Barcode),c("PatientCode","PAM50","HER2_low_status_RE")]
colnames(clinicalData_for_maf)[1]<-c("Tumor_Sample_Barcode")



laml1 = read.maf(maf = maf_for_maf1,
                 clinicalData = clinicalData_for_maf ,
                 verbose = FALSE)


maf_for_maf2=FUSCC_HER2_low_project_WES_Somatic[FUSCC_HER2_low_project_WES_Somatic$Tumor_Sample_Barcode %in% HR_neg_HER2_pos_list,]

clinicalData_for_maf=FUSCC_HER2_low_project_Cohort.Info[unique(maf_for_maf2$Tumor_Sample_Barcode),c("PatientCode","PAM50","HER2_low_status_RE")]
colnames(clinicalData_for_maf)[1]<-c("Tumor_Sample_Barcode")



laml2 = read.maf(maf = maf_for_maf2,
                 clinicalData = clinicalData_for_maf ,
                 verbose = FALSE)


pt.vs.rt <- mafCompare(m1 = laml1, m2 = laml2, m1Name = 'HR_neg_HER2_low_list_nonbasal', m2Name = 'HR_neg_HER2_pos_list', minMut = 5)
print(pt.vs.rt)
# $results
# Hugo_Symbol HR_neg_HER2_low_list_nonbasal HR_neg_HER2_pos_list       pval        or     ci.up     ci.low   adjPval
# 1:      PIK3CA                             8                   18 0.02415511 4.1770150 18.537229 1.04630564 0.1932409
# 2:       ALMS1                             0                    5 0.58427715 0.0000000  5.752461 0.00000000 0.9348434
# 3:      ARID1A                             0                    5 0.58427715 0.0000000  5.752461 0.00000000 0.9348434
# 4:       EPHA5                             0                    5 0.58427715 0.0000000  5.752461 0.00000000 0.9348434
# 5:        UTRN                             0                    5 0.58427715 0.0000000  5.752461 0.00000000 0.9348434
# 6:        TP53                             7                   39 0.76535620 0.8099105  3.268436 0.20674886 1.0000000
# 7:       DNAH6                             1                    5 1.00000000 1.0164538 10.360953 0.01986078 1.0000000
# 8:      SPHKAP                             1                    6 1.00000000 0.8351631  7.940409 0.01674410 1.0000000
# 
# $SampleSummary
# Cohort SampleSize
# 1: HR_neg_HER2_low_list_nonbasal         13
# 2:          HR_neg_HER2_pos_list         66

forestPlot(mafCompareRes = pt.vs.rt, pVal = 0.05, color = c('royalblue', 'maroon'), geneFontSize = 0.8)

## internal compartion of HER2-low nonbasal (gene set)
load("../../04 molecular landscape/GSVA/RNA_C2.Rdata")
exprSet<-esmicro

HER2_low_nonbasal_WT<-intersect(HR_neg_HER2_low_list_nonbasal,rownames(histplot_dataframe)[histplot_dataframe[,"PIK3CA"] %in% c("WT")])
HER2_low_nonbasal_MT<-intersect(HR_neg_HER2_low_list_nonbasal,rownames(histplot_dataframe)[histplot_dataframe[,"PIK3CA"] %in% c("MT")])

plot_list<-c(rep("HER2_low_nonbasal_WT",length(HER2_low_nonbasal_WT)),
             rep("HER2_low_nonbasal_MT",length(HER2_low_nonbasal_MT)))
names(plot_list)<-c(HER2_low_nonbasal_WT,HER2_low_nonbasal_MT)
names(plot_list) <- paste(names(plot_list),"_RNA_T",sep = "")
plot_list<-factor(plot_list)
plot_list<-plot_list[names(plot_list) %in% colnames(exprSet)]

my_comparisons <- list( c("HER2_low_nonbasal_WT", "HER2_low_nonbasal_MT"))

gene<-c("REACTOME_PI3K_AKT_ACTIVATION")
col <- c("#E41A1C","#4DAF4A","#984EA3","#377EB8")
my_comparisons <- list( c("HER2_low_nonbasal_WT", "HER2_low_nonbasal_MT"))
dataset=data.frame(mRNA=as.matrix(exprSet)[gene,names(plot_list)],group=plot_list)
library(ggplot2)
library(ggpubr)
pdf(paste("boxplot_GSVA_",gene,"_HER2_low_nonbasal_WT_vs_MT.pdf",sep=""),width =2,height = 5)
ggboxplot(dataset, x = "group", y = "mRNA",ylab = "GSVA score",
          color = "group", palette = col,
          add = "jitter", title = gene,xlab  = FALSE)+stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")
dev.off()

gene<-c("REACTOME_PI3K_EVENTS_IN_ERBB2_SIGNALING")
col <- c("#E41A1C","#4DAF4A","#984EA3","#377EB8")
my_comparisons <- list( c("HER2_low_nonbasal_WT", "HER2_low_nonbasal_MT"))
dataset=data.frame(mRNA=as.matrix(exprSet)[gene,names(plot_list)],group=plot_list)
library(ggplot2)
library(ggpubr)
pdf(paste("boxplot_GSVA_",gene,"_HER2_low_nonbasal_WT_vs_MT.pdf",sep=""),width =2,height = 5)
ggboxplot(dataset, x = "group", y = "mRNA",ylab = "GSVA score",
          color = "group", palette = col,
          add = "jitter", title = gene,xlab  = FALSE)+stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")
dev.off()


