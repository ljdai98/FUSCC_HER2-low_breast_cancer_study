table(FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE)
# HER2_0      HER2_low HER2_positive 
# 91           434           182 


list_a_samplenames<-rownames(FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_0") & 
                                                                  FUSCC_HER2_low_project_Cohort.Info$HR_status %in% c("Positive"),])
list_a_typename<-c("HR_pos_HER2_0")

list_b_samplenames<-rownames(FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_low") & 
                                                                  FUSCC_HER2_low_project_Cohort.Info$HR_status %in% c("Positive"),])
list_b_typename<-c("HR_pos_HER2_low")

list_c_samplenames<-rownames(FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_positive") & 
                                                                  FUSCC_HER2_low_project_Cohort.Info$HR_status %in% c("Positive"),])
list_c_typename<-c("HR_pos_HER2_pos")


group_list<-c(rep(list_a_typename,length(list_a_samplenames)),
              rep(list_b_typename,length(list_b_samplenames)),
              rep(list_c_typename,length(list_c_samplenames)))
names(group_list)<-c(list_a_samplenames,list_b_samplenames,list_c_samplenames)
group_list<-factor(group_list,levels = c(list_a_typename,list_b_typename,list_c_typename))
table(group_list)
# HR_pos_HER2_0 HR_pos_HER2_low HR_pos_HER2_pos 
#  63             361             101



## forest HER2-0 vs HER2-low

library(maftools)
maf_for_maf1=FUSCC_HER2_low_project_WES_Somatic[FUSCC_HER2_low_project_WES_Somatic$Tumor_Sample_Barcode %in% list_a_samplenames,]
genes=read.csv("CAG.csv")[,1]
genes<-genes[genes %in% unique(maf_for_maf1$Hugo_Symbol)]
maf_for_maf1<-maf_for_maf1[maf_for_maf1$Hugo_Symbol %in% genes,]

clinicalData_for_maf=FUSCC_HER2_low_project_Cohort.Info[unique(maf_for_maf1$Tumor_Sample_Barcode),c("PatientCode","PAM50","HER2_low_status_RE")]
colnames(clinicalData_for_maf)[1]<-c("Tumor_Sample_Barcode")


laml1 = read.maf(maf = maf_for_maf1,
                 clinicalData = clinicalData_for_maf ,
                 verbose = FALSE)

maf_for_maf2=FUSCC_HER2_low_project_WES_Somatic[FUSCC_HER2_low_project_WES_Somatic$Tumor_Sample_Barcode %in% list_b_samplenames,]
genes=read.csv("CAG.csv")[,1]
genes<-genes[genes %in% unique(maf_for_maf2$Hugo_Symbol)]
maf_for_maf2<-maf_for_maf2[maf_for_maf2$Hugo_Symbol %in% genes,]


clinicalData_for_maf=FUSCC_HER2_low_project_Cohort.Info[unique(maf_for_maf2$Tumor_Sample_Barcode),c("PatientCode","PAM50","HER2_low_status_RE")]
colnames(clinicalData_for_maf)[1]<-c("Tumor_Sample_Barcode")


laml2 = read.maf(maf = maf_for_maf2,
                 clinicalData = clinicalData_for_maf ,
                 verbose = FALSE)

pdf(file=paste("Mutation_comparation_forest_plot_",list_a_typename,"_vs_",list_b_typename,"_cancer_associated_genes.pdf",sep = ""))
pt.vs.rt <- mafCompare(m1 = laml1, m2 = laml2, m1Name = list_a_typename, m2Name = list_b_typename, minMut = 5)
print(pt.vs.rt)
# $results
# Hugo_Symbol HR_pos_HER2_0 HR_pos_HER2_low       pval        or     ci.up      ci.low adjPval
# 1:        TP53            20              77 0.09699726 1.7188435  3.272548 0.884697318       1
# 2:       SPTA1             3               7 0.18031096 2.4808138 11.301386 0.401293306       1
# $SampleSummary
# Cohort SampleSize
# 1:   HR_pos_HER2_0         55
# 2: HR_pos_HER2_low        309

forestPlot(mafCompareRes = pt.vs.rt, pVal = 0.05, color = c('royalblue', 'maroon'), geneFontSize = 0.8)
dev.off()

sig_genes1 <- as.data.frame(pt.vs.rt[["results"]])
rownames(sig_genes1) <- sig_genes1$Hugo_Symbol
sig_genes1 <- rownames(sig_genes1[sig_genes1$pval < 0.05,])


## forest HER2-pos vs HER2-low

library(maftools)
maf_for_maf1=FUSCC_HER2_low_project_WES_Somatic[FUSCC_HER2_low_project_WES_Somatic$Tumor_Sample_Barcode %in% list_c_samplenames,]
genes=read.csv("CAG.csv")[,1]
genes<-genes[genes %in% unique(maf_for_maf1$Hugo_Symbol)]
maf_for_maf1<-maf_for_maf1[maf_for_maf1$Hugo_Symbol %in% genes,]

clinicalData_for_maf=FUSCC_HER2_low_project_Cohort.Info[unique(maf_for_maf1$Tumor_Sample_Barcode),c("PatientCode","PAM50","HER2_low_status_RE")]
colnames(clinicalData_for_maf)[1]<-c("Tumor_Sample_Barcode")


laml1 = read.maf(maf = maf_for_maf1,
                 clinicalData = clinicalData_for_maf ,
                 verbose = FALSE)

maf_for_maf2=FUSCC_HER2_low_project_WES_Somatic[FUSCC_HER2_low_project_WES_Somatic$Tumor_Sample_Barcode %in% list_b_samplenames,]
genes=read.csv("CAG.csv")[,1]
genes<-genes[genes %in% unique(maf_for_maf2$Hugo_Symbol)]
maf_for_maf2<-maf_for_maf2[maf_for_maf2$Hugo_Symbol %in% genes,]


clinicalData_for_maf=FUSCC_HER2_low_project_Cohort.Info[unique(maf_for_maf2$Tumor_Sample_Barcode),c("PatientCode","PAM50","HER2_low_status_RE")]
colnames(clinicalData_for_maf)[1]<-c("Tumor_Sample_Barcode")


laml2 = read.maf(maf = maf_for_maf2,
                 clinicalData = clinicalData_for_maf ,
                 verbose = FALSE)

pdf(file=paste("Mutation_comparation_forest_plot_",list_c_typename,"_vs_",list_b_typename,"_cancer_associated_genes.pdf",sep = ""))
pt.vs.rt <- mafCompare(m1 = laml1, m2 = laml2, m1Name = list_c_typename, m2Name = list_b_typename, minMut = 5)
print(pt.vs.rt)
# $results
# Hugo_Symbol HR_pos_HER2_pos HR_pos_HER2_low         pval        or      ci.up      ci.low      adjPval
# 1:        TP53              45              77 4.847529e-07 3.6500111  6.2712427 2.139779395 2.035962e-05
# 2:      MAP3K1               0              42 3.499923e-05 0.0000000  0.3083588 0.000000000 7.349838e-04
# 3:         NF1               7               6 8.099117e-03 4.6881077 17.4219886 1.306750522 1.133876e-01
# 4:        AKT1               1              24 3.846443e-02 0.1470185  0.9302323 0.003526143 4.038765e-0
# $SampleSummary
# Cohort SampleSize
# 1: HR_pos_HER2_pos         82
# 2: HR_pos_HER2_low        309
forestPlot(mafCompareRes = pt.vs.rt, pVal = 0.05, color = c('royalblue', 'maroon'), geneFontSize = 0.8)
dev.off()

sig_genes2 <- as.data.frame(pt.vs.rt[["results"]])
rownames(sig_genes2) <- sig_genes2$Hugo_Symbol
sig_genes2 <- rownames(sig_genes2[sig_genes2$pval < 0.05,])


##########
vc_cols = RColorBrewer::brewer.pal(n = 9, name = 'Paired')
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del',"Translation_Start_Site"
)


##waterfall

library(maftools)

table(group_list)

maf_for_maf=FUSCC_HER2_low_project_WES_Somatic[FUSCC_HER2_low_project_WES_Somatic$Tumor_Sample_Barcode %in% names(group_list),]
genes=read.csv("CAG.csv")[,1]
genes<-genes[genes %in% unique(maf_for_maf$Hugo_Symbol)]
maf_for_maf<-maf_for_maf[maf_for_maf$Hugo_Symbol %in% genes,]


clinicalData_for_maf=FUSCC_HER2_low_project_Cohort.Info[unique(maf_for_maf$Tumor_Sample_Barcode),c("PatientCode","PAM50","HER2_low_status_RE")]
colnames(clinicalData_for_maf)[1]<-c("Tumor_Sample_Barcode")
clinicalData_for_maf$Category<-group_list[clinicalData_for_maf[,1]]


laml = read.maf(maf = maf_for_maf,
                clinicalData = clinicalData_for_maf ,
                verbose = FALSE)

##waterfall
pdf(file=paste("Mutation_waterfall_HR+_HER2_0_VS_HER2_LOW_VS_HER2_POS.pdf"),width = 8,height = 6)
oncoplot(colors = vc_cols,maf = laml, draw_titv = TRUE,keepGeneOrder=T,writeMatrix = T,clinicalFeatures = c('Category'),annotationOrder=c("HR_pos_HER2_0","HR_pos_HER2_low","HR_pos_HER2_pos"),sortByAnnotation = TRUE,top = 20,removeNonMutated=F)#writeMatrix可自动产生一个对应的txt文件
dev.off()

file.rename(from = "onco_matrix.txt",to="onco_matrix_all_HR+.txt")

##waterfallseparate
gene_fall_list<-rownames(read.csv("onco_matrix_all_HR+.txt",sep = "\t")[1])[1:10]

other_sig_genes <- unique(c(sig_genes1,sig_genes2))
other_sig_genes <- other_sig_genes[!other_sig_genes %in% gene_fall_list]

gene_fall_list <- c(gene_fall_list[1:10],other_sig_genes)

oncoplot(colors = vc_cols,maf = laml, genes =gene_fall_list,writeMatrix = T, draw_titv = TRUE,keepGeneOrder=T,removeNonMutated=F,clinicalFeatures = c('Category'),annotationOrder=c("HR_pos_HER2_0","HR_pos_HER2_low","HR_pos_HER2_pos"),sortByAnnotation = TRUE)#writeMatrix可自动产生一个对应的txt文件


##1
maf_for_maf=FUSCC_HER2_low_project_WES_Somatic[FUSCC_HER2_low_project_WES_Somatic$Tumor_Sample_Barcode %in% list_a_samplenames,]

clinicalData_for_maf=FUSCC_HER2_low_project_Cohort.Info[unique(maf_for_maf$Tumor_Sample_Barcode),c("PatientCode","PAM50","HER2_low_status_RE")]
colnames(clinicalData_for_maf)[1]<-c("Tumor_Sample_Barcode")

laml = read.maf(maf = maf_for_maf,
                clinicalData = clinicalData_for_maf ,
                verbose = FALSE)

pdf(file=paste("Mutation_waterfall_HR+_HER2_0_in_total_order.pdf.pdf"),width = 8,height = 6)
oncoplot(colors = vc_cols,maf = laml, genes =gene_fall_list, draw_titv = TRUE,keepGeneOrder=T,removeNonMutated=F)#writeMatrix可自动产生一个对应的txt文件
dev.off()


##2
maf_for_maf=FUSCC_HER2_low_project_WES_Somatic[FUSCC_HER2_low_project_WES_Somatic$Tumor_Sample_Barcode %in% list_b_samplenames,]

clinicalData_for_maf=FUSCC_HER2_low_project_Cohort.Info[unique(maf_for_maf$Tumor_Sample_Barcode),c("PatientCode","PAM50","HER2_low_status_RE")]
colnames(clinicalData_for_maf)[1]<-c("Tumor_Sample_Barcode")



laml = read.maf(maf = maf_for_maf,
                clinicalData = clinicalData_for_maf ,
                verbose = FALSE)

pdf(file=paste("Mutation_waterfall_HR+_HER2_low_in_total_order.pdf.pdf"),width = 8,height = 6)
oncoplot(colors = vc_cols,maf = laml, genes =gene_fall_list, draw_titv = TRUE,keepGeneOrder=T,removeNonMutated=F)#writeMatrix可自动产生一个对应的txt文件
dev.off()

##3
maf_for_maf=FUSCC_HER2_low_project_WES_Somatic[FUSCC_HER2_low_project_WES_Somatic$Tumor_Sample_Barcode %in% list_c_samplenames,]

clinicalData_for_maf=FUSCC_HER2_low_project_Cohort.Info[unique(maf_for_maf$Tumor_Sample_Barcode),c("PatientCode","PAM50","HER2_low_status_RE")]
colnames(clinicalData_for_maf)[1]<-c("Tumor_Sample_Barcode")



laml = read.maf(maf = maf_for_maf,
                clinicalData = clinicalData_for_maf ,
                verbose = FALSE)

pdf(file=paste("Mutation_waterfall_HR+_HER2_pos_in_total_order.pdf.pdf"),width = 8,height = 6)
oncoplot(colors = vc_cols,maf = laml, genes =gene_fall_list, draw_titv = TRUE,keepGeneOrder=T,removeNonMutated=F)#writeMatrix可自动产生一个对应的txt文件
dev.off()

