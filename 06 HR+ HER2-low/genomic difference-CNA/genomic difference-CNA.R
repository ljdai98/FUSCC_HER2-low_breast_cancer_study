table(FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE)
# HER2_0      HER2_low HER2_positive 
# 91           434           182 



##group_list

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
# group_list
# HR_pos_HER2_0 HR_pos_HER2_low HR_pos_HER2_pos 
# 63             361             101 

CBCGA_GISTICpeaks<-rbind(FUSCC_HER2_low_project_GISTICpeaks.amp.thre,-FUSCC_HER2_low_project_GISTICpeaks.del.thre)

cnv.lista.peaks = CBCGA_GISTICpeaks[,list_a_samplenames[list_a_samplenames %in% colnames(CBCGA_GISTICpeaks)]]
cnv.listb.peaks = CBCGA_GISTICpeaks[,list_b_samplenames[list_b_samplenames %in% colnames(CBCGA_GISTICpeaks)]]
cnv.listc.peaks = CBCGA_GISTICpeaks[,list_c_samplenames[list_c_samplenames %in% colnames(CBCGA_GISTICpeaks)]]

##common genes
commonPeak = intersect(intersect(rownames(cnv.listb.peaks),rownames(cnv.lista.peaks)),rownames(cnv.listc.peaks))
cnv.lista.peaks = cnv.lista.peaks[commonPeak,] ; cnv.listb.peaks = cnv.listb.peaks[commonPeak,];cnv.listc.peaks = cnv.listc.peaks[commonPeak,]

###
ls = list(cnv.lista.peaks,cnv.listb.peaks,cnv.listc.peaks)
names(ls) = c(list_a_typename,list_b_typename,list_c_typename)
for (i in 1:3){
  cnv = ls[[i]]
  
  del = -apply(cnv,1, function(x){table(x)["-2"]})
  del[is.na(del)] = 0
  loss = -apply(cnv,1, function(x){table(x)["-1"]})
  loss[is.na(loss)] = 0
  #norm = apply(METABRIC_cnv_f4chr,1, function(x){table(x)["0"]})
  #norm[is.na(norm)] = 0
  gain = apply(cnv,1, function(x){table(x)["1"]})
  gain[is.na(gain)] = 0
  amp = apply(cnv,1, function(x){table(x)["2"]})
  amp[is.na(amp)] = 0
  
  cnv_locus_order = cbind(gain,amp,
                          #norm,
                          loss,del)
  cnv_locus_order_freq = cnv_locus_order/ncol(cnv)
  cnv_locus_order_freq = t(cnv_locus_order_freq)
  
  names<-colnames(cnv_locus_order_freq)
  pdf(paste0("CNV_histplot_peaks_",names(ls)[i],".pdf"),height =  4, width =  10)
  barplot(cnv_locus_order_freq[1:2,],col= c("#E6B3B3","#E21A21"),border = NA,ylim = c(-1,1), lwd  = 2,names.arg=as.character(as.data.frame(strsplit(names,"_"))[3,]),las=3,cex.names=0.5);
  par(new= T)
  barplot(cnv_locus_order_freq[3:4,],col= c("#B9C2D8","#3653A5"),border = NA,ylim = c(-1,1), lwd  = 2,axisnames = F)
  par(new= F)
  dev.off()
  
}


## GISTIC peaks DCG HER2-0 vs HER2-low

mat1 = as.data.frame(t(abs(cnv.lista.peaks[,]) ))
mat2 = as.data.frame(t(abs(cnv.listb.peaks[,])))

mat1$group <- list_a_typename
mat2$group <- list_b_typename
mat_combine <- rbind(mat1,mat2)
dim(mat_combine)

mat_combine_gather <- reshape2::melt(mat_combine, id.vars = c("group"), 
                                     measure.vars = colnames(mat_combine)[-ncol(mat_combine)],
                                     variable.name = c('Peak'),
                                     value.name = 'Copy_number')
dim(mat_combine_gather)
mat_combine_xtabs<- xtabs(  ~ group+Copy_number+Peak,data=mat_combine_gather)


fisher.res = as.data.frame(matrix(nrow = length(commonPeak),ncol = 8))
rownames(fisher.res) = commonPeak
colnames(fisher.res)= c(paste(list_a_typename,"GISTIC_0",sep="_"),paste(list_a_typename,"GISTIC_1",sep="_"),paste(list_a_typename,"GISTIC_2",sep="_"),paste(list_b_typename,"GISTIC_0",sep="_"),paste(list_b_typename,"GISTIC_1",sep="_"),paste(list_b_typename,"GISTIC_2",sep="_"),"p.value","adj.p")


for ( i in 1:length(commonPeak) ){

  fisher.res[commonPeak[i],c(1:6)] = c(mat_combine_xtabs[list_a_typename,,commonPeak[i]],mat_combine_xtabs[list_b_typename,,commonPeak[i]])
  fisher.res[commonPeak[i],"p.value"] = fisher.test(mat_combine_xtabs[,,commonPeak[i]])[["p.value"]]
  if (i/10 == round(i/10)) {print(i)}
}

fisher.res$adj.p <- p.adjust(fisher.res$p.value ,method = "fdr")

# write.csv(fisher.res, file = paste("CNA_comparation_peaks_fisher.exact_",list_a_typename,"_vs_",list_b_typename,".csv",sep = ""))

CNV_FDR_criteria<-0.05
names<-rownames(fisher.res)
pdf(paste("CNA_comparation_peaks_histplot_",list_a_typename,"_vs_",list_b_typename,"_FDR_",CNV_FDR_criteria,".pdf",sep = ""),height =  4, width =  10)
barplot(-log10(fisher.res$adj.p),col = "#EAAC53",border = NA, axisnames = F,lwd  = 2,ylim = c(-4,4))
par(new= T)
barplot(-log10(1/fisher.res$adj.p),col = "#B3D486",border = NA, lwd  = 2,ylim = c(-4,4),names.arg=as.character(as.data.frame(strsplit(names,"_"))[3,]),las=3,cex.names=0.5)

abline(h = -log10(CNV_FDR_criteria), lty = 1,lwd  = 1) 
abline(h = -log10(1/CNV_FDR_criteria), lty = 1,lwd  = 1) 
par(new= F)
dev.off()


##GISTIC peaks DCG HER2-pos vs HER2-low 

mat1 = as.data.frame(t(abs(cnv.listc.peaks[,]) ))
mat2 = as.data.frame(t(abs(cnv.listb.peaks[,])))

mat1$group <- list_c_typename
mat2$group <- list_b_typename
mat_combine <- rbind(mat1,mat2)
dim(mat_combine)

mat_combine_gather <- reshape2::melt(mat_combine, id.vars = c("group"), 
                                     measure.vars = colnames(mat_combine)[-ncol(mat_combine)],
                                     variable.name = c('Peak'),
                                     value.name = 'Copy_number')
dim(mat_combine_gather)
mat_combine_xtabs<- xtabs(  ~ group+Copy_number+Peak,data=mat_combine_gather)


fisher.res = as.data.frame(matrix(nrow = length(commonPeak),ncol = 8))
rownames(fisher.res) = commonPeak
colnames(fisher.res)= c(paste(list_c_typename,"GISTIC_0",sep="_"),paste(list_c_typename,"GISTIC_1",sep="_"),paste(list_c_typename,"GISTIC_2",sep="_"),paste(list_b_typename,"GISTIC_0",sep="_"),paste(list_b_typename,"GISTIC_1",sep="_"),paste(list_b_typename,"GISTIC_2",sep="_"),"p.value","adj.p")


for ( i in 1:length(commonPeak) ){
  
  fisher.res[commonPeak[i],c(1:6)] = c(mat_combine_xtabs[list_c_typename,,commonPeak[i]],mat_combine_xtabs[list_b_typename,,commonPeak[i]])
  fisher.res[commonPeak[i],"p.value"] = fisher.test(mat_combine_xtabs[,,commonPeak[i]])[["p.value"]]
  if (i/10 == round(i/10)) {print(i)}
}

fisher.res$adj.p <- p.adjust(fisher.res$p.value ,method = "fdr")

# write.csv(fisher.res, file = paste("CNA_comparation_peaks_fisher.exact_",list_c_typename,"_vs_",list_b_typename,".csv",sep = ""))

CNV_FDR_criteria<-0.05
names<-rownames(fisher.res)
pdf(paste("CNA_comparation_peaks_histplot_",list_c_typename,"_vs_",list_b_typename,"_FDR_",CNV_FDR_criteria,".pdf",sep = ""),height =  4, width =  10)
barplot(-log10(fisher.res$adj.p),col = "#EAAC53",border = NA, axisnames = F,lwd  = 2,ylim = c(-4,4))
par(new= T)
barplot(-log10(1/fisher.res$adj.p),col = "#B3D486",border = NA, lwd  = 2,ylim = c(-4,4),names.arg=as.character(as.data.frame(strsplit(names,"_"))[3,]),las=3,cex.names=0.5)

abline(h = -log10(CNV_FDR_criteria), lty = 1,lwd  = 1) 
abline(h = -log10(1/CNV_FDR_criteria), lty = 1,lwd  = 1) 
par(new= F)
dev.off()



