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

##ERBB2 CNV
Subtype = t(FUSCC_HER2_low_project_GISTICgene.thre)[names(group_list[names(group_list) %in% colnames(FUSCC_HER2_low_project_GISTICgene.thre)]),] ; Subtype<-as.data.frame(Subtype);Subtype$cluster = as.character(group_list[names(group_list) %in% colnames(FUSCC_HER2_low_project_GISTICgene.thre)])
ERBB2_Table <- xtabs(~ ERBB2 + cluster, data = Subtype)

Color_ERBB2<-c("#C2C2C4","#EC7878","#E21A21","#88A3D8","#3653A5")
names(Color_ERBB2)<-c("0","1","2","-1","-2")

pdf("ERBB2_amp_percentage_HR_neg_HER2_low_basal_vs_nonbasal.pdf",width = 5,height = 5)
barplot(apply(ERBB2_Table, 2, function(x){x/sum(x)}), col = Color_ERBB2[row.names(ERBB2_Table)] ,horiz =T,border = F)
dev.off()

fisher.test(cbind(ERBB2_Table[,3],ERBB2_Table[,1]))
# p-value = 0.245
fisher.test(cbind(ERBB2_Table[,3],ERBB2_Table[,2]))
# p-value = 0.9021
fisher.test(cbind(ERBB2_Table[,3],ERBB2_Table[,4]))
# p-value = 9.064e-07


##PTK6 CNV
Subtype = t(FUSCC_HER2_low_project_GISTICgene.thre)[names(group_list[names(group_list) %in% colnames(FUSCC_HER2_low_project_GISTICgene.thre)]),] ; Subtype<-as.data.frame(Subtype);Subtype$cluster = as.character(group_list[names(group_list) %in% colnames(FUSCC_HER2_low_project_GISTICgene.thre)])
PTK6_Table <- xtabs(~ PTK6 + cluster, data = Subtype)

Color_PTK6<-c("#C2C2C4","#EC7878","#E21A21","#88A3D8","#3653A5")
names(Color_PTK6)<-c("0","1","2","-1","-2")

pdf("PTK6_amp_percentage_HR_neg_HER2_low_basal_vs_nonbasal.pdf",width = 5,height = 5)
barplot(apply(PTK6_Table, 2, function(x){x/sum(x)}), col = Color_PTK6[row.names(PTK6_Table)] ,horiz =T,border = F)
dev.off()

fisher.test(cbind(PTK6_Table[,3],PTK6_Table[,1]))
# p-value = 0.1817
fisher.test(cbind(PTK6_Table[,3],PTK6_Table[,2]))
# p-value = 0.2129
fisher.test(cbind(PTK6_Table[,3],PTK6_Table[,4]))
# p-value = 0.4522

##FGFR4 CNV
Subtype = t(FUSCC_HER2_low_project_GISTICgene.thre)[names(group_list[names(group_list) %in% colnames(FUSCC_HER2_low_project_GISTICgene.thre)]),] ; Subtype<-as.data.frame(Subtype);Subtype$cluster = as.character(group_list[names(group_list) %in% colnames(FUSCC_HER2_low_project_GISTICgene.thre)])
FGFR4_Table <- xtabs(~ FGFR4 + cluster, data = Subtype)

Color_FGFR4<-c("#C2C2C4","#EC7878","#E21A21","#88A3D8","#3653A5")
names(Color_FGFR4)<-c("0","1","2","-1","-2")

pdf("FGFR4_amp_percentage_HR_neg_HER2_low_basal_vs_nonbasal.pdf",width = 5,height = 5)
barplot(apply(FGFR4_Table, 2, function(x){x/sum(x)}), col = Color_FGFR4[row.names(FGFR4_Table)] ,horiz =T,border = F)
dev.off()

fisher.test(cbind(FGFR4_Table[,3],FGFR4_Table[,1]))
# p-value = 0.1615
fisher.test(cbind(FGFR4_Table[,3],FGFR4_Table[,2]))
# p-value = 0.2873
fisher.test(cbind(FGFR4_Table[,3],FGFR4_Table[,4]))
# p-value = 0.825


##ERBB4 CNV
Subtype = t(FUSCC_HER2_low_project_GISTICgene.thre)[names(group_list[names(group_list) %in% colnames(FUSCC_HER2_low_project_GISTICgene.thre)]),] ; Subtype<-as.data.frame(Subtype);Subtype$cluster = as.character(group_list[names(group_list) %in% colnames(FUSCC_HER2_low_project_GISTICgene.thre)])
ERBB4_Table <- xtabs(~ ERBB4 + cluster, data = Subtype)

Color_ERBB4<-c("#C2C2C4","#EC7878","#E21A21","#88A3D8","#3653A5")
names(Color_ERBB4)<-c("0","1","2","-1","-2")

pdf("ERBB4_amp_percentage_HR_neg_HER2_low_basal_vs_nonbasal.pdf",width = 5,height = 5)
barplot(apply(ERBB4_Table, 2, function(x){x/sum(x)}), col = Color_ERBB4[row.names(ERBB4_Table)] ,horiz =T,border = F)
dev.off()

fisher.test(cbind(ERBB4_Table[,3],ERBB4_Table[,1]))
# p-value = 1
fisher.test(cbind(ERBB4_Table[,3],ERBB4_Table[,2]))
# p-value = 0.5106
fisher.test(cbind(ERBB4_Table[,3],ERBB4_Table[,4]))
# p-value = 0.448
