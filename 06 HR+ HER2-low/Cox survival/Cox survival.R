table(FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE)
# HER2_0      HER2_low HER2_positive 
# 91           434           182 

##single var

for_three <- data.frame(Del_Peak.42_17q11.2=as.character(FUSCC_HER2_low_project_GISTICpeaks.del.thre["Del_Peak.42_17q11.2",]),
                        Del_Peak.43_17q21.31=as.character(FUSCC_HER2_low_project_GISTICpeaks.del.thre["Del_Peak.43_17q21.31",]),
                        Amp_Peak.21_17q12=as.character(FUSCC_HER2_low_project_GISTICpeaks.amp.thre["Amp_Peak.21_17q12",]))
rownames(for_three) <- colnames(FUSCC_HER2_low_project_GISTICpeaks.del.thre)
for_three <- for_three[FUSCC_HER2_low_project_Cohort.Info[rownames(for_three),"HR_status"] %in% c("Positive") & FUSCC_HER2_low_project_Cohort.Info[rownames(for_three),"HER2_low_status_RE"] %in% c("HER2_low","HER2_0") ,]
for_three[,c("HER2_low_status_RE","DMFS_status","DMFS_time_Months")] <- FUSCC_HER2_low_project_Cohort.Info[rownames(for_three),c("HER2_low_status_RE","DMFS_status","DMFS_time_Months")]
for_three$Del_Peak.42_17q11.2 <- gsub("0","Others",for_three$Del_Peak.42_17q11.2)
for_three$Del_Peak.42_17q11.2 <- gsub("1|2","Ld",for_three$Del_Peak.42_17q11.2)
for_three$Del_Peak.42_17q11.2 <- factor(for_three$Del_Peak.42_17q11.2,levels = c("Others","Ld"))
for_three$Del_Peak.43_17q21.31 <- gsub("0","Others",for_three$Del_Peak.43_17q21.31)
for_three$Del_Peak.43_17q21.31 <- gsub("1|2","Ld",for_three$Del_Peak.43_17q21.31)
for_three$Del_Peak.43_17q21.31 <- factor(for_three$Del_Peak.43_17q21.31,levels = c("Others","Ld"))
for_three$Amp_Peak.21_17q12 <- gsub("0","Others",for_three$Amp_Peak.21_17q12)
for_three$Amp_Peak.21_17q12 <- gsub("1|2","Ga",for_three$Amp_Peak.21_17q12)
for_three$Amp_Peak.21_17q12 <- factor(for_three$Amp_Peak.21_17q12,levels = c("Others","Ga"))
phe <- for_three

library("survival")
library("survminer")


pdf("02_survival_HR+_cox_Del_Peak.42_17q11.2.pdf",width = 5.5,height = 2.5)
cox_model <- NA;cox_model = coxph(Surv  (DMFS_time_Months, DMFS_status) ~ Del_Peak.42_17q11.2, data = phe)
ggforest(cox_model)
cox_model <- NA;cox_model = coxph(Surv  (DMFS_time_Months, DMFS_status) ~ HER2_low_status_RE + Del_Peak.42_17q11.2 , data = phe)
ggforest(cox_model)
dev.off()

pdf("02_survival_HR+_cox_Del_Peak.43_17q21.31.pdf",width = 5.5,height = 2.5)
cox_model <- NA;cox_model = coxph(Surv  (DMFS_time_Months, DMFS_status) ~ Del_Peak.43_17q21.31, data = phe)
ggforest(cox_model)
cox_model <- NA;cox_model = coxph(Surv  (DMFS_time_Months, DMFS_status) ~ HER2_low_status_RE + Del_Peak.43_17q21.31 , data = phe)
ggforest(cox_model)
dev.off()

pdf("02_survival_HR+_cox_Amp_Peak.21_17q12.pdf",width = 5.5,height = 2.5)
cox_model <- NA;cox_model = coxph(Surv  (DMFS_time_Months, DMFS_status) ~ Amp_Peak.21_17q12, data = phe)
ggforest(cox_model)
cox_model <- NA;cox_model = coxph(Surv  (DMFS_time_Months, DMFS_status) ~ HER2_low_status_RE +Amp_Peak.21_17q12, data = phe)
ggforest(cox_model)
dev.off()
