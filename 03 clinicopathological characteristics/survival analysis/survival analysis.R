table(FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE)
# HER2_0      HER2_low HER2_positive 
# 91           434           182 

##follow-up time
prognosis_time <- FUSCC_HER2_low_project_Cohort.Info[,"OS_time_months"]/12
median(prognosis_time,na.rm = T)
# 6.899384
summary(prognosis_time)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.01917 5.89733 6.89938 6.26582 7.45791 9.00000       2 

##follow-up time
prognosis_time <- FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_low"),"OS_time_months"]/12
median(prognosis_time,na.rm = T)
# 6.90486
summary(prognosis_time)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.01917 5.96578 6.90486 6.37969 7.54004 9.00000       1

#HR_all HER2-0 VS HER2-LOW
her2_0_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_0")& FUSCC_HER2_low_project_Cohort.Info$HR_status %in% c("Negative","Positive"),"PatientCode"]
her2_low_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_low")& FUSCC_HER2_low_project_Cohort.Info$HR_status %in% c("Negative","Positive"),"PatientCode"]

group_list<-c(rep("HER2_0",length(her2_0_list)),rep("HER2_low",length(her2_low_list)))
names(group_list)<-c(her2_0_list,her2_low_list)
table(group_list)

phe<-FUSCC_HER2_low_project_Cohort.Info[names(group_list),c("OS_status","OS_time_months")]
phe$cluster<-(group_list)
head(phe)

colnames(phe)=c('event','time',"cluster")

library(survival)
library(survminer)

pdf("survival_HR_all_0_vs_low_OS(risk table).pdf",width = 5,height = 5)
sfit <- survfit(Surv(time, event)~cluster, data=phe)
sfit
summary(sfit)
ggsurvplot(sfit,palette = c("jco"),risk.table =TRUE,pval =TRUE,break.time.by = 12,xlab ="Time(months)",risk.table.y.text = FALSE,xlim = c(0,108),risk.table.height = 0.3)
dev.off()

surv_diff <- survdiff(Surv(time, event)~cluster, data=phe)
p.value <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) -1)
p.value

#HR_all HER2-0 VS HER2-LOW
her2_0_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_0")& FUSCC_HER2_low_project_Cohort.Info$HR_status %in% c("Negative","Positive"),"PatientCode"]
her2_low_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_low")& FUSCC_HER2_low_project_Cohort.Info$HR_status %in% c("Negative","Positive"),"PatientCode"]

group_list<-c(rep("HER2_0",length(her2_0_list)),rep("HER2_low",length(her2_low_list)))
names(group_list)<-c(her2_0_list,her2_low_list)
table(group_list)

phe<-FUSCC_HER2_low_project_Cohort.Info[names(group_list),c("DMFS_status","DMFS_time_Months")]
phe$cluster<-(group_list)
head(phe)

colnames(phe)=c('event','time',"cluster")

library(survival)
library(survminer)

pdf("survival_HR_all_0_vs_low_DMFS(risk table).pdf",width = 5,height = 5)
sfit <- survfit(Surv(time, event)~cluster, data=phe)
sfit
summary(sfit)
ggsurvplot(sfit,palette = c("jco"),risk.table =TRUE,pval =TRUE,break.time.by = 12,xlab ="Time(months)",risk.table.y.text = FALSE,xlim = c(0,108),risk.table.height = 0.3)
dev.off()

surv_diff <- survdiff(Surv(time, event)~cluster, data=phe)
p.value <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) -1)
p.value

#HR+ HER2-0 VS HER2-LOW
her2_0_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_0")& FUSCC_HER2_low_project_Cohort.Info$HR_status %in% c("Positive"),"PatientCode"]
her2_low_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_low")& FUSCC_HER2_low_project_Cohort.Info$HR_status %in% c("Positive"),"PatientCode"]

group_list<-c(rep("HER2_0",length(her2_0_list)),rep("HER2_low",length(her2_low_list)))
names(group_list)<-c(her2_0_list,her2_low_list)
table(group_list)

phe<-FUSCC_HER2_low_project_Cohort.Info[names(group_list),c("OS_status","OS_time_months")]
phe$cluster<-(group_list)
head(phe)

colnames(phe)=c('event','time',"cluster")

library(survival)
library(survminer)

pdf("survival_HR+_0_vs_low_OS(risk table).pdf",width = 5,height = 5)
sfit <- survfit(Surv(time, event)~cluster, data=phe)
sfit
summary(sfit)
ggsurvplot(sfit,palette = c("jco"),risk.table =TRUE,pval =TRUE,break.time.by = 12,xlab ="Time(months)",risk.table.y.text = FALSE,xlim = c(0,108),risk.table.height = 0.3)
dev.off()

surv_diff <- survdiff(Surv(time, event)~cluster, data=phe)
p.value <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) -1)
p.value

#HR+ HER2-0 VS HER2-LOW
her2_0_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_0")& FUSCC_HER2_low_project_Cohort.Info$HR_status %in% c("Positive"),"PatientCode"]
her2_low_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_low")& FUSCC_HER2_low_project_Cohort.Info$HR_status %in% c("Positive"),"PatientCode"]

group_list<-c(rep("HER2_0",length(her2_0_list)),rep("HER2_low",length(her2_low_list)))
names(group_list)<-c(her2_0_list,her2_low_list)
table(group_list)

phe<-FUSCC_HER2_low_project_Cohort.Info[names(group_list),c("DMFS_status","DMFS_time_Months")]
phe$cluster<-(group_list)
head(phe)

colnames(phe)=c('event','time',"cluster")

library(survival)
library(survminer)

pdf("survival_HR+_0_vs_low_DMFS(risk table).pdf",width = 5,height = 5)
sfit <- survfit(Surv(time, event)~cluster, data=phe)
sfit
summary(sfit)
ggsurvplot(sfit,palette = c("jco"),risk.table =TRUE,pval =TRUE,break.time.by = 12,xlab ="Time(months)",risk.table.y.text = FALSE,xlim = c(0,108),risk.table.height = 0.3)
dev.off()

surv_diff <- survdiff(Surv(time, event)~cluster, data=phe)
p.value <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) -1)
p.value

#HR- HER2-0 VS HER2-LOW
her2_0_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_0")& FUSCC_HER2_low_project_Cohort.Info$HR_status %in% c("Negative"),"PatientCode"]
her2_low_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_low")& FUSCC_HER2_low_project_Cohort.Info$HR_status %in% c("Negative"),"PatientCode"]

group_list<-c(rep("HER2_0",length(her2_0_list)),rep("HER2_low",length(her2_low_list)))
names(group_list)<-c(her2_0_list,her2_low_list)
table(group_list)

phe<-FUSCC_HER2_low_project_Cohort.Info[names(group_list),c("OS_status","OS_time_months")]
phe$cluster<-(group_list)
head(phe)

colnames(phe)=c('event','time',"cluster")

library(survival)
library(survminer)

pdf("survival_HR-_0_vs_low_OS(risk table).pdf",width = 5,height = 5)
sfit <- survfit(Surv(time, event)~cluster, data=phe)
sfit
summary(sfit)
ggsurvplot(sfit,palette = c("jco"),risk.table =TRUE,pval =TRUE,break.time.by = 12,xlab ="Time(months)",risk.table.y.text = FALSE,xlim = c(0,108),risk.table.height = 0.3)
dev.off()

surv_diff <- survdiff(Surv(time, event)~cluster, data=phe)
p.value <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) -1)
p.value

#HR- HER2-0 VS HER2-LOW
her2_0_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_0")& FUSCC_HER2_low_project_Cohort.Info$HR_status %in% c("Negative"),"PatientCode"]
her2_low_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_low")& FUSCC_HER2_low_project_Cohort.Info$HR_status %in% c("Negative"),"PatientCode"]

group_list<-c(rep("HER2_0",length(her2_0_list)),rep("HER2_low",length(her2_low_list)))
names(group_list)<-c(her2_0_list,her2_low_list)
table(group_list)

phe<-FUSCC_HER2_low_project_Cohort.Info[names(group_list),c("DMFS_status","DMFS_time_Months")]
phe$cluster<-(group_list)
head(phe)

colnames(phe)=c('event','time',"cluster")

library(survival)
library(survminer)

pdf("survival_HR-_0_vs_low_DMFS(risk table).pdf",width = 5,height = 5)
sfit <- survfit(Surv(time, event)~cluster, data=phe)
sfit
summary(sfit)
ggsurvplot(sfit,palette = c("jco"),risk.table =TRUE,pval =TRUE,break.time.by = 12,xlab ="Time(months)",risk.table.y.text = FALSE,xlim = c(0,108),risk.table.height = 0.3)
dev.off()

surv_diff <- survdiff(Surv(time, event)~cluster, data=phe)
p.value <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) -1)
p.value



