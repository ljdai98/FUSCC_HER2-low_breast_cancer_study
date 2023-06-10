table(FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE)
# HER2_0      HER2_low HER2_positive 
# 91           434           182 


FUSCC_HER2_low_project_Cohort.Info$PAM50 <- factor(FUSCC_HER2_low_project_Cohort.Info$PAM50,levels = rev(c("LumA","LumB","Her2","Basal","Normal")))


##HER2-low vs HER2-0 all patients
her2_0_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_0"),"PatientCode"]
her2_low_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_low"),"PatientCode"]

group_list<-c(rep("HER2_0",length(her2_0_list)),rep("HER2_low",length(her2_low_list)))
names(group_list)<-c(her2_0_list,her2_low_list)
table(group_list)
# group_list
# HER2_0 HER2_low 
# 91      434

col <- c("#000286","#87CCFA","#FF67B8","#FE0302","#36E200")
names(col) <- c("LumA","LumB","Her2","Basal","Normal")

table_temp <- table(FUSCC_HER2_low_project_Cohort.Info[names(group_list),c("PAM50","HER2_low_status_RE")])
table_temp
# PAM50    HER2_0 HER2_low
# Normal      4       27
# Basal      31       58
# Her2        4       36
# LumB       24      141
# LumA       25      159

fisher.test(table_temp)
# p-value = 0.000354

table_temp <- as.data.frame.matrix(prop.table(table_temp,2))
table_temp <- table_temp[,c(1,2)]
round(table_temp*100,1)
# HER2_0 HER2_low
# Normal    4.5      6.4
# Basal    35.2     13.8
# Her2      4.5      8.6
# LumB     27.3     33.5
# LumA     28.4     37.8

pdf("bar chart HER2-low vs HER2-0 all patients.pdf",width = 3,height = 5)
barplot(as.matrix(table_temp),width=2,space = 1,col = rev(col),border = NA,ylab = "Percentage",xlab = "All patients")
dev.off()



##HER2-low vs HER2-0 HR pos
her2_0_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_0") & FUSCC_HER2_low_project_Cohort.Info$HR_status %in% "Positive","PatientCode"]
her2_low_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_low")& FUSCC_HER2_low_project_Cohort.Info$HR_status %in% "Positive","PatientCode"]

group_list<-c(rep("HER2_0",length(her2_0_list)),rep("HER2_low",length(her2_low_list)))
names(group_list)<-c(her2_0_list,her2_low_list)
table(group_list)
# group_list
# HER2_0 HER2_low 
# 63      361

col <- c("#000286","#87CCFA","#FF67B8","#FE0302","#36E200")
names(col) <- c("LumA","LumB","Her2","Basal","Normal")

table_temp <- table(FUSCC_HER2_low_project_Cohort.Info[names(group_list),c("PAM50","HER2_low_status_RE")])
table_temp
# PAM50    HER2_0 HER2_low
# Normal      4       24
# Basal       5       12
# Her2        4       22
# LumB       23      140
# LumA       25      157

fisher.test(table_temp)
# p-value = 0.5336

table_temp <- as.data.frame.matrix(prop.table(table_temp,2))
table_temp <- table_temp[,c(1,2)]
round(table_temp*100,1)
# HER2_0 HER2_low
# Normal    6.6      6.8
# Basal     8.2      3.4
# Her2      6.6      6.2
# LumB     37.7     39.4
# LumA     41.0     44.2

pdf("bar chart HER2-low vs HER2-0 HR pos.pdf",width = 3,height = 5)
barplot(as.matrix(table_temp),width=2,space = 1,col = rev(col),border = NA,ylab = "Percentage",xlab = "HR−positive")
dev.off()


##HER2-low vs HER2-0 HR neg
her2_0_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_0") & FUSCC_HER2_low_project_Cohort.Info$HR_status %in% "Negative","PatientCode"]
her2_low_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_low")& FUSCC_HER2_low_project_Cohort.Info$HR_status %in% "Negative","PatientCode"]

group_list<-c(rep("HER2_0",length(her2_0_list)),rep("HER2_low",length(her2_low_list)))
names(group_list)<-c(her2_0_list,her2_low_list)
table(group_list)
# group_list
# HER2_0 HER2_low 
# 28       73

col <- c("#000286","#87CCFA","#FF67B8","#FE0302","#36E200")
names(col) <- c("LumA","LumB","Her2","Basal","Normal")

table_temp <- table(FUSCC_HER2_low_project_Cohort.Info[names(group_list),c("PAM50","HER2_low_status_RE")])
table_temp
# HER2_low_status_RE
# PAM50    HER2_0 HER2_low
# Normal      0        3
# Basal      26       46
# Her2        0       14
# LumB        1        1
# LumA        0        2

fisher.test(table_temp)
#  0.01454

table_temp <- as.data.frame.matrix(prop.table(table_temp,2))
table_temp <- table_temp[,c(1,2)]
round(table_temp*100,1)
# HER2_0 HER2_low
# Normal    0.0      4.5
# Basal    96.3     69.7
# Her2      0.0     21.2
# LumB      3.7      1.5
# LumA      0.0      3.0

pdf("bar chart HER2-low vs HER2-0 HR neg.pdf",width = 3,height = 5)
barplot(as.matrix(table_temp),width=2,space = 1,col = rev(col),border = NA,ylab = "Percentage",xlab = "HR−negative")
dev.off()





