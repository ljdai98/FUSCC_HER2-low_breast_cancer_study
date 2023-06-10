table(FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE)
# HER2_0      HER2_low HER2_positive 
# 91           434           182 

#HER2-low vs HER2-0

her2_0_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_0") & FUSCC_HER2_low_project_Cohort.Info$HR_status %in% "Negative","PatientCode"]
her2_low_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_low")& FUSCC_HER2_low_project_Cohort.Info$HR_status %in% "Negative","PatientCode"]

group_list<-c(rep("HER2_0",length(her2_0_list)),rep("HER2_low",length(her2_low_list)))
names(group_list)<-c(her2_0_list,her2_low_list)
table(group_list)
# group_list
# HER2_0 HER2_low 
# 28       73


FUSCC_HER2_low_project_Cohort.Info$PAM50 <- factor(FUSCC_HER2_low_project_Cohort.Info$PAM50,levels = c("Normal","LumA","LumB","Her2","Basal"))



## FUSCC HER2-low vs HER2-0
col <- c("#73AAF9","#73AAF9","#73AAF9","#73AAF9","#D9D9D9")
names(col) <- c("Normal","LumA","LumB","Her2","Basal")

table_temp <- table(FUSCC_HER2_low_project_Cohort.Info[names(group_list),c("PAM50","HER2_low_status_RE")])
table_temp
# PAM50    HER2_0 HER2_low
# Normal      0        3
# LumA        0        2
# LumB        1        1
# Her2        0       14
# Basal      26       46

fisher.test(cbind(c(sum(table_temp[,1])-table_temp["Basal",1],table_temp["Basal",1]),c(sum(table_temp[,2])-table_temp["Basal",2],table_temp["Basal",2])))
# p-value = 0.005309

table_temp <- as.data.frame.matrix(prop.table(table_temp,2))
table_temp <- table_temp[,c(2,1)]
round(table_temp*100,1)
#        HER2_low HER2_0
# Normal      4.5    0.0
# LumA        3.0    0.0
# LumB        1.5    3.7
# Her2       21.2    0.0
# Basal      69.7   96.3
pdf("bar chart FUSCC HER2-low vs HER2-0.pdf",width = 3,height = 5)
barplot(as.matrix(table_temp),width=2,space = 1,col = col,border = NA,ylab = "Percentage",xlab = "FUSCC")
dev.off()

## NPJ-2021 HER2-low vs HER2-0
col <- c("#73AAF9","#73AAF9","#73AAF9","#73AAF9","#D9D9D9")
names(col) <- c("Normal","LumA","LumB","Her2","Basal")

table_temp <- cbind(c(12,5,1,28,265),c(10,2,0,9,105))
rownames(table_temp) <- c("Normal","LumA","LumB","Her2","Basal")
colnames(table_temp) <- c("HER2_0","HER2_low")
table_temp
#        HER2_0 HER2_low
# Normal     12       10
# LumA        5        2
# LumB        1        0
# Her2       28        9
# Basal     265      105

fisher.test(cbind(c(sum(table_temp[,1])-table_temp["Basal",1],table_temp["Basal",1]),c(sum(table_temp[,2])-table_temp["Basal",2],table_temp["Basal",2])))
# p-value = 0.6607

table_temp <- as.data.frame.matrix(prop.table(table_temp,2))
table_temp <- table_temp[,c(2,1)]
round(table_temp*100,1)
#        HER2_low HER2_0
# Normal      7.9    3.9
# LumA        1.6    1.6
# LumB        0.0    0.3
# Her2        7.1    9.0
# Basal      83.3   85.2
pdf("bar chart npj2021 HER2-low vs HER2-0.pdf",width = 3,height = 5)
barplot(as.matrix(table_temp),width=2,space = 1,col = col,border = NA,ylab = "Percentage",xlab = "NPJ-2021")
dev.off()

## FUSCC HER2-low vs NPJ-2021 HER2-low
table_temp <- cbind(c(3,2,1,14,46),c(10,2,0,9,105))
rownames(table_temp) <- c("Normal","LumA","LumB","Her2","Basal")
colnames(table_temp) <- c("FUSCC","NPJ-2021")
table_temp
 
#        FUSCC NPJ-2021
# Normal     3       10
# LumA       2        2
# LumB       1        0
# Her2      14        9
# Basal     46      105
fisher.test(cbind(c(sum(table_temp[,1])-table_temp["Basal",1],table_temp["Basal",1]),c(sum(table_temp[,2])-table_temp["Basal",2],table_temp["Basal",2])))
# p-value = 0.04057


#HER2-1 vs HER2-2

her2_1_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_long_RE %in% c("HER2_1") & FUSCC_HER2_low_project_Cohort.Info$HR_status %in% "Negative","PatientCode"]
her2_2_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_long_RE %in% c("HER2_2")& FUSCC_HER2_low_project_Cohort.Info$HR_status %in% "Negative","PatientCode"]

group_list<-c(rep("HER2_1",length(her2_1_list)),rep("HER2_2",length(her2_2_list)))
names(group_list)<-c(her2_1_list,her2_2_list)
table(group_list)
# group_list
# HER2_1 HER2_2 
# 48     25


FUSCC_HER2_low_project_Cohort.Info$PAM50 <- factor(FUSCC_HER2_low_project_Cohort.Info$PAM50,levels = c("Normal","LumA","LumB","Her2","Basal"))

## FUSCC HER2-2 vs HER2-1
col <- c("#73AAF9","#73AAF9","#73AAF9","#73AAF9","#D9D9D9")
names(col) <- c("Normal","LumA","LumB","Her2","Basal")

table_temp <- table(FUSCC_HER2_low_project_Cohort.Info[names(group_list),c("PAM50","HER2_low_status_long_RE")])
table_temp
# PAM50    HER2_1 HER2_2
# Normal      3      0
# LumA        2      0
# LumB        0      1
# Her2        3     11
# Basal      37      9

fisher.test(cbind(c(sum(table_temp[,1])-table_temp["Basal",1],table_temp["Basal",1]),c(sum(table_temp[,2])-table_temp["Basal",2],table_temp["Basal",2])))
# p-value = 0.003148

table_temp <- as.data.frame.matrix(prop.table(table_temp,2))
round(table_temp*100,1)
#        HER2_1 HER2_2
# Normal    6.7    0.0
# LumA      4.4    0.0
# LumB      0.0    4.8
# Her2      6.7   52.4
# Basal    82.2   42.9
pdf("bar chart FUSCC HER2-2 vs HER2-1.pdf",width = 3,height = 5)
barplot(as.matrix(table_temp),width=2,space = 1,col = col,border = NA,ylab = "Percentage",xlab = "FUSCC")
dev.off()

## NPJ-2021 HER2-2 vs HER2-1
col <- c("#73AAF9","#73AAF9","#73AAF9","#73AAF9","#D9D9D9")
names(col) <- c("Normal","LumA","LumB","Her2","Basal")

table_temp <- cbind(c(6,0,0,7,76),c(4,2,0,2,29))
rownames(table_temp) <- c("Normal","LumA","LumB","Her2","Basal")
colnames(table_temp) <- c("HER2_1","HER2_2")
table_temp
#        HER2_1 HER2_2
# Normal      6      4
# LumA        0      2
# LumB        0      0
# Her2        7      2
# Basal      76     29

fisher.test(cbind(c(sum(table_temp[,1])-table_temp["Basal",1],table_temp["Basal",1]),c(sum(table_temp[,2])-table_temp["Basal",2],table_temp["Basal",2])))
# p-value = 0.4311

table_temp <- as.data.frame.matrix(prop.table(table_temp,2))
round(table_temp*100,1)
#        HER2_1 HER2_2
# Normal    6.7   10.8
# LumA      0.0    5.4
# LumB      0.0    0.0
# Her2      7.9    5.4
# Basal    85.4   78.4
pdf("bar chart npj2021 HER2-2 vs HER2-1.pdf",width = 3,height = 5)
barplot(as.matrix(table_temp),width=2,space = 1,col = col,border = NA,ylab = "Percentage",xlab = "NPJ-2021")
dev.off()

## FUSCC HER2-2 vs NPJ-2021 HER2-2
table_temp <- cbind(c(0,0,1,11,9),c(4,2,0,2,29))
rownames(table_temp) <- c("Normal","LumA","LumB","Her2","Basal")
colnames(table_temp) <- c("FUSCC","NPJ-2021")
table_temp

#        FUSCC NPJ-2021
# Normal     0        4
# LumA       0        2
# LumB       1        0
# Her2      11        2
# Basal      9       29
fisher.test(cbind(c(sum(table_temp[,1])-table_temp["Basal",1],table_temp["Basal",1]),c(sum(table_temp[,2])-table_temp["Basal",2],table_temp["Basal",2])))
# p-value =0.009656


