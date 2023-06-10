table(FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE)
# HER2_0      HER2_low HER2_positive 
# 91           434           182 

##pheno matrix

phe=FUSCC_HER2_low_project_Cohort.Info[,c("Age_at_surgery","Menopause","Laterality",
                                                           "HER2_low_status_long_RE","Histology","Grade","Ki67_percentage","HR_status","sTILs","iTILs",
                                                           "Stage_pT","Stage_pN",
                                                           "TMB","HRD","ACF","PAM50","FUSCCTNBC_mRNA_Subtype",
                                                           "Surgery_breast","Surgery_axilla","Adjuvant_chemotherapy","Taxane_usage","Anthracyclines_usage","Platinum_usage","Capecitabine_usage","Adjuvant_radiotherapy","Adjuvant_endocrine_therapy")]

phe$Age_cla <- as.numeric(phe$Age_at_surgery)
phe[phe$Age_at_surgery >=40 & phe$Age_at_surgery <60, "Age_cla"] <- c("40-59")
phe[phe$Age_at_surgery >=60, "Age_cla"] <- c("≥60")
phe[phe$Age_at_surgery < 40,"Age_cla" ] <- c("<40")
phe$Age_cla <- factor(phe$Age_cla,levels = c("<40","40-59","≥60"))
table(phe$Age_cla)

phe$Ki67_percentage_cla   <- as.numeric(phe$Ki67_percentage)
phe[phe$Ki67_percentage >=20, "Ki67_percentage_cla"] <- c("≥20")
phe[phe$Ki67_percentage <20,"Ki67_percentage_cla" ] <- c("<20")
phe$Ki67_percentage_cla <- factor(phe$Ki67_percentage_cla,levels = c("<20","≥20"))
table(phe$Ki67_percentage_cla)

phe[phe$ACF == max(phe$ACF,na.rm = T) & !is.na(phe$ACF == max(phe$ACF,na.rm = T)),"HRD"] <- NA

## HER2_0 VS HER2_LOW VS HER2_POS
her2_0_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_0"),"PatientCode"]
her2_low_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_low"),"PatientCode"]
her2_pos_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_positive"),"PatientCode"]

group_list<-c(rep("HER2_0",length(her2_0_list)),rep("HER2_low",length(her2_low_list)),rep("HER2_positive",length(her2_pos_list)))
names(group_list)<-c(her2_0_list,her2_low_list,her2_pos_list)
table(group_list)
# group_list
# HER2_0      HER2_low HER2_positive 
# 91           434           182 

phe2 <-phe[names(group_list),c("Age_at_surgery","Age_cla","Menopause","Laterality",
              "HER2_low_status_long_RE","Histology","Grade","Ki67_percentage_cla","HR_status","sTILs","iTILs",
              "Stage_pT","Stage_pN",
              "TMB","HRD","PAM50")]


phe2$cluster <- group_list

vars<-colnames(phe2)[!colnames(phe2) %in% c("cluster")]

library("tableone")
tableOne<-CreateTableOne(vars=vars,strata = c("cluster"),data=phe2 ,addOverall = T,includeNA = F)

tab2<-print(tableOne,showAllLevels = T,exact =c(colnames(phe)), quote = FALSE, noSpaces = TRUE)
tab3<-as.matrix(tab2)

tableOne<-CreateTableOne(vars=vars,strata = c("cluster"),data=phe2 ,addOverall = T,includeNA = T)

tab2<-print(tableOne,showAllLevels = T,exact =c(colnames(phe)), quote = FALSE, noSpaces = TRUE)
tab3_2<-as.matrix(tab2)

tab3_combined <- tab3_2
tab3_combined[!is.na(tab3_combined[,"level"] ),] <- tab3

tab3_combined[is.na(tab3_combined[,"level"] ),c("Overall",unique(group_list))] <- stringr::str_match(string = tab3_combined[is.na(tab3_combined[,"level"] ),c("Overall",unique(group_list))],
                                                                                                     pattern = ".*(?= )")

write.csv(tab3_combined,"clinical_features_her2_0_vs_her2_low_vs_her2_pos.csv",quote=F)


## HER2_0 VS HER2_LOW
group_list<-c(rep("HER2_0",length(her2_0_list)),rep("HER2_low",length(her2_low_list)))
names(group_list)<-c(her2_0_list,her2_low_list)
table(group_list)
# HER2_0      HER2_low 
#     91      434      

phe2 <-phe[names(group_list),c("Age_at_surgery","Age_cla","Menopause","Laterality",
                               "HER2_low_status_long_RE","Histology","Grade","Ki67_percentage_cla","HR_status","sTILs","iTILs",
                               "Stage_pT","Stage_pN",
                               "TMB","HRD","PAM50")]


phe2$cluster <- group_list

vars<-colnames(phe2)[!colnames(phe2) %in% c("cluster")]

library("tableone")

tableOne<-CreateTableOne(vars=vars,strata = c("cluster"),data=phe2 ,addOverall = T,includeNA = F)

tab2<-print(tableOne,showAllLevels = T,exact =c(colnames(phe)), quote = FALSE, noSpaces = TRUE)
tab3<-as.matrix(tab2)

tableOne<-CreateTableOne(vars=vars,strata = c("cluster"),data=phe2 ,addOverall = T,includeNA = T)

tab2<-print(tableOne,showAllLevels = T,exact =c(colnames(phe)), quote = FALSE, noSpaces = TRUE)
tab3_2<-as.matrix(tab2)

tab3_combined <- tab3_2
tab3_combined[!is.na(tab3_combined[,"level"] ),] <- tab3

tab3_combined[is.na(tab3_combined[,"level"] ),c("Overall",unique(group_list))] <- stringr::str_match(string = tab3_combined[is.na(tab3_combined[,"level"] ),c("Overall",unique(group_list))],
                                                                                                     pattern = ".*(?= )")
write.csv(tab3_combined,"clinical_features_her2_0_vs_her2_low.csv",quote=F)



## HER2_LOW VS HER2_POS
group_list<-c(rep("HER2_low",length(her2_low_list)),rep("HER2_positive",length(her2_pos_list)))
names(group_list)<-c(her2_low_list,her2_pos_list)
table(group_list)
# HER2_low HER2_positive 
# 434           182

phe2 <-phe[names(group_list),c("Age_at_surgery","Age_cla","Menopause","Laterality",
                               "HER2_low_status_long_RE","Histology","Grade","Ki67_percentage_cla","HR_status","sTILs","iTILs",
                               "Stage_pT","Stage_pN",
                               "TMB","HRD","PAM50")]


phe2$cluster <- group_list

vars<-colnames(phe2)[!colnames(phe2) %in% c("cluster")]

library("tableone")

tableOne<-CreateTableOne(vars=vars,strata = c("cluster"),data=phe2 ,addOverall = T,includeNA = F)

tab2<-print(tableOne,showAllLevels = T,exact =c(colnames(phe)), quote = FALSE, noSpaces = TRUE)
tab3<-as.matrix(tab2)

tableOne<-CreateTableOne(vars=vars,strata = c("cluster"),data=phe2 ,addOverall = T,includeNA = T)

tab2<-print(tableOne,showAllLevels = T,exact =c(colnames(phe)), quote = FALSE, noSpaces = TRUE)
tab3_2<-as.matrix(tab2)

tab3_combined <- tab3_2
tab3_combined[!is.na(tab3_combined[,"level"] ),] <- tab3

tab3_combined[is.na(tab3_combined[,"level"] ),c("Overall",unique(group_list))] <- stringr::str_match(string = tab3_combined[is.na(tab3_combined[,"level"] ),c("Overall",unique(group_list))],
                                                                                                     pattern = ".*(?= )")



write.csv(tab3_combined,"clinical_features_her2_low_vs_her2_pos.csv",quote=F)



## (HR+)HER2_0 VS HER2_LOW
her2_0_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_0") & FUSCC_HER2_low_project_Cohort.Info$Clinical_Subtype %in% c("HR+HER2-"),"PatientCode"]
her2_low_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_low")& FUSCC_HER2_low_project_Cohort.Info$Clinical_Subtype %in% c("HR+HER2-"),"PatientCode"]

group_list<-c(rep("HER2_0",length(her2_0_list)),rep("HER2_low",length(her2_low_list)))
names(group_list)<-c(her2_0_list,her2_low_list)
table(group_list)

phe2 <-phe[names(group_list),c("Age_at_surgery","Age_cla","Menopause","Laterality",
                               "HER2_low_status_long_RE","Histology","Grade","Ki67_percentage_cla","HR_status","sTILs","iTILs",
                               "Stage_pT","Stage_pN",
                               "TMB","HRD","PAM50","FUSCCTNBC_mRNA_Subtype",
                               "Surgery_breast","Surgery_axilla","Adjuvant_chemotherapy","Taxane_usage","Anthracyclines_usage","Platinum_usage","Capecitabine_usage","Adjuvant_radiotherapy","Adjuvant_endocrine_therapy")]


phe2$cluster <- group_list

vars<-colnames(phe2)[!colnames(phe2) %in% c("cluster")]

library("tableone")

tableOne<-CreateTableOne(vars=vars,strata = c("cluster"),data=phe2 ,addOverall = T,includeNA = F)

tab2<-print(tableOne,showAllLevels = T,exact =c(colnames(phe)), quote = FALSE, noSpaces = TRUE)
tab3<-as.matrix(tab2)

tableOne<-CreateTableOne(vars=vars,strata = c("cluster"),data=phe2 ,addOverall = T,includeNA = T)

tab2<-print(tableOne,showAllLevels = T,exact =c(colnames(phe)), quote = FALSE, noSpaces = TRUE)
tab3_2<-as.matrix(tab2)

tab3_combined <- tab3_2
tab3_combined[!is.na(tab3_combined[,"level"] ),] <- tab3

tab3_combined[is.na(tab3_combined[,"level"] ),c("Overall",unique(group_list))] <- stringr::str_match(string = tab3_combined[is.na(tab3_combined[,"level"] ),c("Overall",unique(group_list))],
                                                                                                     pattern = ".*(?= )")

write.csv(tab3_combined,"clinical_features_her2_0_vs_her2_low(HR+).csv",quote=F)


##(HR-)HER2_0 VS HER2_LOW
her2_0_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_0") & FUSCC_HER2_low_project_Cohort.Info$Clinical_Subtype %in% c("TNBC"),"PatientCode"]
her2_low_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_low")& FUSCC_HER2_low_project_Cohort.Info$Clinical_Subtype %in% c("TNBC"),"PatientCode"]

group_list<-c(rep("HER2_0",length(her2_0_list)),rep("HER2_low",length(her2_low_list)))
names(group_list)<-c(her2_0_list,her2_low_list)
table(group_list)

phe2 <-phe[names(group_list),c("Age_at_surgery","Age_cla","Menopause","Laterality",
                               "HER2_low_status_long_RE","Histology","Grade","Ki67_percentage_cla","HR_status","sTILs","iTILs",
                               "Stage_pT","Stage_pN",
                               "TMB","HRD","PAM50","FUSCCTNBC_mRNA_Subtype",
                               "Surgery_breast","Surgery_axilla","Adjuvant_chemotherapy","Taxane_usage","Anthracyclines_usage","Platinum_usage","Capecitabine_usage","Adjuvant_radiotherapy","Adjuvant_endocrine_therapy")]


phe2$cluster <- group_list

vars<-colnames(phe2)[!colnames(phe2) %in% c("cluster")]

library("tableone")

tableOne<-CreateTableOne(vars=vars,strata = c("cluster"),data=phe2 ,addOverall = T,includeNA = F)

tab2<-print(tableOne,showAllLevels = T,exact =c(colnames(phe)), quote = FALSE, noSpaces = TRUE)
tab3<-as.matrix(tab2)

tableOne<-CreateTableOne(vars=vars,strata = c("cluster"),data=phe2 ,addOverall = T,includeNA = T)

tab2<-print(tableOne,showAllLevels = T,exact =c(colnames(phe)), quote = FALSE, noSpaces = TRUE)
tab3_2<-as.matrix(tab2)

tab3_combined <- tab3_2
tab3_combined[!is.na(tab3_combined[,"level"] ),] <- tab3

tab3_combined[is.na(tab3_combined[,"level"] ),c("Overall",unique(group_list))] <- stringr::str_match(string = tab3_combined[is.na(tab3_combined[,"level"] ),c("Overall",unique(group_list))],
                                                                                                     pattern = ".*(?= )")


write.csv(tab3_combined,"clinical_features_her2_0_vs_her2_low(HR-).csv",quote=F)

