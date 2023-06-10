table(FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE)
# HER2_0      HER2_low HER2_positive 
# 91           434           182 
her2_0_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_0") & FUSCC_HER2_low_project_Cohort.Info$Clinical_Subtype %in% c("TNBC"),"PatientCode"]
her2_low_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_low")& FUSCC_HER2_low_project_Cohort.Info$Clinical_Subtype %in% c("TNBC"),"PatientCode"]

group_list<-c(rep("HER2_0",length(her2_0_list)),rep("HER2_low",length(her2_low_list)))
names(group_list)<-c(her2_0_list,her2_low_list)
table(group_list)

common_tnbc_subtype<-FUSCC_HER2_low_project_Cohort.Info[names(group_list),c("PAM50","FUSCCTNBC_mRNA_Subtype","TNBC_Burstein_subtype","TNBC_Lehmann_Subtype","TNBC_Quist_subtype")]
common_tnbc_subtype$Category<-group_list

common_tnbc_subtype$PAM50 <- factor(common_tnbc_subtype$PAM50,levels = c("LumA","LumB","Her2","Basal","Normal"))
common_tnbc_subtype$TNBC_Burstein_subtype <- factor(common_tnbc_subtype$TNBC_Burstein_subtype,levels = c("LAR","MES","BLIA","BLIS"))
common_tnbc_subtype$TNBC_Lehmann_Subtype <- factor(common_tnbc_subtype$TNBC_Lehmann_Subtype,levels = c("LAR","UNS","IM","M","BL2","BL1"))
common_tnbc_subtype$FUSCCTNBC_mRNA_Subtype <- factor(common_tnbc_subtype$FUSCCTNBC_mRNA_Subtype,levels = c("LAR","IM","BLIS","MES"))


library("tableone")
vars<-colnames(common_tnbc_subtype)[!colnames(common_tnbc_subtype) %in% c("Category")]
tableOne<-CreateTableOne(vars=vars,strata = c("Category"),data=common_tnbc_subtype )

tab2<-print(tableOne,showAllLevels = T, exact = vars, quote = FALSE, noSpaces = TRUE)
tab3<-as.matrix(tab2)


tableOne<-CreateTableOne(vars=vars,strata = c("Category"),data=common_tnbc_subtype,includeNA = T )

tab2<-print(tableOne,showAllLevels = T, exact = vars, quote = FALSE, noSpaces = TRUE)
tab3_2<-as.matrix(tab2)


tab3_combined <- tab3_2
tab3_combined[!is.na(tab3_combined[,"level"] ),] <- tab3

tab3_combined[is.na(tab3_combined[,"level"] ),unique(group_list)] <- stringr::str_match(string = tab3_combined[is.na(tab3_combined[,"level"] ),unique(group_list)],
                                                                                    pattern = ".*(?= )")
write.csv(tab3_combined,"Molecular_subtype_in_TNBCs_combined.csv",quote=F)

