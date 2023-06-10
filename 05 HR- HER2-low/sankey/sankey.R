table(FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE)
# HER2_0      HER2_low HER2_positive 
# 91           434           182 

FUSCC_HER2_low_project_Cohort.Info$TNBC_Burstein_subtype <- factor(FUSCC_HER2_low_project_Cohort.Info$TNBC_Burstein_subtype,levels = c("MES","LAR","BLIA","BLIS"))
FUSCC_HER2_low_project_Cohort.Info$TNBC_Lehmann_Subtype <- factor(FUSCC_HER2_low_project_Cohort.Info$TNBC_Lehmann_Subtype,levels = c("UNS","LAR","IM","M","BL2","BL1"))
FUSCC_HER2_low_project_Cohort.Info$TNBC_Quist_subtype <- factor(FUSCC_HER2_low_project_Cohort.Info$TNBC_Quist_subtype,levels = c("MC1","MC2","MC3","MC4","MC5","MC6"))


HR_neg_HER2_low_list_basal<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HR_status %in% c("Negative") & FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_low") & FUSCC_HER2_low_project_Cohort.Info$PAM50 %in% c("Basal"),"PatientCode"]
HR_neg_HER2_low_list_nonbasal<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HR_status %in% c("Negative") & FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_low") & FUSCC_HER2_low_project_Cohort.Info$PAM50 %in% c("Her2","LumA","LumB","Normal"),"PatientCode"]

group_list<-c(rep("HER2_low_basal",length(HR_neg_HER2_low_list_basal)),
              rep("HER2_low_nonbasal",length(HR_neg_HER2_low_list_nonbasal)))

names(group_list)<-c(HR_neg_HER2_low_list_basal,HR_neg_HER2_low_list_nonbasal)
group_list<-factor(group_list,levels=c("HER2_low_basal","HER2_low_nonbasal"))
table(group_list)
# group_list
# HER2_low_basal HER2_low_nonbasal 
# 46                20

##FUSCCTNBC_mRNA_Subtype
common_tnbc_subtype<-FUSCC_HER2_low_project_Cohort.Info[names(group_list),c("FUSCCTNBC_mRNA_Subtype","PAM50")]

common_tnbc_subtype$Category<-group_list
common_tnbc_subtype<-na.omit(common_tnbc_subtype)
x<-data.frame(table(common_tnbc_subtype))

library(ggalluvial)

pdf("sankey_FUSCCTNBC_mRNA_Subtype.pdf",width = 5.5,height = 4)
ggplot(x,
       aes(y =Freq,
           axis1 = FUSCCTNBC_mRNA_Subtype, axis2 = PAM50, axis3 =Category ))+
  geom_alluvium(aes(fill = Category),width = 0, reverse = FALSE)+ 
  geom_stratum(width = 1/3, reverse = FALSE) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),reverse = FALSE, size = 4,angle=0)+ 
  scale_x_continuous(breaks = 1:3, labels = c("PAM50", "FUSCCTNBC_mRNA_Subtype","Category"))+
  theme(legend.position = "none") +
  coord_flip()
dev.off()


##TNBC_Lehmann_Subtype
common_tnbc_subtype<-FUSCC_HER2_low_project_Cohort.Info[names(group_list),c("TNBC_Lehmann_Subtype","PAM50")]

common_tnbc_subtype$Category<-group_list
common_tnbc_subtype<-na.omit(common_tnbc_subtype)
x<-data.frame(table(common_tnbc_subtype))

x$TNBC_Lehmann_Subtype <- factor(x$TNBC_Lehmann_Subtype ,levels = c("UNS","IM","BL2","BL1","LAR","M"))


library(ggalluvial)

pdf("sankey_TNBC_Lehmann_Subtype.pdf",width = 5.5,height = 4)
ggplot(x,
       aes(y =Freq,
           axis1 = TNBC_Lehmann_Subtype, axis2 = PAM50, axis3 =Category ))+
  geom_alluvium(aes(fill = Category),width = 0, reverse = FALSE)+ 
  geom_stratum(width = 1/3, reverse = FALSE) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),reverse = FALSE, size = 4,angle=0)+ 
  scale_x_continuous(breaks = 1:3, labels = c("PAM50", "Lehmann_Subtype","Category"))+
  theme(legend.position = "none") +
  coord_flip()
dev.off()


##TNBC_Burstein_subtype
common_tnbc_subtype<-FUSCC_HER2_low_project_Cohort.Info[names(group_list),c("TNBC_Burstein_subtype","PAM50")]

common_tnbc_subtype$Category<-group_list
common_tnbc_subtype<-na.omit(common_tnbc_subtype)
x<-data.frame(table(common_tnbc_subtype))

x$TNBC_Burstein_subtype <- factor(x$TNBC_Burstein_subtype ,levels = rev(c("LAR","MES","BLIA","BLIS")))


library(ggalluvial)

pdf("sankey_TNBC_Burstein_subtype.pdf",width = 5.5,height = 4)
ggplot(x,
       aes(y =Freq,
           axis1 = TNBC_Burstein_subtype, axis2 = PAM50, axis3 =Category ))+
  geom_alluvium(aes(fill = Category),width = 0, reverse = FALSE)+ 
  geom_stratum(width = 1/3, reverse = FALSE) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),reverse = FALSE, size = 4,angle=0)+ 
  scale_x_continuous(breaks = 1:3, labels = c("PAM50", "TNBC_Burstein_subtype","Category"))+
  theme(legend.position = "none") +
  coord_flip()
dev.off()

