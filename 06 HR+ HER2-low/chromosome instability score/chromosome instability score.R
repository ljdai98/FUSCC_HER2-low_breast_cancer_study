table(FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE_RE)
# HER2_0      HER2_low HER2_positive 
# 91           434           182 


##
her2_0_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_0"),"PatientCode"]
her2_low_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_low"),"PatientCode"]
her2_pos_list<-FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE %in% c("HER2_positive"),"PatientCode"]

group_list<-c(rep("HER2_0",length(her2_0_list)),rep("HER2_low",length(her2_low_list)),rep("HER2_positive",length(her2_pos_list)))
names(group_list)<-c(her2_0_list,her2_low_list,her2_pos_list)
table(group_list)

##################################
dataset <- FUSCC_HER2_low_project_Cohort.Info[,c("HER2_low_status_RE","HRD","LOH","TAI","LST","AiCNA")]
dataset <- dataset[FUSCC_HER2_low_project_Cohort.Info[FUSCC_HER2_low_project_Cohort.Info$Clinical_Subtype %in% "HR+HER2-","PatientCode"],]
dataset <- dataset[rownames(dataset) %in% colnames(FUSCC_HER2_low_project_GISTICpeaks.del.thre),]

dataset$copy_number <- as.character(FUSCC_HER2_low_project_GISTICpeaks.del.thre["Del_Peak.43_17q21.31",rownames(dataset)])


dataset[,ncol(dataset)] <- gsub("1|2","Loss/del",dataset[,ncol(dataset)])
dataset[,ncol(dataset)] <- gsub("0","Others",dataset[,ncol(dataset)])

dataset$Combine_group <- paste(dataset$HER2_low_status_RE,dataset[,ncol(dataset)],sep = "_")
table(dataset$Combine_group)


my_comparisons <- list( c("HER2_0_Loss/del","HER2_0_Others"),
                        c("HER2_low_Loss/del","HER2_low_Others"))


palette <- c("#EE0000","#FF9076","#003889","#5290FF")
names(palette) <-  c("HER2_0_Loss/del","HER2_0_Others","HER2_low_Loss/del","HER2_low_Others")
  
pdf("HR_pos_chromosome_instability_score_HER2_low_status_17q21.31.pdf",width = 2,height = 5)
for (i in c("HRD","LOH","TAI","LST","AiCNA")){
  library(ggplot2)
  library(ggpubr)

  tempplot <- ggboxplot(dataset, x = "Combine_group", y = i,
                        color = "Combine_group",palette = palette,order=c("HER2_0_Loss/del","HER2_0_Others","HER2_low_Loss/del","HER2_low_Others"),
                        add = "jitter",title = i)+stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")
  print(tempplot)

}
dev.off()

