table(FUSCC_HER2_low_project_Cohort.Info$HER2_low_status_RE)
# HER2_0      HER2_low HER2_positive 
# 91           434           182 

for_three <- data.frame(Peak_17q11.2_val=as.numeric(FUSCC_HER2_low_project_GISTICpeaks.val["Del_Peak.42_17q11.2",]),
                        Peak_17q21.31_val=as.numeric(FUSCC_HER2_low_project_GISTICpeaks.val["Del_Peak.43_17q21.31",]),
                        Peak_17q12_val=as.numeric(FUSCC_HER2_low_project_GISTICpeaks.val["Amp_Peak.21_17q12",]))
rownames(for_three) <- colnames(FUSCC_HER2_low_project_GISTICpeaks.val)
for_three <- for_three[FUSCC_HER2_low_project_Cohort.Info[rownames(for_three),"HR_status"] %in% c("Positive") & FUSCC_HER2_low_project_Cohort.Info[rownames(for_three),"HER2_low_status_RE"] %in% c("HER2_low","HER2_0") ,]


library(ggplot2)
library(ggpubr)

dat=data.frame(for_three)
pdf(file="02_correlation_dotplot_of_Peak_17q11.2_val_Peak_17q21.31_val.pdf",width = 4,height = 4)
ggplot(data=dat, aes(x=Peak_17q11.2_val, y=Peak_17q21.31_val),main="LOH")+geom_point(color="red")+stat_smooth(method="lm",se=FALSE)+stat_cor(data=dat, method = "spearman",cor.coef.name="rho") #有序等级资料用spearman(rho)，连续资料用pearson(R)
dev.off()

dat=data.frame(for_three)
pdf(file="02_correlation_dotplot_of_Peak_17q11.2_val_Peak_17q12_val.pdf",width = 4,height = 4)
ggplot(data=dat, aes(x=Peak_17q11.2_val, y=Peak_17q12_val),main="LOH")+geom_point(color="red")+stat_smooth(method="lm",se=FALSE)+stat_cor(data=dat, method = "spearman",cor.coef.name="rho") #有序等级资料用spearman(rho)，连续资料用pearson(R)
dev.off()

dat=data.frame(for_three)
pdf(file="02_correlation_dotplot_of_Peak_17q21.31_val_Peak_17q12_val.pdf",width = 4,height = 4)
ggplot(data=dat, aes(x=Peak_17q21.31_val, y=Peak_17q12_val),main="LOH")+geom_point(color="red")+stat_smooth(method="lm",se=FALSE)+stat_cor(data=dat, method = "spearman",cor.coef.name="rho") #有序等级资料用spearman(rho)，连续资料用pearson(R)
dev.off()


