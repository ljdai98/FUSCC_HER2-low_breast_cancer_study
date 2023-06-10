library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)

DEG <- read.csv("HR_neg_HER2_0_low_deg.csv")
# DEG <- read.csv("HR_pos_HER2_0_low_deg.csv")

DEG <- DEG[is.na(DEG$padj)==F,]

up <- DEG[DEG$log2FoldChange >= 1 & DEG$padj < 0.05,1]
down <- DEG[DEG$log2FoldChange <= -1 & DEG$padj < 0.05,1]

gene <- up
# gene <- down

gene <- mapIds(org.Hs.eg.db,gene,column = 'ENTREZID', keytype = 'SYMBOL',multiVals = 'filter')

GO <- enrichGO(gene=na.omit(gene),OrgDb=org.Hs.eg.db,ont="ALL",pAdjustMethod = "fdr",
               pvalueCutoff=0.05,qvalueCutoff  = 0.05) %>% as.data.frame()

GO <- GO[order(-GO$Count),][1:10,]
GO$number <- factor(rev(1:nrow(GO)))

ggplot(data=GO, aes(x=number,y=Count))+
  geom_bar(stat="identity", width=0.8) + coord_flip()+ 
  theme_test()+ 
  scale_x_discrete(labels=rev(GO$Description)) +
  theme(axis.text=element_text(size = 10,face = "bold", color="black")) +            
  labs(title = "Enriched GO terms")
