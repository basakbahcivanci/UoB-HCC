rm(list = ls())
library(limma)
library(EnhancedVolcano)

setwd('/Users/user/Documents/MSc BIOINFORMATICS/Individual Project/HCCDB4')
data = read.csv('recent_data4.csv', header = TRUE,row.names = 1)
sample_info = read.csv('recent_sample_info4.csv', header = TRUE, row.names = 1)



ct <- factor(sample_info$TYPE)
table(ct)
# Create Design Matrix 
design <- model.matrix(~0+ct)
# Need renaming the column
colnames(design) <- levels(ct)

fit <- lmFit(data,design)

contrasts <- makeContrasts(Diff = HCC-Adjacent, levels=design)


fit2 <- contrasts.fit(fit, contrasts)
fit3 <- eBayes(fit2, trend=TRUE)
summary(decideTests(fit3, method="global"))
deg = topTable(fit3, coef='Diff',number = nrow(data))



pdf(file = "HCCDB4 - Volcano.pdf")
EnhancedVolcano(deg,
                lab = rownames(deg),
                x = 'logFC',
                y = 'adj.P.Val',
                title = 'Adjacent versus HCC',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 1.5,
                labSize = 3.0)
dev.off()

deg = topTable(fit3, coef='Diff',p.value = 0.05, lfc = 1, number = nrow(data)) # Min fold change is 2, padjcutoff=0.05
write.csv(deg,'HCCDB4_FC2_DEGS.csv', quote = F,row.names = TRUE)

down = deg[which(deg$logFC<=-1),]
write.csv(down,'HCCDB4_FC2_DEGS_DOWN.csv', quote = F,row.names = TRUE)

up= deg[which(deg$logFC>=1),]
write.csv(up,'HCCDB4_FC2_DEGS_UP.csv', quote = F,row.names = TRUE)

degs = data[rownames(data) %in% rownames(deg),]
write.csv(degs,'DEGS_HCCDB4.csv', quote = F,row.names = TRUE)
