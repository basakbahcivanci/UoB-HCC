rm(list = ls())
library(limma)
library(EnhancedVolcano)

# reading and storing data HCCDB6
setwd('/Users/user/Documents/MSc BIOINFORMATICS/Individual Project/HCCDB6/GSE14520(GPL3721 Subset)_Affymetrix')
data = read.csv('prepared_data6.csv', header = TRUE,row.names = 1)
sample_info = read.csv('prepared_sample_info6.csv', header = TRUE, row.names = 1)



# assigning type as a foctor for experimental design matrix
ct <- factor(sample_info$TYPE)

# create Design Matrix 
design <- model.matrix(~0+ct)
# renaming the column
colnames(design) <- levels(ct)
# fitting linear model for each gene
fit <- lmFit(data,design)
# setting up comparisons of interest
contrasts <- makeContrasts(Diff = HCC-Adjacent, levels=design)
# to estimate original fit and the difference between groups
fit2 <- contrasts.fit(fit, contrasts)
# calculates moderated t-statistics
fit3 <- eBayes(fit2, trend=TRUE)
# summary(decideTests(fit3, method="global"))


deg = topTable(fit3, coef='Diff',number = nrow(data)) 


# quartz() # for polygon edge not found error
pdf(file = "HCCDB6 - Volcano.pdf")
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

deg = topTable(fit3, coef='Diff',p.value = 0.05, lfc = 1, number = nrow(data)) # Min log fold change is 2, padjcutoff=0.05
write.csv(deg,'HCCDB6_FC2_DEGS.csv', quote = F,row.names = TRUE)

down = deg[which(deg$logFC<=-1),]
write.csv(down,'HCCDB6_FC2_DEGS_DOWN.csv', quote = F,row.names = TRUE)

up= deg[which(deg$logFC>=1),]
write.csv(up,'HCCDB6_FC2_DEGS_UP.csv', quote = F,row.names = TRUE)

# filtering DEGs from data and writing it out
degs = data[rownames(data) %in% rownames(deg),]
write.csv(degs,'DEGS_HCCDB6.csv', quote = F,row.names = TRUE)

