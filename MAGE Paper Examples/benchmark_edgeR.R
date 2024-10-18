# Matthew Beltran 11/3/2023
# Import expression profile from MATLAB and find differentially expressed genes using edgeR
#
# Starting with: ...
#
# cts.csv       ---> profile           (genes X samples)
# geneNames.csv ---> gene names        (genes X 0      )
# colData.csv   ---> sample names      (0     X samples)
# batch.csv     ---> batch labels      (0     X samples)
# condition.csv ---> condition labels  (0     X samples)



# read input from csv created in MATLAB

setwd("E:/MAGE paper/MAGE MATLAB/workspaces/edgeR")
colData <- read.csv("colData.csv",header = FALSE, sep = ",")
cts <- read.csv("cts.csv",header = FALSE, sep = ",")
batch <- read.csv("batch.csv",header = FALSE, sep = ",")
condition <- read.csv("condition.csv",header = FALSE, sep = ",")
geneName <- read.csv("geneName.csv",header = FALSE, sep = ",")



# create foactor of group condition labels
grps <- as.factor(condition$V1)




# create edgeR object DGEList
y <- DGEList(counts=cts,group=grps)

# filtering
#keep <- filterByExpr(y)
#y <- y[keep,,keep.lib.sizes=FALSE]

y <- normLibSizes(y)
design <- model.matrix(~grps)
y <- estimateDisp(y,design)

logcpm <- cpm(y, log=TRUE)
plotMDS(y)
plotBCV(y)


# perform likelihood ratio tests:
fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef=2)
topTags(lrt)

plotMD(lrt)


# perform quasi-likelihood F-tests:
# more conservative p-vals/FDR
fit <- glmQLFit(y,design)
qlf <- glmQLFTest(fit,coef=2)
topTags(qlf)

plotQLDisp(fit)
plotMD(qlf)

# export results to MATLAB (.csv)
exportResults <- data.frame(gene_name = geneName[row.names(lrt[["coefficients"]]),1],
                            logFC = lrt[["table"]][["logFC"]],
                            logCPM = lrt[["table"]][["logCPM"]],
                            L_ratio = lrt[["table"]][["LR"]],
                            p_val = lrt[["table"]][["PValue"]],
                            qF_test = qlf[["table"]][["F"]],
                            qF_p_val = qlf[["table"]][["PValue"]])
write.csv(exportResults,"benchmark_edgeR.csv",row.names = FALSE)


