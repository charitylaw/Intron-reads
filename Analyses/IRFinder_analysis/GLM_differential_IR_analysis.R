# 11 Sep 2018 (Last updated 15 July 2019)
# Charity Law
# IRFinder analysis using GLM test from DESeq2 (local paths removed)

# Set up
library(DESeq2)
library(limma)
source("./bin/DESeq2Constructor.R")  #Load IRFinder-related function

# GLM on IRFinder results tables
results = read.table("filePaths.txt")
paths = as.vector(results$V1)
experiment = read.table("experiment.txt",header=T)                       
experiment$Group = as.factor(experiment$Group)
rownames(experiment)=NULL       
metaList=DESeqDataSetFromIRFinder(filePaths=paths, designMatrix=experiment, designFormula=~1)
dds = metaList$DESeq2Object
colData(dds)    
design(dds) = ~Group + Group:IRFinder
dds = DESeq(dds)                      
resultsNames(dds)                     
res.diff = results(dds, contrast=list("GroupHCC827_PolyA.IRFinderIR","GroupNCI.H1975_PolyA.IRFinderIR"))

# Order results by adjusted p-value
o <- order(res.diff$padj, decreasing=FALSE)
res.top <- res.diff[o,]
write.table(res.top, file="IR_results.txt", sep="\t")