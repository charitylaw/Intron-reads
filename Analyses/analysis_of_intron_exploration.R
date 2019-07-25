# 10 Jan 2017 (Last updated 7 Jun 2018)
# Charity Law and Albert Zhang

# Fields to define
# exon.anno <- "Exon.txt" # an exon saf file
# genebody.anno <- "Genebody.txt" # a genebody saf file
# bam <- c("a.bam", "b.bam", "c.bam") # bam file names
# is.paired.end <- FALSE # seq protocol for read summarisation
# sample.names <- NULL # sample names (unique)
# groups <- sample.names # groups for samples (non-unique)

# Load libraries
library(limma)
library(edgeR)
library(Rsubread)
library(stringr)

# Create folders
dir.create("Counts", showWarnings=FALSE)
dir.create("Results", showWarnings=FALSE)
dir.create("Plots", showWarnings=FALSE)


### SUMMARISE READS TO COUNTS

# Exon
anno <- read.delim(exon.anno)
counts <- featureCounts(annot.ext=anno, files=bam, isPairedEnd=is.paired.end, useMetaFeatures=TRUE, allowMultiOverlap=FALSE, strandSpecific=0, nthreads=8)
colnames(counts$counts) <- colnames(counts$stat)[-1] <- sample.names
save(counts, file="Counts/Exon.RData")  
# Genebody
anno <- read.delim(genebody.anno)
counts <- featureCounts(annot.ext=anno, files=bam, isPairedEnd=is.paired.end, useMetaFeatures=TRUE, allowMultiOverlap=FALSE, strandSpecific=0, nthreads=8)
colnames(counts$counts) <- colnames(counts$stat)[-1] <- sample.names
save(counts, file="Counts/Genebody.RData")    
# Intron  
load("Counts/Exon.RData")
exon <- counts
load("Counts/Genebody.RData")
genebody <- counts
counts$counts <- pmax(genebody$counts-exon$counts,0)
save(counts, file="Counts/Intron.RData")


### CALC READ PERCENTAGES

# Get library size and number of exon counts
load("Counts/Exon.RData")
ntotal <- colSums(counts$stat[,-1])
nexon <- colSums(counts$counts)
# Get total intron reads
load("Counts/Intron.RData")
nintron <- colSums(counts$counts) 
# Calculate read percentages
read.percentages <- rbind(Exon=nexon, Intron=nintron)/rep(ntotal, each=2)
colnames(read.percentages) <- groups
# Save
write.table(read.percentages, file="Results/read_percentage.txt", sep="\t")
# Plot
group.f <- as.factor(groups)
pdf("Plots/read_percentage.pdf")
par(mar=c(10,5,1,1))
plot(as.numeric(group.f)+rnorm(length(group.f), sd=0.05), read.percentages["Exon",],
     ylim=c(0,1), col="dodgerblue", ylab="Proportion of reads", xlab="", xaxt="n")
points(as.numeric(group.f)+rnorm(length(group.f), sd=0.05), read.percentages["Intron",], col="magenta")
axis(side=1, at=1:nlevels(group.f), labels=levels(group.f), las=2)
legend("topleft", fill=c("dodgerblue", "magenta"), legend=rownames(read.percentages), bty="n")
dev.off()


# MDS PLOTS

# Create and save MDS plots
pdf("Plots/sample_clusters.pdf", width=8, height=4)
par(mar=c(5,5,2,1))
par(mfrow=c(1,2))
# Make exon MDS plot
load("Counts/Exon.RData")
lcpm <- cpm(counts$counts, log=TRUE)
mds.exon <- plotMDS(lcpm, labels=groups, col=as.numeric(as.factor(groups)), main="Exon")
# Make intron MDS plot
load("Counts/Intron.RData")
lcpm <- cpm(counts$counts, log=TRUE)
mds.intron <- plotMDS(lcpm, labels=groups, col=as.numeric(as.factor(groups)), main="Intron")
dev.off()
# Save 
save(mds.exon, mds.intron, file="Results/sample_clusters.RData")


# GENE SIGNAL

# Exon counts
load("Counts/Exon.RData")
exon <- DGEList(counts$counts, genes=counts$annotation)
# Intron counts
load("Counts/Intron.RData")
intron <- DGEList(counts$counts, genes=counts$annotation)
m <- match(exon$genes$GeneID, intron$genes$GeneID)
intron <- intron[m,] 
# Genes with signal in exon and intron regions
exon.signal <- exon$counts>=3
intron.signal <- intron$counts>=3
signal <- matrix(NA, nrow=nrow(exon), ncol=ncol(exon))
signal[exon.signal & intron.signal] <- "Exon and intron signal"
signal[!exon.signal & !intron.signal] <- "No signal"
signal[exon.signal & !intron.signal] <- "Exon signal"
signal[!exon.signal & intron.signal] <- "Intron signal"
# Mark genes with a single exon
nexons <- str_count(exon$genes$Chr, ";")+1
single.exon <- rep("", length(nexons))
single.exon[nexons==1] <- " (no introns)"
signal <- matrix(paste0(signal, single.exon), nrow=length(single.exon))
colnames(signal) <- sample.names
rownames(signal) <- exon$genes$GeneID
# Summarise 
signal.summarised <- apply(signal, 2, function(x) table(x))
signal.summarised <- signal.summarised/nrow(signal)
# Save
write.table(signal, file="Results/gene_signal.txt", sep="\t")
write.table(signal.summarised, file="Results/gene_signal_summarised.txt", sep="\t")
# Plot
dir.create("Plots", showWarnings=FALSE)
group <- as.factor(groups)
pdf("Plots/gene_signal.pdf")
par(mar=c(10,5,1,1))
MAX <- max(signal.summarised)+0.01
plot(as.numeric(group)+rnorm(length(group), sd=0.05), signal.summarised["Exon and intron signal",],
     ylim=c(0,MAX), col="grey", ylab="Proportion of genes", xlab="", xaxt="n")
points(as.numeric(group)+rnorm(length(group), sd=0.05), signal.summarised["Exon signal",], col="dodgerblue")
points(as.numeric(group)+rnorm(length(group), sd=0.05), signal.summarised["Intron signal",], col="magenta")
points(as.numeric(group)+rnorm(length(group), sd=0.05), signal.summarised["Exon signal (no introns)",], col="limegreen")
axis(side=1, at=1:nlevels(group), labels=levels(group), las=2)
legend("topleft", fill=c("grey", "dodgerblue", "limegreen", "magenta"), legend=c("Exon and intron signal", "Exon signal", "Exon signal (no introns)", "Intron signal"), bty="n")
dev.off()


# EXON VS INTRON LOG-COUNTS

load("Counts/Exon.RData")
exon <- DGEList(counts$counts, genes=counts$annotation)
# Intron counts
load("Counts/Intron.RData")
intron <- DGEList(counts$counts, genes=counts$annotation)
m <- match(exon$genes$GeneID, intron$genes$GeneID)
intron <- intron[m,] 
signal <- read.delim("Results/gene_signal.txt", header=TRUE)
logcounts.signal <- rep(NA, ncol(exon))
names(logcounts.signal) <- groups
for (i in 1:length(sample.names)){
  exon.i <- log2(exon$counts[,i]+1)
  intron.i <- log2(intron$counts[,i]+1)
  signal.i <- signal[,i]=="Exon and intron signal"
  logcounts.signal[i] <- cor(exon.i[signal.i], intron.i[signal.i])
}
pdf("Plots/gene_signal_correlation.pdf", height=4, width=4)
par(mar=c(10,4,1,1))
boxplot(logcounts.signal~groups, las=2)
dev.off()
# Save
write.table(logcounts.signal, file="Results/gene_signal_correlation.txt", row.names=FALSE)
