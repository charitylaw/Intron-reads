# Charity Law
# 12 Sep 2017 (Last updated 6 June 2018)


# Define fields from targets file
targets <- read.delim("targets_clean.txt")
targets <- targets[targets$Cellline=="HCC827",]
bam <- targets$File
is.paired.end <- FALSE
sample.names <- as.character(targets$SampleName)
groups <- as.character(targets$RNAPrep)

# Create folders
dir.create("coverage_for_HCC827", showWarnings=FALSE)

# Lood libraries
library(limma)
library(edgeR)
library(Rsubread)
library(R.utils)
source("FUN_create_binned_annotation.R")
source("FUN_plot_coverage_profile.R")

### SELECT GENES TO ANALYSE

# Select genes with exon and intron signal in HCC827 poly(A) RNA samples
genes.signal <- read.table("Results/gene_signal.txt", header=TRUE)
genes.signal <- genes.signal[,grep("HCC827$", colnames(genes.signal))]
keep.genes.signal <- rowSums(genes.signal=="Exon and intron signal")==ncol(genes.signal)
genes.signal <- rownames(genes.signal)[keep.genes.signal]
length(genes.signal)
# Select genes with no overlap with other genes
genes.solitary <- read.table("Genebody_nonoverlapping.txt", header=TRUE)
# Genes excluding those on chrM
genes.solitary.chr <- genes.solitary$GeneID[genes.solitary$Chr!="chrM"]
# Genes that are protein coding
genes.type <- read.table("Gene_information.txt", header=TRUE)
genes.proteincoding <- genes.type$GeneID[genes.type$Type=="protein_coding"]
# Select genes
genes <- intersect(genes.signal, genes.solitary.chr)
genes <- intersect(genes, genes.proteincoding)
length(genes)
write(genes, file="coverage_for_HCC827/select_genes.txt")


# CREATE BINNED ANNOTATION 

# Roughly 7 mins to create annotation for Exons; 74 mins for Introns. 
# Load function that makes binned annotation
# Create binned annotation for Exon and Intron
for (feature in c("Exon", "Intron")){
  cat("Creating binned annotation for", feature, "features.\n")
  # Load feature annotation
  anno <- read.table(paste0(feature, ".txt"), header=TRUE)
  # Select genes
  genes <- read.table("coverage_for_HCC827/select_genes.txt")[,1]
  anno <- anno[anno$GeneID %in% genes,]
  cat("\t", "Number of features to process:", nrow(anno), "\n")
  cat("\t", "Number of genes to process:", length(unique(anno$GeneID)), "\n")
  # Create annotation
  start.time <- Sys.time()
  anno.binned <- createSubanno(anno, feature=toupper(feature), binsize=30)
  run.time <- difftime(Sys.time(), start.time, units="secs")
  # Order new annotation by position
  o <- order(anno.binned$Chr, anno.binned$Start)
  anno.binned <- anno.binned[o,]
  rownames(anno.binned) <- NULL
  # Save 
  write.table(anno.binned, file=paste0("coverage_for_HCC827/", feature, "_binned.txt"))
  cat("\t Time taken to create binned annotation:", round(run.time/60,2), "mins \n")
}


# SUMMARISE READS

for (feature in c("Exon", "Intron")){
  cat("Creating subcounts for", feature, "features.\n")
  # Binned annotation
  anno <- read.table(paste0("coverage_for_HCC827/", feature, "_binned.txt"))
  anno$FeatureID <- anno$GeneID
  genes <- strsplit2(anno$GeneID, split=paste0("_", toupper(feature), "_"))
  anno$GeneID <- genes[,1]
  anno$FeatureNumber <- genes[,2]
  anno$FeatureType <- feature
  anno <- anno[,c("GeneID", "FeatureID", "FeatureType", "FeatureNumber", "Chr", "Start", "End", "Strand", "Length")]
  # Summarise reads (overlapping reads are multi-counted)
  start.time <- Sys.time()
  counts <- featureCounts(files=bam, annot.ext=anno, useMetaFeatures=FALSE, isPairedEnd=is.paired.end, allowMultiOverlap=TRUE, nthreads=8)
  run.time <- difftime(Sys.time(), start.time, units="hours")
  colnames(counts$counts) <- colnames(counts$stat)[-1] <- sample.names
  # Clean up as DGEList object
  x <- DGEList(counts$counts)
  x$genes <- anno
  x$genes$Length <- x$genes$Length+1
  # Save results to file
  savefile <- paste0("coverage_for_HCC827/", feature, "_coverage.RData")
  save(x, file=savefile)
  # Print time
  cat(paste("Time taken to run featureCounts:", round(run.time, 2), "hours."))
}



# DISTRIBUTION OF COVERAGE

# Load binned coverage for exon and intron 
load("coverage_for_HCC827/Exon_coverage.RData")
x.exon <- x
y.exon <- log2(x.exon$counts+1)
geneid.exon <- x.exon$genes$GeneID
load("coverage_for_HCC827/Intron_coverage.RData")
x.intron <- x
y.intron <- log2(x.intron$counts+1)
geneid.intron <- x.intron$genes$GeneID
nsamples <- length(sample.names)
# Gene trimmed mean for log-coverage
tmean.exon <- matrix(NA, nrow=length(unique(geneid.exon)), ncol=nsamples)
rownames(tmean.exon) <- sort(unique(geneid.exon))
colnames(tmean.exon) <- sample.names
for (j in 1:nsamples) tmean.exon[,j] <- tapply(y.exon[,j], geneid.exon, function(x) mean(x, trim=0.1))
tmean.intron <- matrix(NA, nrow=length(unique(geneid.intron)), ncol=nsamples)
rownames(tmean.intron) <- sort(unique(geneid.intron))
colnames(tmean.intron) <- sample.names
for (j in 1:nsamples) tmean.intron[,j] <- tapply(y.intron[,j], geneid.intron, function(x) mean(x, trim=0.1))
# Gene number of bins
nbins.exon <- as.numeric(table(geneid.exon))
names(nbins.exon) <- sort(unique(geneid.exon))
Nbins.exon <- rep("regular", length(nbins.exon))
Nbins.exon[nbins.exon<quantile(nbins.exon, 0.33)] <- "short"
Nbins.exon[nbins.exon>quantile(nbins.exon,0.67)] <- "long"
names(Nbins.exon) <- names(nbins.exon)
table(Nbins.exon)
Nbins.exon <- factor(Nbins.exon, levels=c("short", "regular", "long"))
nbins.intron <- as.numeric(table(geneid.intron))
names(nbins.intron) <- sort(unique(geneid.intron))
Nbins.intron <- rep("regular", length(nbins.intron))
Nbins.intron[nbins.intron<quantile(nbins.intron, 0.33)] <- "short"
Nbins.intron[nbins.intron>quantile(nbins.intron,0.67)] <- "long"
names(Nbins.intron) <- names(nbins.intron)
table(Nbins.intron)
Nbins.intron <- factor(Nbins.intron, levels=c("short", "regular", "long"))
# Median length of categories
tapply(nbins.exon, Nbins.exon, median)*30
tapply(nbins.intron, Nbins.intron, median)*30
# Plot by replicate
for (replicate in c("R1", "R2", "R3")){
pdf.name <- paste0("coverage_for_HCC827/coverage_distribution_and_length_", replicate, ".pdf")
pdf(pdf.name, height=2.5, width=7.5)
l <- layout(
  matrix(c(1,2,3), nrow=1, ncol=3, byrow=FALSE),
  width=c(1,2,0.5)
)
par(mar=c(2,4,1,1))
# Select samples
i.select <- grep(replicate, sample.names)
# Plot of coverage
boxplot(y.exon[,i.select[1]], main="Coverage of bins", 
        ylim=c(-0.8,13), ylab="Log-coverage", xlim=c(1-0.05,1.35+0.05), outline=0, boxwex=0.1)
axis(side=1, at=c(1.05, 1.3), labels=c("Exon bin", "Intron bin"), tick=FALSE, cex=0.7, line=-1)
boxplot(y.intron[,i.select[1]], add=TRUE, at=1.25, outline=0, boxwex=0.1)
boxplot(y.exon[,i.select[2]], add=TRUE, at=1.1, col="grey", outline=0, boxwex=0.1)
boxplot(y.intron[,i.select[2]], add=TRUE, at=1.35, col="grey", outline=0, boxwex=0.1)
# Plot of average coverage
bp <- boxplot(tmean.exon[,i.select[1]]~Nbins.exon, main="Average coverage of  bins",
              ylab="Average gene log-coverage", xlab="", ylim=c(-0.5,10), xlim=c(0,17),
              xaxt="n", outline=0, border=c("indianred1", "indianred3", "indianred4"))
boxplot(tmean.exon[,i.select[2]]~Nbins.exon, add=TRUE, at=5:7, 
        xaxt="n", col="grey", outline=0, border=c("indianred1", "indianred3", "indianred4"))
bp <- boxplot(tmean.intron[,i.select[1]]~Nbins.intron, add=TRUE, at=10:12, 
              xaxt="n", outline=0, border=c("indianred1", "indianred3", "indianred4"))
boxplot(tmean.intron[,i.select[2]]~Nbins.intron, add=TRUE, at=14:16, 
        xaxt="n", col="grey", outline=0, border=c("indianred1", "indianred3", "indianred4"))
axis(side=1, at=c(4,13), labels=c("Exon bin", "Intron bin"), tick=FALSE, cex=0.7, line=-1)
# Legend
par(mar=c(2,0,1,0))
plot(1:5, rep(1,5), 
     type="n", xaxt="n", yaxt="n", ylab="", xlab="", bty="n")
legend("topright", 
       legend=c("short region", "regular region", "long region", "poly(A) RNA", "Total RNA"),
       fill=c(NA, NA, NA, "white", "grey"), border=rep(c("white", "black"), c(3,2)),
       text.col=c("indianred1", "indianred3", "indianred4", "black", "black"), cex=1,
       bty="n")
dev.off()
}


# COVERAGE ALONG GENEBODY
# Coverage over genebody separated into nsections
nsections <- 20
genes <- sort(unique(geneid.exon))
ngenes <- length(genes)
for(feature in c("exon", "intron")){
  # Set up objects
  coverage <- as.list(rep(NA, nsamples))
  for(s in 1:nsamples) {
    coverage[[s]] <- matrix(NA, nrow=ngenes, ncol=nsections)
    rownames(coverage[[s]]) <- genes
  }
  names(coverage) <- sample.names
  geneid <- y <- strand <- Nbins <- NA
  assign("geneid", get(paste0("geneid.", feature)))
  assign("y", get(paste0("y.", feature)))
  assign("strand", get(paste0("x.", feature))$genes$Strand)
  assign("Nbins", get(paste0("Nbins.", feature)))
  
  # Get adjusted log-coverage for each gene
  for(i in 1:ngenes){
    cat(i, "of", ngenes, "\n")
    # Log-coverage for this gene
    gene.i <- genes[i]
    select <- geneid==gene.i
    strand.i <- as.character(strand[select][1])
    y.i <- y[select,,drop=FALSE]
    # Adjust log-counts relative to max log-count in the gene
    y.max <- apply(y.i, 2, max)
    for (s in 1:nsamples) if(y.max[s]>0) y.i[,s] <- y.i[,s]/y.max[s]
    # Split bins into nsections (symmetric bin lengths) and get the mean log-coverage of each section
    n <- nrow(y.i)
    if(nsections<=n){#} & all(y.max>1)){
      sectionlength <- floor(n/nsections)
      l <- rep(sectionlength, nsections)
      midbin <- floor(nsections/2)
      extra <- n-sectionlength*nsections
      if(extra>0) l[midbin] <- l[midbin]+extra
      end <- cumsum(l)
      start <- c(1,end[1:(nsections-1)]+1)
      for(s in 1:nsamples){
        coverage.s <- rep(NA, nsections)  
        for (ns in 1:nsections) coverage.s[ns] <- mean(y.i[start[ns]:end[ns],s])
        if(strand.i=="-") coverage.s <- coverage.s[nsections:1]
        coverage[[sample.names[s]]][gene.i,] <- coverage.s
      }  
    }
  }
  # Get means for each section separating short, regular and long genes
  coverage.means.short <- lapply(coverage, function(x) colMeans(x[Nbins=="short",], na.rm=TRUE))
  coverage.means.regular <- lapply(coverage, function(x) colMeans(x[Nbins=="regular",], na.rm=TRUE))
  coverage.means.long <- lapply(coverage, function(x) colMeans(x[Nbins=="long",], na.rm=TRUE))
  # Save results
  savefile <- paste0("coverage_for_HCC827/coverage_across_body_of_", feature, ".RData")
  save(coverage, coverage.means.short, coverage.means.regular, coverage.means.long, nsections, file=savefile)
  # Plot results
  pdfname <- paste0("coverage_for_HCC827/coverage_across_body_", feature, ".pdf")
  pdf(pdfname, height=1.5, width=7.5)
  cols <- c("indianred1", "indianred3", "indianred4")
  names(cols) <- c("short", "regular", "long")
  par(mfrow=c(1,4))
  if(feature=="exon") YLIM <- c(0.15, 0.8)
  if(feature=="intron") YLIM <- c(0.01,0.50)
  for(length in c("short", "regular", "long")){
    if(length=="short") par(mar=c(2,2,1,0))
    if(length=="regular") par(mar=c(2,1,1,1))
    if(length=="long") par(mar=c(2,0,1,2))
    plot(nsections, nsections, ylim=YLIM, xaxt="n", yaxt="n", xlim=c(1,nsections), xlab="", ylab="")
    if(length=="short" & feature=="exon") axis(side=2, at=c(0.2,0.4,0.6,0.8), labels=c(0.2,0.4,0.6,0.8))
    if(length=="short" & feature=="intron") axis(side=2, at=c(0, 0.2,0.4), labels=c(0,0.2,0.4))
    title(xlab="5' ------Genebody------ 3'", line=0.5, cex=1)
    for(sn in sample.names){
      if(length(grep("_T", sn))>0) {
        col <- "grey"} else {
          col <- "black"
        }
      lines(1:nsections, get(paste0("coverage.means.", length))[[sn]], col=col)
    }
    legend("topleft", paste0(length, " regions"), bty="n", text.col=cols[length])
    # if(length=="short") title(ylab="Adjusted log-coverage", line=2)
  }
  # legend
  plot(1:5, rep(1,5), type="n", xaxt="n", yaxt="n", ylab="", xlab="", bty="n")
  legend("topleft", legend=c("poly(A) RNA", "Total RNA"), col=c("black", "grey"), lty=1, bty="n")
  dev.off()
}


# COVERAGE PROFILES FOR TWO SHORT GENES AND TWO LONG GENES

pdf("coverage_for_HCC827/coverage_profile_of_two_short_and_two_long_genes.pdf", height=6)
par(mfrow=c(4,2))
cutoff.logcoverage <- 3
G <- c("ENSG00000136997.17", "ENSG00000117525.13", "ENSG00000196937.10", "ENSG00000196428.12")
m <- match(G, genes.type$GeneID)
G <- genes.type[m,]
m <- match(G$GeneID, x.exon$genes$GeneID)
G$Strand <- x.exon$genes$Strand[m]
G
# r <- c(FALSE, TRUE, TRUE, FALSE)
y <- c(10,12,8,8)
for (i in 1:nrow(G)){
  g <- as.character(G$GeneID[i])
  if(G$Strand[i]=="-") {
    r <- TRUE
  } else {
      r <- FALSE
      }
  plotSubcounts(gene=g, exon=x.exon[,groups=="PolyA"], intron=x.intron[,groups=="PolyA"], intron.high=cutoff.logcoverage, groups=rep("",3), save.to.file=FALSE, cols = c("grey", "dodgerblue", "gray40"), reverse=r, ylim=y[i])
  plotSubcounts(gene=g, exon=x.exon[,groups=="Total"], intron=x.intron[,groups=="Total"], intron.high=cutoff.logcoverage, groups=rep("",3), save.to.file=FALSE, cols = c("grey", "dodgerblue", "gray40"), reverse=r, ylim=y[i])
}
dev.off()


# COVERAGE PROFILE OF IR-LIKE GENES IN R1 POLYA SAMPLE

# Remove genes with very high median log-coverage in introns
cutoff.intronbackground <- 2
remove <- rownames(tmean.intron)[tmean.intron[,"R1.HCC827"]>cutoff.intronbackground]
length(remove)
y.intron.adj <- y.intron[,"R1.HCC827"]
y.intron.adj[rownames(y.intron) %in% remove][] <- 0
# Select genes with extreme intron bins
cutoff.logcoverage <- 6
y.extreme <- y.intron.adj>cutoff.logcoverage
genes.extreme <- names(y.extreme)[y.extreme]
# Genes ordered by most number of extreme bins to least
genes.extreme <- sort(table(genes.extreme), decreasing=TRUE)
head(genes.extreme)
# Excludes genes that have a few extreme intron bins
cutoff.minnbins <- 50
table(genes.extreme>=cutoff.minnbins)
genes.extreme <- genes.extreme[genes.extreme>=cutoff.minnbins]
length(genes.extreme)
head(genes.extreme)
# Get gene information
m <- match(names(genes.extreme), genes.type$GeneID)
genes.extreme <- genes.type[m,]
m <- match(genes.extreme$GeneID, x.exon$genes$GeneID)
genes.extreme$Strand <- x.exon$genes$Strand[m]
genes.extreme
# Plot genes with the most number of extreme intron bins
pdf("coverage_for_HCC827/coverage_profile_of_IRlike_genes_in_R1_polyA.pdf", height=5)
par(mfrow=c(3,2))
for (i in 1:nrow(genes.extreme)){
  g <- genes.extreme$GeneID[i]
  name.i <- genes.extreme$Symbol[i]
  strand.i <- genes.extreme$Strand[i]
  if(strand.i=="-") {
    rev <- TRUE
  } else {
    rev <- FALSE
  }
  plotSubcounts(gene=g, exon=x.exon[,"R1.HCC827"], intron=x.intron[,"R1.HCC827"], intron.high=cutoff.logcoverage, groups="", save.to.file=FALSE, cols = c("grey", "dodgerblue", "black"), reverse=rev)
  title(xlab=name.i, line=2)
  plotSubcounts(gene=g, exon=x.exon[,"R1.HCC827_T"], intron=x.intron[,"R1.HCC827_T"], intron.high=cutoff.logcoverage, groups="", save.to.file=FALSE, cols = c("grey", "dodgerblue", "black"), reverse=rev)
  title(xlab=name.i, line=2)
}
dev.off()


# COVERAGE PROFILES FOR GENES IN RASKO PAPER
rasko <- c("LMNB1", "LBR", "PYGL", "FXYD5", "EIF1", "NRM")
m <- match(rasko, genes.type$Symbol)
rasko <- genes.type[m,]
m <- match(rasko$GeneID, x.exon$genes$GeneID)
rasko$Strand <- x.exon$genes$Strand[m]
rasko <- rasko[rasko$GeneID %in% genes,]
rasko <- rasko[!is.na(rasko$Strand),]
pdf("coverage_for_HCC827/coverage_profile_of_rasko_genes.pdf", height=7)
par(mfrow=c(4,2))
for (i in 1:nrow(rasko)){
  g <- rasko$GeneID[i]
  name.i <- rasko$Symbol[i]
  strand.i <- rasko$Strand[i]
  if(strand.i=="-") {
    rev <- TRUE
  } else {
    rev <- FALSE
  }
  plotSubcounts(gene=g, exon=x.exon[,groups=="PolyA"], intron=x.intron[,groups=="PolyA"], groups=rep("",3), save.to.file=FALSE, cols = c("grey", "dodgerblue", "gray40"), reverse=rev)
  title(xlab=name.i, line=2)
  plotSubcounts(gene=g, exon=x.exon[,groups=="Total"], intron=x.intron[,groups=="Total"], groups=rep("",3), save.to.file=FALSE, cols = c("grey", "dodgerblue", "gray40"), reverse=rev)
  title(xlab=name.i, line=2)
}
dev.off()