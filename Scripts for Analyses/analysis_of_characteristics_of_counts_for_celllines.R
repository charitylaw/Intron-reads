# 19 March 2019 (Last edited 15 July 2019)
# Charity Law

# Load packages
library(limma)
library(edgeR)
library(Homo.sapiens)

# Load data from Human cell lines dataset keeping "pure" samples only 
samp.total <- c(paste0("R", 1:3, ".000.Total"), paste0("R", 1:3, ".100.Total"))
samp.poly <- gsub("Total", "mRNA", samp.total)
keep <- c(samp.total, samp.poly)

# Genebody counts
load("Genebody_genelevel.RData")
x <- DGEList(counts$counts)
x$genes <- counts$annotation 
geneid <- strsplit2(rownames(x),split="\\.")[,1]
genes <- select(Homo.sapiens, keys=geneid, columns="SYMBOL", keytype="ENSEMBL")
m <- match(geneid, genes$ENSEMBL)
x$genes$Symbol <- genes$SYMBOL[m]
x <- x[,keep]
genebody <- x

# Exon counts
load("Exon_genelevel.RData")
x <- DGEList(counts$counts)
x$genes <- counts$annotation # Length is actually the length of individual exons
nfeatures <- nchar(gsub(";", "", x$genes$Strand))
x$genes <- x$genes[,c("GeneID", "Length")]
x$genes$Length <- as.numeric(as.character(x$genes$Length))
x$genes$NFeatures <- nfeatures
x <- x[,keep]
exon <- x

# Intron counts
load("Delta_genelevel.RData")
x <- DGEList(counts$counts)
x$genes <- as.data.frame(cbind(GeneID=rownames(x), Length=genebody$genes$Length-exon$genes$Length+1))
x$genes$Length <- as.numeric(as.character(x$genes$Length))
x$genes$NFeatures <- exon$genes$NFeatures-1
x <- x[,keep]
intron <- x

# Sample information
nsamples <- ncol(exon)

# Gene information 
multiexonic <- exon$genes$NFeatures>1
exon.loglength <- log2(exon$genes$Length+0.001)
intron.loglength <- log2(intron$genes$Length+0.001)
relative.loglength <- intron.loglength-exon.loglength

# Library sizes
lib.exon <- colSums(exon$counts)
lib.intron <- colSums(intron$counts)
lib.total <- lib.exon + lib.intron
genebody$samples$lib.size <- exon$samples$lib.size <- intron$samples$lib.size <- lib.total

# Log-cpm
genebody.lcpm <- cpm(genebody, log=TRUE)
exon.lcpm <- cpm(exon, log=TRUE)
intron.lcpm <- cpm(intron, log=TRUE)

# Log-rpkm
genebody.lrpkm <- rpkm(genebody, log=TRUE, length=genebody$genes$Length)
exon.lrpkm <- rpkm(exon, log=TRUE, length=exon$genes$Length)
intron.lrpkm <- rpkm(intron, log=TRUE, length=intron$genes$Length)
relative.lrpkm <- exon.lrpkm-intron.lrpkm

# Colors for plots
colpal <- c("white","lightgrey", "darkgrey", "blue","yellow","red")



# BETWEEN LIBRARIES
# LogCPM cutoff of 0.1 since total library sizes for samples are roughly 30M, 
# so if I want at least 3 counts in each gene then that is a value of roughly -3 (=log2(3/30))
# I want at least 3 counts in each gene. Library sizes are approx 30M, so count>3 is cpm>0.1, and log-CPM of roughly -3. 
pdf("between_libraries.pdf", height=3, width=6)
par(mfrow=c(1,2))
par(mar=c(4,4,1,0.1))
for (i in samp.total){
  # Plot 1
  y <- exon.lcpm[,i]
  x <- exon.lcpm[, gsub("Total", "mRNA", i)]
  keep <- y>-3 | x>-3
  remove <- !keep
  cor.val <- round(cor(x[!remove],y[!remove]), 4)
  smoothScatter(y=y[!remove], x=x[!remove], 
                main="Exon log-CPM",
                ylab="Total RNA", xlab="poly(A) RNA",
                colramp=colorRampPalette(colpal))
  legend("topleft", legend=paste0("Cor=",cor.val), bty="n", text.col="darkgreen")
  abline(0,1)
  # Plot 2
  y <- intron.lcpm[,i]
  x <- intron.lcpm[, gsub("Total", "mRNA", i)]
  keep <- y>-3 | x>-3
  remove <- !keep
  y <- y[!remove]
  x <- x[!remove]
  i.genes <- genes$Symbol[!remove]
  cor.val <- round(cor(x,y), 4)
  smoothScatter(y=y, x=x, 
                main="Intron log-CPM",
                ylab="Total RNA", xlab="poly(A) RNA",
                colramp=colorRampPalette(colpal))
  legend("topleft", legend=paste0("Cor=",cor.val), bty="n", text.col="darkgreen")
  # highlight <- i.genes %in% histone
  # points(x[highlight], y[highlight], col="green")
  abline(0,1)
}
dev.off()

# Plot of intron length
pdf("intron_length.pdf", height=3, width=3)
par(mar=c(4,4,1,0.1))
y <- relative.loglength[multiexonic]
x <- intron.loglength[multiexonic]
cor.val <- round(cor(x,y), 4)
smoothScatter(y=y, x=x, 
              main="Length",
              ylab="Log( intron length / exon length )", xlab="Intron log-length",
              colramp=colorRampPalette(colpal))
legend("topleft", legend=paste0("Cor=",cor.val), bty="n", text.col="darkgreen")
dev.off()


# WITHIN LIBRARIES of expressed genes
pdf("within_libraries.pdf", height=4, width=10)
par(mfrow=c(2,5))
par(mar=c(4,4,1,0.1))
for (i in samp.poly){
  # Plot 1a
  y <- intron.lcpm[,i]
  x <- exon.lcpm[,i]
  keep <- y>-3 & x>-3
  cor.val <- round(cor(x[keep],y[keep]), 4)
  smoothScatter(y=y[keep], x=x[keep], 
                # main="poly(A) RNA",
                ylab="Intron log-CPM", xlab="Exon log-CPM",
                colramp=colorRampPalette(colpal))
  legend("topleft", legend=paste0("Cor=",cor.val), bty="n", text.col="darkgreen")
  abline(0,1)
  # Plot 2a
  y <- intron.lrpkm[,i]
  x <- exon.lrpkm[,i]
  cor.val <- round(cor(x[keep],y[keep]), 4)
  smoothScatter(y=y[keep], x=x[keep], 
                # main="poly(A) RNA",
                ylab="Intron log-RPKM", xlab="Exon log-RPKM",
                colramp=colorRampPalette(colpal))
  legend("topleft", legend=paste0("Cor=",cor.val), bty="n", text.col="darkgreen")
  abline(0,1)
  cat(i, median(x[keep]-y[keep]), "\n")
  # Plot 3a
  y <- intron.lrpkm[,i]
  x <- intron.loglength
  cor.val <- round(cor(x[keep],y[keep]), 4)
  smoothScatter(y=y[keep], x=x[keep], 
                # main="poly(A) RNA",
                ylab="Intron log-RPKM", xlab="Intron log-length",
                colramp=colorRampPalette(colpal))
  legend("topleft", legend=paste0("Cor=",cor.val), bty="n", text.col="darkgreen")
  # Plot 4a
  y <- relative.lrpkm[,i]
  x <- intron.loglength
  cor.val <- round(cor(x[keep],y[keep]), 4)
  smoothScatter(y=y[keep], x=x[keep], 
                # main="poly(A) RNA",
                ylab="Log( exon RPKM / intron RPKM )", xlab="Intron log-length",
                colramp=colorRampPalette(colpal))
  legend("topleft", legend=paste0("Cor=",cor.val), bty="n", text.col="darkgreen")
  # Plot 5a
  y <- relative.lrpkm[,i]
  x <- exon.lrpkm[,i]
  cor.val <- round(cor(x[keep],y[keep]), 4)
  smoothScatter(y=y[keep], x=x[keep], 
                # main="poly(A) RNA",
                ylab="Log( exon RPKM / intron RPKM )", xlab="Exon log-RPKM",
                colramp=colorRampPalette(colpal))
  legend("topleft", legend=paste0("Cor=",cor.val), bty="n", text.col="darkgreen")
  
  # Plot 1b
  j <- gsub("mRNA", "Total", i)
  y <- intron.lcpm[,j]
  x <- exon.lcpm[,j]
  keep <- y>-3 & x>-3
  cor.val <- round(cor(x[keep],y[keep]), 4)
  smoothScatter(y=y[keep], x=x[keep], 
                # main="Total RNA",
                ylab="Intron log-CPM", xlab="Exon log-CPM",
                colramp=colorRampPalette(colpal))
  legend("topleft", legend=paste0("Cor=",cor.val), bty="n", text.col="darkgreen")
  abline(0,1)
  # Plot 2b
  y <- intron.lrpkm[,j]
  x <- exon.lrpkm[,j]
  cor.val <- round(cor(x[keep],y[keep]), 4)
  smoothScatter(y=y[keep], x=x[keep], 
                # main="Total RNA",
                ylab="Intron log-RPKM", xlab="Exon log-RPKM",
                colramp=colorRampPalette(colpal))
  legend("topleft", legend=paste0("Cor=",cor.val), bty="n", text.col="darkgreen")
  abline(0,1)
  cat(j, median(x[keep]-y[keep]), "\n")
  # Plot 3b
  y <- intron.lrpkm[,j]
  x <- intron.loglength
  cor.val <- round(cor(x[keep],y[keep]), 4)
  smoothScatter(y=y[keep], x=x[keep], 
                # main="Total RNA",
                ylab="Intron log-RPKM", xlab="Intron log-length",
                colramp=colorRampPalette(colpal))
  legend("topleft", legend=paste0("Cor=",cor.val), bty="n", text.col="darkgreen")
  # Plot 4b
  y <- relative.lrpkm[,j]
  x <- intron.loglength
  cor.val <- round(cor(x[keep],y[keep]), 4)
  smoothScatter(y=y[keep], x=x[keep], 
                # main="Total RNA",
                ylab="Log( exon RPKM / intron RPKM )", xlab="Intron log-length",
                colramp=colorRampPalette(colpal))
  legend("topleft", legend=paste0("Cor=",cor.val), bty="n", text.col="darkgreen")
  # Plot 5b
  y <- relative.lrpkm[,j]
  x <- exon.lrpkm[,j]
  cor.val <- round(cor(x[keep],y[keep]), 4)
  smoothScatter(y=y[keep], x=x[keep], 
                # main="Total RNA",
                ylab="Log( exon RPKM / intron RPKM )", xlab="Exon log-RPKM",
                colramp=colorRampPalette(colpal))
  legend("topleft", legend=paste0("Cor=",cor.val), bty="n", text.col="darkgreen")

}
dev.off()