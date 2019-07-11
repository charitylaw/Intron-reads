# 25 Jan 2018 (Last updated 19 April 2018)
# Function for plotting binned log-coverage in genes

plotSubcounts <- function(
  gene=NULL, exon=NULL, intron=NULL, intron.high=Inf, groups=NULL, select.samples=NULL, cols=c("grey", "dodgerblue", "red"), ylim=NULL,
  loaddata=FALSE, filename=NULL, print.genename=TRUE, print.samplename=TRUE, print.strandlabel=TRUE, reverse=FALSE, save.to.file=TRUE
  )
  {
  # Load packages
  library(limma)
  library(edgeR)
  
  # Check objects
  if(any(c(is.null(gene), is.null(exon), is.null(intron)))) stop
  if(class(exon)!="DGEList") stop("exon is not a DGEList object. \n")
  if(class(intron)!="DGEList") stop("intron is not a DGEList object. \n")
  
  # Select data for gene
  exon.g <- exon[exon$genes$GeneID==gene,]
  intron.g <- intron[intron$genes$GeneID==gene,]
  
  # Combine exon and intron data
  x <- DGEList(counts=rbind(exon.g$counts, intron.g$counts))
  x$genes <- rbind(exon.g$genes, intron.g$genes)
  
  # Order subregions 
  o <- order(x$genes$Start)
  x <- x[o,]
    
  # Select samples of interest
  if(is.null(select.samples)) select.samples <- 1:ncol(x)
  if(is.null(groups)) groups <- rep(1, ncol(x))
  x <- x[,select.samples]
  groups <- groups[select.samples]
  
  # Calculate average log-coverage of groups
  logcounts <- log2(x$counts+1)
  logcounts <- t( rowsum(t(logcounts), group=groups) / rowsum(rep(1,ncol(logcounts)), group=groups)[,1] )
  if(reverse) logcounts <- logcounts[nrow(logcounts):1,,drop=FALSE]
  
  # Strand information
  strand <- x$genes$Strand[1]
  if(strand=="+") strand.label <- "5' ------ +ve strand ------ 3'"
  if(strand=="-") strand.label <- "3' ------ -ve strand ------ 5'"
  if(reverse & strand=="-") strand.label <- "5' ------ -ve strand ------ 3'"
  if(reverse & strand=="+") strand.label <- "3' ------ +ve strand ------ 5'"
  
  # Filename
  if(is.null(filename)) filename <- gene
  filename <- paste0(filename, ".pdf")
  
  # Plot
  ngroups <- length(unique(groups))
  if(save.to.file) {
    pdf(filename, height=ngroups*2)
    par(mfrow=c(ngroups,1))
  }
  par(mar=c(5,5,2,1))
  for(i in 1:ngroups){
    # Base plot
    if(is.null(ylim)) {
      barplot(logcounts[,i], col=cols[1], border=NA, xlab="", ylab="Log-coverage", xaxt="n")
    }else {
      barplot(logcounts[,i], col=cols[1], border=NA, xlab="", ylab="Log-coverage", xaxt="n", ylim=c(0,ylim))
    }
    add.zeros <- add.exons <- add.highlights <- rep(0, nrow(logcounts))
    # Add exons
    select.feature <- x$genes$FeatureType=="Exon"
    if(reverse) select.feature <- select.feature[length(select.feature):1]
    add.exons[select.feature] <- logcounts[select.feature,i]
    barplot(add.exons, col=cols[2], border=NA, add=TRUE, yaxt="n")
    # Add highlights for high expression introns
    select.feature <- x$genes$FeatureType=="Intron"
    if(reverse) select.feature <- select.feature[length(select.feature):1]
    select.feature <- select.feature & logcounts[,i]>intron.high
    add.highlights[select.feature] <- logcounts[select.feature,i]
    barplot(add.highlights, col=cols[3], border=NA, add=TRUE, yaxt="n")
    # Fix baseline color
    barplot(add.zeros, col="white", border=NA, add=TRUE, yaxt="n")
    if(print.genename) title(xlab=strand.label, line=1)
    if(print.genename) title(main=gene) 
    if(print.genename) title(xlab=colnames(logcounts)[i], line=3) 
  }
  if(save.to.file) dev.off()
}