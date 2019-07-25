# 24 Feb 2017 (Last updated 4 June 2018)
# Charity Law and Albert Zhang

# Fields to define
exon.filename <- "Exon.txt"

# Load libraries
library(limma)

# Load exon annotation
anno <- read.delim(exon.filename)
o <- order(anno$GeneID, anno$Start)
anno <- anno[o,]

# Check for genes with non-unique chromosomes
id.chrs <- paste(anno$GeneID, anno$Chr, sep="-")
id.chrs <- unique(id.chrs)
id.chrs <- strsplit2(id.chrs, split="_")
remove.id <- unique(id.chrs[,1][duplicated(id.chrs[,1])])
remove.id

# Function for making genebody annotation
makeGenebodyAnno <- function(anno){
  # Check that anno is in the right format
  if(class(anno)!="data.frame") stop("anno is not a data.frame object. \n")
  if(any(!c("GeneID","Chr","Start","End","Strand") %in% colnames(anno))) stop("anno does not contain the correct columns/ column names. \n")
  # Order annotation by GeneID and Start
  o <- order(anno$GeneID, anno$Start)
  anno <- anno[o,]
  # Create genebody annotation by taking the min Start and max End of each gene in the exon annotation.
  start <- unlist(tapply(anno$Start, anno$GeneID, min))
  end <- unlist(tapply(anno$End, anno$GeneID, max))
  anno.new <- anno[!duplicated(anno$GeneID),]
  anno.new$Start <- as.integer(start)
  anno.new$End <- as.integer(end)
  rownames(anno.new) <- NULL
  # Return
  return(anno.new)
}

# Function for making intron annotation
makeIntronicAnno <- function(anno){
  # Check that anno is in the right format
  if(class(anno)!="data.frame") stop("anno is not a data.frame object. \n")
  if(any(!c("GeneID","Chr","Start","End","Strand") %in% colnames(anno))) stop("anno does not contain the correct columns/ column names. \n")
  # Order annotation by GeneID and Start
  o <- order(anno$GeneID, anno$Start)
  anno <- anno[o,]
  # Get number of exons for each gene
  nexons <- rowsum(x=rep(1, nrow(anno)), group=anno$GeneID, reorder=FALSE)[,1]
  cat("Out of", length(nexons), "genes,", sum(nexons==1), "contain 1 exon only and do not contain intronic regions. \n")
  # Keep genes with multiple exons
  nexons <- nexons[nexons>1]
  anno <- anno[anno$GeneID %in% names(nexons),]
  # Adjust exon positions to get intron positions
  positions <- anno[,c("Start", "End")]
  positions$Start <- positions$Start-1
  positions$End <- positions$End+1
  positions <- c(t(positions))
  # The first and last position coordinate within each gene should be removed
  last <- cumsum(nexons)*2
  first <- last+1
  first <- c(1, first[-length(first)])
  positions <- positions[-c(first, last)]
  positions <- matrix(positions, ncol=2, byrow=TRUE)
  # Create intron annotation by removing one row for each gene in the exon annotation, and inserting intron positions
  anno.new <- anno[-cumsum(nexons),]
  # Check that dimensions are correct
  if(nrow(positions)!=nrow(anno.new)) stop("Check the number of intron positions calculated. \n")
  anno.new$Start <- as.integer(positions[,1])
  anno.new$End <- as.integer(positions[,2])
  rownames(anno.new) <- NULL
  # Return
  return(anno.new)
}

# Create genebody and intron annotation
genebody <- makeGenebodyAnno(anno)
write.table(genebody, file="Genebody.txt", row.names=FALSE, sep="\t")
intron <- makeIntronicAnno(anno)
write.table(intron, file="Intron.txt", row.names=FALSE, sep="\t")
