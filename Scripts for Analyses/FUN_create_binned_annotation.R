# 12 Sep 2017 (Last updated 27 Sep 2017)
# Create subcounts for single and multi genes
# Charity Law

# FUNCTION: CHOPS ANNOTATION FOR A SINGLE GENE INTO SMALL BIN SIZES
createSubannoForSingleGene <- function(anno=NULL, binsize=20, featurename="FEATURE", genename=anno$GeneID[1], chr=anno$Chr[1], strand=anno$Strand[1]){
  # Cleaning up
  if(any(is.null(anno$Start), is.null(anno$End), is.null(anno$GeneID), is.null(anno$Chr), is.null(anno$Strand))) stop("Some or all slots in anno are missing. \n")
  o <- order(anno$Start)
  anno <- anno[o,]
  length <- anno$End-anno$Start
  # Create subanno object
  nexons <- nrow(anno)
  nsubexons <- sum(ceiling((length+1)/binsize))
  out <- as.data.frame(matrix(NA, nrow=nsubexons, ncol=ncol(anno)))
  colnames(out) <- colnames(anno)
  features <- paste(genename, featurename, 1:nexons, sep="_")
  out$GeneID <- rep(features, ceiling((length+1)/binsize))
  out$Chr <- chr
  out$Strand <- strand
  for(i in 1:nexons){
    i.feature <- features[i]
    pos.feature <- which(out$GeneID==i.feature)
    start <- anno$Start[i]
    start.add <- floor(length[i]/binsize)
    if(start.add>0){
      start <- c(start, start+(1:start.add)*binsize)  
    }
    if(start.add==0){
      end <- anno$End[i]
    } else {
      end <- start[-1]-1
      end <- c(end, anno$End[i])
    }
    if(length(pos.feature)!=length(start)) stop("pos.feature length not the same as start length")
    if(length(pos.feature)!=length(end)) stop("pos.feature length not the same as end length")
    out$Start[pos.feature] <- start
    out$End[pos.feature] <- end
  }
  out$Length <- out$End-out$Start
  out
}



# FUNCTION: CHOPS ANNOTATION FOR MULTIPLE GENES INTO SMALL BIN SIZES
createSubanno <- function(anno=NULL, binsize=20, featurename="FEATURE"){
  if(any(is.null(anno$Start), is.null(anno$End), is.null(anno$GeneID), is.null(anno$Chr), is.null(anno$Strand))) stop("Some or all slots in anno are missing. \n")
  genes <- unique(anno$GeneID)
  genes <- genes[!is.na(genes)]
  out <- NULL
  for(g in genes){
    anno.g <- anno[anno$GeneID==g,]
    anno.g.new <- createSubannoForSingleGene(anno.g, featurename=featurename, binsize=binsize)
    out <- rbind(out, anno.g.new)
  }
  out
}