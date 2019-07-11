# 22 Feb 2017 (Last updated 25 May 2018)
# Charity Law

# Fields to define
genome.filename <- "gencode.v27.annotation.gtf"

# Load libraries
library(limma)
library(GenomicRanges)
library(rtracklayer)

# Modify gtf function
modGtf <- function(gtf, gtf.new, method="fatten", remove.gene.overlap=TRUE, ignore.strand=FALSE){
  # import gtf file and remove non-"exon" features
  gtf1 <- import(gtf)
  idx <- mcols(gtf1)$type=="exon"
  gtf1 <- gtf1[idx]
  # split gtf by geneid
  gtf1.g <- split(gtf1, mcols(gtf1)$gene_id)
  # fatten or flatten within gene and combine
  if(!method %in% c("fatten", "flatten")) stop("Unknown method \n")
  if(method=="fatten") gtf2.g <- reduce(gtf1.g)
  if(method=="flatten") gtf2.g <- disjoin(gtf1.g)
  gtf2 <- unlist(gtf2.g)
  mcols(gtf2)$gene_id <- names(gtf2)	
  # remove regions that are overlapping between genes
  if(remove.gene.overlap){
    # find overlapping regions and remove them
    gtf2.disjoin <- disjoin(gtf2, ignore.strand=ignore.strand)
    co <- countOverlaps(gtf2.disjoin, gtf2)
    gtf2.disjoin <- gtf2.disjoin[co==1]	
    # reassign geneids after overlap removal
    fo <- findOverlaps(gtf2.disjoin, gtf2)
    mcols(gtf2.disjoin)$gene_id <- mcols(gtf2)$gene_id[subjectHits(fo)]
    gtf2 <- gtf2.disjoin
  }
  # add annotation for transcriptid and type
  mcols(gtf2)$transcript_id <- mcols(gtf2)$gene_id
  mcols(gtf2)$type <- "exon"
  # order new gtf by geneid
  gtf2.sort <- gtf2[order(mcols(gtf2)$gene_id, decreasing=FALSE)]
  # add annotation for exonnumber and exonid
  exonid <- split(mcols(gtf2.sort)$gene_id, mcols(gtf2.sort)$gene_id)
  exonid <- unlist(lapply(exonid, function(g){seq(1, length(g))}))
  mcols(gtf2.sort)$exon_number <- exonid
  mcols(gtf2.sort)$exon_id <- paste0(mcols(gtf2.sort)$gene_id, ":", sprintf("%03d", exonid))
  # save new gtf to file
  export(gtf2.sort, gtf.new, format = "gtf")
  invisible(gtf2.sort)
}

# Save gene information from gtf file
gtf <- read.delim(genome.filename, skip=5, header=FALSE)
genes <- gtf[gtf$V3=="gene",9]
genes <- strsplit2(genes, split=";")
genes <- as.data.frame(genes[,1:4])
colnames(genes) <- c("GeneID", "Type", "Symbol", "Level")
genes$GeneID <- gsub("gene_id ", "", genes$GeneID)
genes$Type <- gsub(" gene_type ", "", genes$Type)
genes$Symbol <- gsub(" gene_name ", "", genes$Symbol)
genes$Level <- gsub(" level ", "", genes$Level)
write.table(genes, file="Gene_information.txt", sep="\t", row.names=F)

# Simplify gtf by taking the union of overlapping exons and save as saf file 
gtf.new <- modGtf(gtf=genome.filename, gtf.new="Exon.gtf", remove.gene.overlap=FALSE)
gc.saf <- as.data.frame(cbind(
	GeneID=gtf.new$gene_id, Chr=as.character(seqnames(gtf.new)),
	Start=as.character(start(gtf.new)), End=as.character(end(gtf.new)), Strand=as.character(strand(gtf.new))
	))
o <- order(gc.saf$Chr, as.numeric(as.character(gc.saf$Start)))
gc.saf <- gc.saf[o,]
write.table(gc.saf, file="Exon.txt", sep="\t", row.names=FALSE)