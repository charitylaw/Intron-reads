# 27 Sep 2017 (Last updated 25 May 2018)

# Fields to define
genebody.filename <- "Genebody.txt"
exon.filename <- "Exon.txt"
intron.filename <- "Intron.txt"

# Load libraries
library(limma)
library(GenomicRanges)

# Find and remove overlapping genes using genebody annotation
genebody <- read.table(genebody.filename, header=TRUE)
genebody.gr <- GRanges(seqnames=as.character(genebody$Chr),
  ranges=IRanges(start=genebody$Start, end=genebody$End, names=genebody$GeneID))
overlaps <- findOverlaps(genebody.gr, genebody.gr)
self <- queryHits(overlaps)==subjectHits(overlaps)
overlaps <- overlaps[!self]
overlaps <- unique(queryHits(overlaps))
nonoverlapping <- setdiff(1:nrow(genebody),overlaps)
genebody.nonoverlapping <- genebody[nonoverlapping,]

# Remove non-standard chromosomes in genebody annotation
genebody.nonoverlapping <- genebody.nonoverlapping[grep("chr", genebody.nonoverlapping$Chr),]
genebody.nonoverlapping$Chr <- as.factor(as.character(genebody.nonoverlapping$Chr))

# Remove overlapping genes in exon and intron annotation
exon <- read.table(exon.filename, header=TRUE)
exon.nonoverlapping <- exon[exon$GeneID %in% genebody.nonoverlapping$GeneID,]
intron <- read.table(intron.filename, header=TRUE)
intron.nonoverlapping <- intron[intron$GeneID %in% genebody.nonoverlapping$GeneID,]

# Save annotations
write.table(genebody.nonoverlapping, file="Genebody_nonoverlapping.txt", row.names=FALSE, sep="\t")
write.table(exon.nonoverlapping, file="Exon_nonoverlapping.txt", row.names=FALSE, sep="\t")
write.table(intron.nonoverlapping, file="Intron_nonoverlapping.txt", row.names=FALSE, sep="\t")