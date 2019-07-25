Steps for processing annotation files:

1. **gtf_to_exon_saf.R** : Convert gtf file to simplified exon locations and save in saf format (also saves GENCODE gene information to txt).
2. **exon_saf_to_genebody_and_intron_saf.R** : Exon saf file is used to create genebody saf and intron saf files.
3. **genebody_saf_to_nonoverlapping_genebody_saf.R** : Create genebody saf file for non-overlapping genes.
