# Supplementary Figures

**Supplementary figure 1** : RNA protocols are compared for gene-level exon log-CPM (left) and gene-level intron log-CPM values (right) for human cell line HCC827 R1, R2, R3, human cell line NCI-H11975 R1, R2, and R3 (in order of plots shown).

**Supplementary figure 2** : From left to right, intron versus exon log-CPM, intron versus exon log-RPKM, intron log-RPKM versus total intron length, relative coverage between exon and intron regions versus total intron length, and relative coverage versus exon log-RPKM; with poly(A) RNA (top row) and Total RNA libraries (bottom row). The plots are displayed for human cell line HCC827 R1, R2, R3, human cell line NCI-H11975 R1, R2, and R3 (in order of plots shown).



# Scripts for Annotations

**gtf_to_exon_saf.R** : Convert gtf file to simplified exon locations and save in saf format (also saves GENCODE gene information to txt).

**exon_saf_to_genebody_and_intron_saf.R** : Exon saf file is used to create genebody saf and intron saf files.

**genebody_saf_to_nonoverlapping_genebody_saf.R** : Create genebody saf file for non-overlapping genes.



# Scripts for Analyses

### Exploring intron reads

**analysis_of_intron_exploration.R** : Main data analysis file that explores basic characteristics of intron reads.

The analysis requires:
-    Annotation files (exon saf and genebody saf)
 -   Bam files
  -  Sample information (sample names, groups, single/paired-end reads)

The analysis includes:
   - Summarising reads into exon, intron and genebody counts
  -  Calculating exon and intron read percentages
  -  Creating MDS plots for exon and intron counts
   - Calculating percentage of reads with exon signal and intron signal
   - Calculating correlations for exon and intron log-counts
   


### Characteristics of counts for celllines
**analysis_of_characteristics_of_counts_for_celllines.R** : Scatter plots for characteristics between and within libraries in cellline data.



### Coverage patterns

**analysis_of_coverage_for_HCC827.R** : Binned coverage analysis for HCC827 human cell line libraries.

The analysis requires:
   - Output from main data analysis
   - Annotation files (exon saf, intron saf, and non-overlapping genes saf)
   - GENCODE gene information
   - Function to create binned annotation: **FUN_create_binned_annotation.R**
   - Function to plot coverage profiles: **FUN_plot_coverage_profile.R**
   
The analysis includes:
   - Selecting a subset of genes to analyse (protein coding and non-overlapping genes in reference chromosomes)
   - Creating binned annotation for exons and introns
   - Summarising reads into exon and intron coverage
   - Examining the distribution of exon and intron coverage
   - Examining the coverage along the genebody for exon and intron regions
   - Creating coverage profile plots for two short and two long genes
   - Creating coverage profile plots of intron retention-like genes
   - Creating coverage profile plots for genes in Rasko paper
   
   
### IRFinder analysis

Scripts in "IRFinder_analysis" folder.

- **IR_detection.txt** : IR detection
- **AC_differential_IR_analysis.txt** : Audic and Claverie method for differential IR analysis
- **GLM_differential_IR_analysis.R** : GLM method for differential IR analysis
- **filePaths.txt** : File paths used in GLM analysis
- **experiment.txt** : Experiment information used in GLM analysis
