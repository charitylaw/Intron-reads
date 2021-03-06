### Exploring intron reads

**analysis_of_intron_exploration.R** : Main data analysis file that explores basic characteristics of intron reads.

The analysis requires:
-    Annotation files (exon saf and genebody saf) from *Scripts for Annotations*
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

This analysis is performed in [**make_fig3.R**](https://github.com/sa-lee/analysis-superintronic/blob/master/scripts/make_fig3.R) 
as part of [superintronic-analysis](https://github.com/sa-lee/analysis-superintronic)

The analysis requires:
   - Coverage GRanges (from *scripts/run_coverage.R*) 
   - GENCODE GTF
   
The analysis includes:
   - Selecting a subset of genes to analyse (protein coding and non-overlapping genes in reference chromosomes)
   - Creating binned annotation for exons and introns
   - Summarising reads into exon and intron coverage
   - Examining the distribution of exon and intron coverage
   - Examining the coverage along the genebody for exon and intron regions
   - Creating coverage profile plots for two short and two long genes
   
   
   
### index analysis

**index_analyses.R** : index analysis of human cell line and immune cells libraries

The analysis requires:
   - DGE objects for introns, exons and genebody counts for cell lines and immune cells

The analysis includes:
   - Creating barplot for index categories
   - Creating scatter plots for t-statistics and logFC comparisons between introns and exons
   - Creating boxplot distribution of intron lengths for each index category


For the analysis on human cell lines, coverage profiles of the top intron+ and intron- genes (by absolute t-statistic) have been created using *superintronic* software. The plots are available at https://github.com/sa-lee/analysis-superintronic/blob/master/img/index_classifications/top_intron_class_genes . An overlap of genes reported as having nonsense mediated decay (NMD) with our results can be found at https://github.com/sa-lee/analysis-superintronic/blob/master/img/index_classifications/overlap_nmd_vs_index.png . 


### superintronic analysis

This is a submodule to [sa-lee/analysis-superintronic](https://github.com/sa-lee/analysis-superintronic) which contains all scripts for elucidating coverage profiles for intron retention signal. It also includes cached data to reproduce the analysis (note this requires git-lfs to be installed).

The analysis requires:

- IRFinder analysis to be run (see below)
- BAM files to be processed for running coverage
- fastq files available if running kallisto and isoformSwitchAnalyzeR

This analysis includes:

- Creating coverage profiles for all IR methods tested in the paper
- Finding IR signal using superintronic for all cell line x kit (RNA library protocol) combinations
- Comparing overlaps between all methods
- Running isoformSwitchAnalyzeR
- A complete write up for all parts of the analysis https://github.com/sa-lee/analysis-superintronic/tree/master/Rmd

Coverage plots for poly(A) RNA HCC827 cell lines can be found at https://github.com/sa-lee/analysis-superintronic/tree/master/img/superintronic-polyA-cov .  

Overlap between *index* and *superintronic* results for Total RNA HCC827 cell lines can be found at https://github.com/sa-lee/analysis-superintronic/blob/master/img/index_classifications/overlap_index_categories_vs_superintronic.png

Overlap between genes with intron split reads of 10 or more (over_ten) with *index* and *superintronic* results for Total RNA HCC827 cell lines can be found at https://github.com/sa-lee/analysis-superintronic/blob/master/img/index_classifications/overlap_intron_split_reads.png . 

Overlap between *superintronic* versus *IRFinder* and *IsoformSwitchAnalyzeR* is available at 
https://github.com/sa-lee/analysis-superintronic/tree/master/img/method_overlaps/overlaps_superintronic_vs_others.png

### IRFinder analysis

Scripts in "IRFinder_analysis" folder.

- **IR_detection.txt** : IR detection
- **AC_differential_IR_analysis.txt** : Audic and Claverie method for differential IR analysis
- **GLM_differential_IR_analysis.R** : GLM method for differential IR analysis
- **filePaths.txt** : File paths used in GLM analysis
- **experiment.txt** : Experiment information used in GLM analysis
