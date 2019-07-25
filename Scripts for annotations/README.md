# Supplementary Figures

**Supplementary figure 1** : RNA protocols are compared for gene-level exon log-CPM (left) and gene-level intron log-CPM values (right) for human cell line HCC827 R1, R2, R3, human cell line NCI-H11975 R1, R2, and R3 (in order of plots shown).

**Supplementary figure 2** : From left to right, intron versus exon log-CPM, intron versus exon log-RPKM, intron log-RPKM versus total intron length, relative coverage between exon and intron regions versus total intron length, and relative coverage versus exon log-RPKM; with poly(A) RNA (top row) and Total RNA libraries (bottom row). The plots are displayed for human cell line HCC827 R1, R2, R3, human cell line NCI-H11975 R1, R2, and R3 (in order of plots shown).

**Supplementary figure 3**: Scatterplot matrix of summary values from running superintronic on the polyA HCC827 cellline. We have show the exon mean (exon_mn) the intron standard deviation (intron_sd) and number of intron bases. More details can be found in the [superintronic vignette](http://htmlpreview.github.io/?https://github.com/sa-lee/analysis-superintronic/blob/master/Rmd/01-superintronic.html)

**Supplementary figure 4**: UpSet plot showing overlap between superintronic and other intron retention methods. More details can be found in the [overlaps vignette](http://htmlpreview.github.io/?https://github.com/sa-lee/analysis-superintronic/blob/master/Rmd/03-overlaps.html)

**Supplementary figure 5**: UpSet plot showing overlap between superintronic and INdEX. More details can be found in the [overlaps vignette](http://htmlpreview.github.io/?https://github.com/sa-lee/analysis-superintronic/blob/master/Rmd/03-overlaps.html).



**Supplmentary Coverage Plots** :
These contain results of running our method (superintronic)
along with IRfinder and IsoformSwitchAnalyzeR on the poly(A) libraries. 
These are links to directories that contain coverage profiles, which are named by the method and GENCODE gene id.

1. [*superintronic*](https://github.com/sa-lee/analysis-superintronic/tree/master/img/superintronic-polyA-cov)
2. [*IRFinder*](https://github.com/sa-lee/analysis-superintronic/tree/master/img/irfinder-cov) 
3. [*IsoformSwitchAnalyzeR*](https://github.com/sa-lee/analysis-superintronic/tree/master/img/isa-cov)




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

This analysis is performed in [**make_fig3.R**](https://github.com/sa-lee/analysis-superintronic/blob/0a0eca18878ffa62c9b88e951724d3b8925c1054/scripts/make_fig3.R) 
as part of [superintronic-analysis](https://github.com/sa-lee/analysis-superintronic)

The analysis requires:
   - Coverage GRanges (from scripts/run_coverage.R) 
   - GENCODE GTF
   
The analysis includes:
   - Selecting a subset of genes to analyse (protein coding and non-overlapping genes in reference chromosomes)
   - Creating binned annotation for exons and introns
   - Summarising reads into exon and intron coverage
   - Examining the distribution of exon and intron coverage
   - Examining the coverage along the genebody for exon and intron regions
   - Creating coverage profile plots for two short and two long genes
   
   
### IRFinder analysis

Scripts in "IRFinder_analysis" folder.

- **IR_detection.txt** : IR detection
- **AC_differential_IR_analysis.txt** : Audic and Claverie method for differential IR analysis
- **GLM_differential_IR_analysis.R** : GLM method for differential IR analysis
- **filePaths.txt** : File paths used in GLM analysis
- **experiment.txt** : Experiment information used in GLM analysis

### index analysis

**index_analyses.R** : index analysis of human cell line and immune cells libraries

The analysis requires:
   - DGE objects for introns, exons and genebody for cell lines and immune cells

The analysis includes:
   - Creating barplot for index categories
   - Creating scatter plots for t-statistics and logFC comparisons between introns and exons
   - Creating boxplot distribution of intron lengths for each index category

## superintronic analysis

This is a submodule to [sa-lee/analysis-superintronic](https://github.com/sa-lee/analysis-superintronic) which contains all scripts for elucidating coverage profiles for intron retention signal. It also includes cached data to reproduce the analysis (note this requires git-lfs to be installed).


The analysis requires:

- IRfinder analysis to be run
- BAM files to be processed for running coverage
- fastq files available if running kallisto and isoformSwitchAnalyzeR

This analysis includes:

- Creating coverage profiles for all IR methods tested in the paper
- Finding IR signal using superintronic for all cellline x kit combintaions
- Comparing overlaps between all methods
- Running isoformSwitchAnalyzeR

