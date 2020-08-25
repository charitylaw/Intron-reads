# Supplementary Figures

**Sup Fig 1 - Between library comparisons.pdf** : Comparing RNA protocols for gene-level exon log-CPM (left) and gene-level intron log-CPM values (right) for human cell line HCC827 R1, R2, R3, human cell line NCI-H11975 R1, R2, and R3 (in order of plots shown).

**Sup Fig 2 - Intron length.pdf** On a log2-scale, the relative length of introns compared to exons are plotted against the length of introns. The exon length of a gene is calculated as the sum of all exon regions as defined in the exon annotation; similarly for intron length. Log-values are calculated using an offset of 0.001. 

**Sup Fig 3 - Within library comparisons.pdf** : From left to right, intron versus exon log-CPM, intron versus exon log-RPKM, intron log-RPKM versus total intron length, relative coverage between exon and intron regions versus total intron length, and relative coverage versus exon log-RPKM; with poly(A) RNA (top row) and Total RNA libraries (bottom row). The plots are displayed for human cell line HCC827 R1, R2, R3, human cell line NCI-H11975 R1, R2, and R3 (in order of plots shown).

**Sup Fig 4 - Superintronic summary values.png**: Scatterplot matrix of summary values from running superintronic on the poly(A) RNA HCC827 cell line, with associated density plots for exon mean (exon_mn), intron mean (intron_mn), intron standard deviation (intron_sd) and number of intron bases above the threshold. More details can be found in the [superintronic vignette](http://htmlpreview.github.io/?https://github.com/sa-lee/analysis-superintronic/blob/master/Rmd/01-superintronic.html)

**Sup Fig 5 - Example DIR genes from IRF and ISA.pdf** Top genes with differential intron retention (DIR) between poly(A) RNA human cell lines as detected by IRFinder and IsoformSwitchAnalyzeR, visualised as coverage plots facetted by cell line. Coverage is orientated from 5' to 3', with exon regions colored green and intron regions colored orange. For IRFinder results, reported DIR regions are highlighted using a light grey rectangle on the annotation track. This is not shown for IsoformSwitchAnalyzeR as a specific region was not reported in its results. In panel (a), HNRNPL is detected by IRFinder’s generalised linear models, with its significant DIR region towards the right, or 3’ end of the gene. We selected HNRNPL since it is the only gene to contain a “clean” only flag out of the top 10 significant regions. Other genes had regions marked by a warning message. In panel (b), NBEAL2 is detected by IRFinder’s Audic and Claverie test, with its significant DIR region slightly right-of-centre of the gene. NBEAL2 contains the highest ranked DIR region (p-value=0.0044).  In panel (c), HLA-B is detected by IsoformSwitchAnalyzeR’s DEXSeq test. HLA-B is the highest ranked gene with DIR (q-value=0). Complete results tables are available at https://github.com/sa-lee/analysis-superintronic/tree/master/data .

<!--
**Overlap in intron retention methods.png**: UpSet plot showing overlap between superintronic and other intron retention methods. More details can be found in the [overlaps vignette](http://htmlpreview.github.io/?https://github.com/sa-lee/analysis-superintronic/blob/master/Rmd/03-overlaps.html)

**Overlap in novel methods.png**: UpSet plot showing overlap between superintronic and index. More details can be found in the [overlaps vignette](http://htmlpreview.github.io/?https://github.com/sa-lee/analysis-superintronic/blob/master/Rmd/03-overlaps.html).

**Coverage Plots** : Directory links for results from superintronic, IRfinder and IsoformSwitchAnalyzeR on poly(A) RNA libraries. Genes are named by GENCODE gene ID.

1. [*superintronic*](https://github.com/sa-lee/analysis-superintronic/tree/master/img/superintronic-polyA-cov)
2. [*IRFinder*](https://github.com/sa-lee/analysis-superintronic/tree/master/img/irfinder-cov) 
3. [*IsoformSwitchAnalyzeR*](https://github.com/sa-lee/analysis-superintronic/tree/master/img/isa-cov)
-->
