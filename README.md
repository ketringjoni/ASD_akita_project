# De novo structural variants in autism spectrum disorder disrupt distal regulatory interactions of neuronal genes

# Abstract
**Background:** Three-dimensional genome organization plays a critical role in gene regulation, and disruption of chromatin structure has been shown to lead to developmental disorders through the changed contact between key genes and their distal regulatory elements. Structural variants (SVs) have the ability to disrupt local genome organization, such as the joining of topologically associated domains upon deletion of a boundary. Testing large numbers of SVs for their effects on chromatin structure and gene expression is time and cost prohibitive. To overcome these experimental limitations, we propose to predict the effects of SVs on genome folding computationally and use the results to prioritize causal hypotheses to test in cells. 

**Results:** We developed a weighted scoring method to measure chromatin contact changes that specifically affect regions of interest, such as cell type-specific regulatory elements or promoters. We implemented this scoring approach in the SuPreMo-Akita software package and used it to  rank hundreds of de novo SVs (dnSVs) from autism spectrum disorder (ASD) individuals and their unaffected siblings based on predictions of how they alter the neuronal regulatory interactions of nearby genes. This revealed that putative cis-regulatory element interactions (CREints) are more disrupted by dnSVs from ASD probands versus unaffected siblings, allowing us to prioritize candidate ASD variants that disrupt CREints of genes involved in the disorder. We experimentally validated our top locus using isogenic induced pluripotent stem cell-derived excitatory neurons with and without the dnSV, finding that our prediction of disrupted contacts in the locus was highly accurate. RNA-seq analysis identified 1,102 mis-regulated genes in cell lines with the dnSV. This in vitro characterization of a candidate causal variant is important, because most ASD patients do not carry a damaging protein-coding variant that alters neurodevelopment or neuronal function. 

**Conclusion:** This study establishes disrupted genome folding as a genetic mechanism worthy of further study in ASD and provides a general strategy for prioritizing variants predicted to disrupt regulatory interactions in any tissue.

***

## In this repo

1. [Simons Simplex Collection dnSVs](https://github.com/ketringjoni/ASD_akita_project/tree/main?tab=readme-ov-file#1-simons-simplex-collection-dnSVs)
2. [Scoring SSC dnSVs with SuPreMo-Akita](https://github.com/ketringjoni/ASD_akita_project/tree/main?tab=readme-ov-file#2-scoring-ssc-dnsvs-with-SuPreMo-Akita)
3. [Scoring SSC snSVs with CREint weights](https://github.com/ketringjoni/ASD_akita_project/tree/main?tab=readme-ov-file#3-scoring-ssc-snsvs-with-creint-weights)
4. [Filtering SSC dnSVs using selection criteria](https://github.com/ketringjoni/ASD_akita_project/tree/main?tab=readme-ov-file#4-filtering-ssc-dnsvs-using-selection-criteria)
5. [HiC data analysis](https://github.com/ketringjoni/ASD_akita_project/tree/main?tab=readme-ov-file#5-hic-data-analysis)
6. [RNAseq data analysis](https://github.com/ketringjoni/ASD_akita_project/tree/main?tab=readme-ov-file#6-rnaseq-data-analysis)


### 1. Simons Simplex Collection dnSVs

De novo structural variants used in this study are from Simons Simplex Collection (SSC), as a part of Simons Foundation Autism Research Initiative (SFARI). Belyeu et al 2021 called dnSVs in hg38 using alignment-based, short-read WGS, and we pulled them from their Supplementary Table 1 into [data/SFARI_SSC_dnSVs.csv](https://github.com/ketringjoni/ASD_akita_project/data/SFARI_SSC_dnSVs.csv).


### 2. Scoring SSC dnSVs with SuPreMo-Akita

We installed [SuPreMo-Akita](https://github.com/ketringjoni/SuPreMo?tab=readme-ov-file#install-supremo-or-supremo-akita) following instructions on the repo. 

We scored SSC dnSVs: 
```
SuPreMo/scripts/SuPreMo.py data/SSC_dnSV_for_SuPreMo.txt \
--get_Akita_scores \
--shift_by -1 0 1 \
--revcomp add_revcomp \
--dir scoring_dnSVs/supremo-akita_output \
--file SSC_dnSV
```
Our steps to process the input and output files are in [scoring_dnSVs.ipynb under Scoring dnSVs using SuPreMo-Akita](https://github.com/ketringjoni/ASD_akita_project/scoring_dnSVs/scoring_dnSVs.ipynb#Scoring-dnSVs-using-SuPreMo-Akita).


### 3.Scoring SSC snSVs with CREint weights

CREints were processed from Song et al 2020 excitatory neuron H3K4me3 PLACseq data, pulled from [NeMO](https://assets.nemoarchive.org/dat-uioqy8b) into [data/eN.MAPS.peaks.txt](https://github.com/ketringjoni/ASD_akita_project/data/eN.MAPS.peaks.txt).

We scored SSC dnSVs near CREints with weighted scoring:
```
SuPreMo/scripts/SuPreMo.py EP_for_SuPreMo.txt \
--get_Akita_scores \
--get_maps \
--get_tracks \
--shifts_file EP_for_SuPreMo_shifts.txt \
--roi data/EP_for_SuPreMo_weights.txt \
--roi_scales 10 1000000 \
--dir data/supremo-akita_output_weighted \
--file EP
```
Our steps to process the input and output files are in [scoring_dnSVs.ipynb under Scoring dnSVs near CREints using SuPreMo-Akita with weighted scoring](https://github.com/ketringjoni/ASD_akita_project/scoring_dnSVs/scoring_dnSVs.ipynb#Scoring-dnSVs-near-CREints-using-SuPreMo-Akita-with-weighted-scoring).


### 5. Filtering SSC dnSVs using selection criteria
*TBD*


### 6. HiC data analysis

To process hic fastq files into mcool files, we use the [4DN pipeline](https://data.4dnucleome.org/resources/data-analysis/hi_c-processing-pipeline). The code in [HiCanalyses/hic_analysis.sh](https://github.com/ketringjoni/ASD_akita_project/HiCanalyses/hic_analysis.sh) has beed adapted from the [4DN HiC Docker GitHub repo](https://github.com/4dn-dcic/docker-4dn-hic/tree/master). We used pacakge versions shown in [HiCanalyses/hic_analyses.yml](https://github.com/ketringjoni/ASD_akita_project/HiCanalyses/hic_analysis.sh) and ran:

```
hic_analyses.sh <nthreads>, <genome_index>, <chrom_sizes>, <fastq1_rep1>, <fastq2_rep1>, <fastq1_rep2>, <fastq2_rep2>, <prefix_rep1>, <prefix_rep2>, <prefix>, <outdir>, <hic_analysis_path>
```

Raw and analyzed HiC data can be found in GEO (accession number TBD).


### 7. RNAseq data analysis
*TBD*

Raw and analyzed RNAseq data can be found in GEO (accession number TBD).


***
*For any questions and/or feedback, please reach out to katie.gjoni at ucsf dot edu*

