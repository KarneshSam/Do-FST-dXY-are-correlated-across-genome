## **Genome-wide Patterns of Divergence and Differentiation Among Populations**

### **Overview**
The project performs a population genomic analysis to investigate differentiation and divergence among multiple populations using variant data.
The workflow includes variant quality assessment, filtering, population structure inference using Principal Component Analysis (PCA), and estimation of genetic differentiation metrics (FST) and absolute divergence (dXY). Finally, understanding how FST and dXY were correlated across the genome between different populations.

The pipeline takes raw variant calling format (VCF) as input, process them and filter the data using vcftools. Then the PLINK is used to perform principal component analysis (PCA) to understand the population structure and finally, the dXY and FST were estimated using pixy.

### **Purpose of the Analysis**
Understanding genetic differentiation and divergence among populations is central to evolutionary biology, conservation genetics, and population genomics. This project aims to:

* Assess the quality and characteristics of raw variant data
* Generate a high-confidence bi-allelic SNPs by considering filtering thresholds
* Identify population structure using PCA
* Quantify FST and dXY between populations
* Correlation test (spearman) between FST and dXY across genome

This type of analysis is commonly used in:

* Population genetics
* Understand the evolutionary changes
* Population structure
* Selection and Gene flow

### **Data Description**
#### *Input Data*
The input data consist of a compressed Variant Call Format (VCF) file containing 3,816,977 variant sites from 16 individuals/subpopulations:

* 8N..
* K00..
* Lesina..
```
No of variant sites:
bcftools view -H ProjTaxa.vcf | wc -l

No of individuals/subpopulation:
bcftools query -l ProjTaxa.over.filter.recode.vcf | wc -l 
```

The VCF file includes only polymorphic sites (no invariant sites) and spans chromosomes 5 and Z.

Population assignments are provided in a separate metadata file.

### **Workflow**
The pipeline consists of the following major steps:

*1. Random Subsampling:*

A small proportion of variants is randomly sampled from the raw VCF to perform exploratory quality assessment.

*2. Variant Quality Assessment:*

The subset is analyzed to evaluate:
* Allele frequency distribution
* Mean sequencing depth per site
* Site quality scores
* Missing data rates
* Variant quality

This step informs filtering thresholds for downstream analyses.

*3. Variant Filtering:*

The full dataset is filtered to retain high-quality biallelic SNPs based on criteria such as:
* Removal of indels
* Missingness
* Site Quality
* Mean depth per site and per individual
* removal of maf

| Filtering criteria              | Threshold        |
|---------------------------------|------------------|
| Minor allele frequency          | > 0.1            |
| Missing data                    | ≤ 0.1 missing    |
| Variant quality                 | > 30             |
| Minimum mean depth (site)       | > 5              |
| Maximum mean depth (site)       | < 30             |
| Individual minimum depth        | > 5              |
| Individual maximum depth        | < 30             |

The above quality measures and filtering were done using vcftools [(O’Leary et al., 2018)](https://pubmed.ncbi.nlm.nih.gov/29987880/) and script was generated using the following link; [(Filtering of variants)](https://speciationgenomics.github.io/filtering_vcfs/) 

*4. Population Structure Analysis (PCA):*

Filtered variants are LD-pruned and used to compute principal components. PCA reveals clustering patterns corresponding to population structure.

The PCA analysis was done using PLINK software [(Sobota et al., 2015)](https://pubmed.ncbi.nlm.nih.gov/25644736/) and script were obtained from the following link; [(PCA)](https://speciationgenomics.github.io/pca/)

*5. Genetic Differentiation and Divergence:*

Genetic differentiation and divergence between population pairs is estimated using sliding-window analyses:
* FST - relative genetic differentiation
* dXY - absolute sequence divergence

These metrics allow identification of genomic regions with elevated differentiation.

The FST and dXY were estimated from pixy software [(Korunes & Samuk, 2021)](https://pubmed.ncbi.nlm.nih.gov/33453139/); The script were generated using the following link; [(FST and dXY)](https://stephenrdoyle.github.io/ancient_trichuris/03_code/ancient_trichuris.13_genomewide_genetic_variation.html)

*6. Visualization:*

Plots like barplot, scatterplot were generated to visualise the correlation between FST and dXY for each population pair across genome to interpret evolutionary processes such as selection, gene flow, and divergence.

### **Required Software**
This workflow uses conda environments for reproducibility.

Main tools:

* Snakemake (9.16.0)
* vcftools (0.1.17)
* PLINK (1.90b6.21)
* bcftools (1.11)
* pixy (2.0.0.beta14)
* htslib (1.23)
* vcflib (1.0.14)
* R-base (4.3.3)

All dependencies are managed through Conda environments. Each tools were specified under the worflow/envs/

### **How to Run the Workflow**
*1. Install Requiremnets:*

Install Conda (Miniconda or Mamba) and Snakemake.

*2. Configure Input Files:*

Edit the configuration file to specify:
* Path to the input VCF file
* Population metadata file

*3. Run the Pipeline:*

```
To perform a dry run:
snakemake -n

To run:
snakemake -p --software-deployment-method conda -j 12
```

### **Output File**
Key outputs include:
* Quality statistics summaries
* Filtered VCF file
* PCA results and plots
* Window-based FST estimates
* Window-based dXY estimates
* FST vs dXY visualization plot

The output files were created under the results directory.

### **People Involved**
* **Name:** Karnesh Sampath
* **Date:** 15.02.2026

### **Notes:**
* The analysis focuses on chromosomes 5 and Z present in the dataset and is not genome-wide.
* Invariant sites were not included in the VCF; appropriate methods were used to account for this during divergence estimation.
* The metadata can be changed based on the study.








