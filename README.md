# **The project works on the Variant Calling Format (VCF) file to study FST and dXY correlation across the genome.**

## **The working directory must follow the tree structure below.**
```text
.
├── config
│   └── config.yaml
├── README.md
├── resources
│   ├── metadata
│   └── vcf_data
├── results
└── workflow
    ├── envs
    │   ├── pixy.yaml
    │   ├── plink.yaml
    │   ├── r.yaml
    │   └── vcftools_bcftools_vcflib.yaml
    ├── R_scripts
    │   ├── fst_dxy.R
    │   ├── pca.R
    │   └── quality_assess.R
    └── Snakefile
```

### *Notes:*
* The resources directory needs to be mentioned by the user with the required VCF file and metadata.
* This Snakemake works for only one VCF file.
* R scripts need to be edited according to the dataset.


