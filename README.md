# Cross-mappability Analysis
This repository contains analyses related to cross-mappability.

### Requirements
* Linux
* R (package: argparser, data.table, MatrixEQTL, igraph, svd, preprocessCore, limma, parallel, vioplot, corrplot)
Note: the pipeline was tested using R v3.5.1 in Linux operating system on slurm computation clusters.

### Data
Data required to run the analyses can be downloaded from [here](https://doi.org/10.6084/m9.figshare.7309625). Note: some data are not publicly sharable. But you will be able to run the pipeline even without those data. In that case, the pipeline will use intermediate results stored in the data folder.

Note:
* GTEx v7 [processed expression](https://storage.googleapis.com/gtex_analysis_v7/single_tissue_eqtl_data/GTEx_Analysis_v7_eQTL_expression_matrices.tar.gz) and [covariate data](https://storage.googleapis.com/gtex_analysis_v7/single_tissue_eqtl_data/GTEx_Analysis_v7_eQTL_covariates.tar.gz) were downloaded from the [GTEx portal](https://gtexportal.org/home/datasets)
* GTEx genotype data are available from dbGap (accession: phs000424.v7.p2)
* DGN data are available through the NIMH Center for Collaborative Genomic Studies on Mental Disorders. Instructions to aceess data are available [here](https://www.nimhgenetics.org/access_data_biomaterial.php), and inquiries should reference the “Depression Genes and Networks study (D. Levinson, PI)”.

### How to run analyses?
The main script to run all analyses is [manuscript_pipeline.sh](https://github.com/battle-lab/crossmap_analysis/blob/master/manuscript_pipeline.sh). Please configure variables in the beginning (configuration section) of the script, and execute the commands one by one on the terminal of a linux computer manually (by copying and pasting).

Note: `sh manuscript_pipeline.sh` will not work correctly, as the script often submits slurm jobs without actually running them, and some commands in the script must run after commands before them.

### Contact
Ashis Saha (ashis@jhu.edu)
