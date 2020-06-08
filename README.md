# Promoter-analysis
These scripts are to 
  1. Search cis-regulatory elements (CREs) within gene promoter regions using position weight matrices (PWMs) obtained from PlantPAN3.0 (http://plantpan.itps.ncku.edu.tw/index.html) and MAST tool from MEME suite (http://meme-suite.org/); 
  2. Combine CRE search results with differential gene expression (DGE) analysis to predict the potential master-regulators among plant transcription factors (TFs);
  3. Infer certain potential TF families responsible for the differential regulation of genes belonging to the particular multigene families within which both up- and downregulated genes were well-represented.
  
**NB**: The extraction of gene promoter regions were performed using extract_promoters.sh shell script (https://github.com/RimGubaev/extract_promoters). You can download promoter regions of Nicotiana tabacum genes to use them as an example input for *run_MAST_parallel.bash*

## System requirements:
* Multi-core CPU (for parallel computations)
* Linux OS is recommended (tested on Ubuntu 14.04)
* MEME suite is to be installed on your system (http://meme-suite.org/)
* R packages: data.table, ggplot2, ggpubr, grid, gridExtra, reshape2, XML

## Quick start
1. Create an empty folder on your machine (name it 'Promoter-analysis' or whatever you like).
2. Download all the folders from this repository into the folder you have created.
