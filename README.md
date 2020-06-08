# Promoter-analysis
  These scripts are to 
  1) search cis-regulatory elements (CREs) within gene promoter regions using position weight matrices (PWMs) obtained from PlantPAN3.0 (http://plantpan.itps.ncku.edu.tw/index.html) and MAST tool from MEME suite (http://meme-suite.org/); 
  2) combine CRE search results with differential gene expression (DGE) analysis to predict the potential master-regulators among plant transcription factors (TFs);
  3) to infer certain potential TF families responsible for the differential regulation of genes belonging to the particular multigene families within which both up- and downregulated genes were well-represented.
  
  The extraction of gene promoter regions is supposed to be performed using extract_promoters.sh shell script (https://github.com/RimGubaev/extract_promoters).

## System requirements:
Multi-core CPU (for parallel computations)
Linux OS is recommended (tested on Ubuntu 14.04)
MEME suite is to be installed on your system (http://meme-suite.org/)

R packages: 
data.table, ggplot2, ggpubr, grid, gridExtra, reshape2, XML
