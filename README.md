# Promoter-analysis
These scripts are to 
  1. Search cis-regulatory elements (CREs) within gene promoter regions using position weight matrices (PWMs) obtained from PlantPAN3.0 (http://plantpan.itps.ncku.edu.tw/index.html) and MAST tool from MEME suite (http://meme-suite.org/). However, you can adapt the scripts to work with other databases and motif search tools.
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

### Example datasets weblinks
 * [Promoters.fa](https://mega.nz/file/CS4RmbxA#eF2pFr8gVK7P05XmTVp6GUJ_Ne27ERF9oT77NRe313w)
 * [MAST_output](https://mega.nz/folder/OepnWDST#2Pw3pp1t0SdNH2ckBfbWtQ)
 * [PlantPAN_meme_motifs](https://mega.nz/folder/zewBGZoZ#vbgjD8kxT81ah6q6YxV67A)
 * [Example_expression_table.tsv](https://mega.nz/file/uOhnAbKY#4mp5yTA-lLanGrGH247M_mLx-7wUEcAKslTrdxaO0u4)
