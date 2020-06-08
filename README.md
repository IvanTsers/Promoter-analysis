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

## Instructions
**NB: If you want to save time and just use PlantPAN3.0 PWMs as the input for MAST, do following steps:**

1. Create an empty folder on your machine (name it 'Promoter-analysis' or whatever you like).
2. Download *'Run_MAST'*, *'MAST_XML_parser'*, *'TF_family_regulons_correlation_analysis'*, and *'TF_regulons_enrichment_analysis'* folders from this repository into the folder you have created.
3. Download [PlantPAN_TF_annotation_filtered.tsv](https://mega.nz/file/eW5jEDzD#5y_PfsgiBfrVan8pgtdImu4P8byE0gH4ztkF7CNlXrE), put it into *'MAST_XML_parser'* folder.

**NB: If you're interested in how *PlantPAN_TF_annotation_filtered.tsv* and chunked *PlantPAN_meme_motifs* were produced, you may perform full analysis by yourself:**

1. Create an empty folder on your machine (name it 'Promoter-analysis' or whatever you like).
2. Download all the folders from this repository into the folder you have created.
3. Download PWMs of TF binding sites (all plants) from [PlantPAN3.0](http://plantpan.itps.ncku.edu.tw/download/home.php). Put the file *'Transcription_factor_weight_matrix.txt'* into into the folder you have created.
4. Download the ID mapping file (all plants) from [PlantPAN3.0](http://plantpan.itps.ncku.edu.tw/download/home.php). Put the file *'ID_mapping_all_plant.txt'* into into the folder you have created.
4. Use [RimGubaev's script](https://github.com/RimGubaev/extract_promoters) to extract promoters of your species' genes. Put the file *'Promoters.fa'* into *'Run_MAST'* directory.

## Example datasets weblinks
 * [Promoters.fa (56.3 MB)](https://mega.nz/file/CS4RmbxA#eF2pFr8gVK7P05XmTVp6GUJ_Ne27ERF9oT77NRe313w)
 * [MAST_output (1.34 GB)](https://mega.nz/folder/OepnWDST#2Pw3pp1t0SdNH2ckBfbWtQ)
 * [PlantPAN_meme_motifs (959 KB)](https://mega.nz/folder/zewBGZoZ#vbgjD8kxT81ah6q6YxV67A)
 * [Example_expression_table.tsv (2 MB)](https://mega.nz/file/uOhnAbKY#4mp5yTA-lLanGrGH247M_mLx-7wUEcAKslTrdxaO0u4)
