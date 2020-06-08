# Promoter-analysis: revealing the potential role of transcription factors in a given physiological process using RNA-Seq data
These scripts are to 
  1. Search cis-regulatory elements (CREs) within gene promoter regions using position weight matrices (PWMs) obtained from PlantPAN3.0 (http://plantpan.itps.ncku.edu.tw/index.html) and MAST tool from MEME suite (http://meme-suite.org/). However, you can adapt the scripts to work with other databases and motif search tools.
  2. Combine CRE search results with differential gene expression (DGE) analysis to predict the potential master-regulators among plant transcription factors (TFs);
  3. Infer certain potential TF families responsible for the differential regulation of genes belonging to the particular multigene families within which both up- and downregulated genes were well-represented.

## System requirements:
* Multi-core CPU (for parallel computations)
* Linux OS is recommended (tested on Ubuntu 14.04)
* MEME suite is to be installed on your system (http://meme-suite.org/)
* R studio (is obligatory for automatic *'setwd()'*)
* R packages: data.table, ggplot2, ggpubr, grid, gridExtra, reshape2, XML

## Instructions
**NB: If you want to save time and just use PlantPAN3.0 PWMs as the input for MAST, do following steps:**

1. Create an empty folder on your machine (name it 'Promoter-analysis' or whatever you like).
2. Download *'Run_MAST'*, *'MAST_XML_parser'*, *'TF_family_regulons_correlation_analysis'*, and *'TF_regulons_enrichment_analysis'* folders from this repository and put them into the folder you have created.
2. Download the ID mapping file (all plants) from [PlantPAN3.0](http://plantpan.itps.ncku.edu.tw/download/home.php). Put the file *'ID_mapping_all_plant.txt'* into into the folder you have created. [Alternative direct link to ID_mapping_all_plant.txt (2.8 MB)](https://mega.nz/file/zbhRiRKC#z9KUmrPrsJmxkAyZhvaZ2JDO5rMO-70mG0a8AotnGvk)
3. Use [RimGubaev's script](https://github.com/RimGubaev/extract_promoters) to extract promoters of your species' genes. Put the output file *'Promoters.fa'* into *'Run_MAST'* directory. Example output: [Promoters.fa (56.3 MB)](https://mega.nz/file/CS4RmbxA#eF2pFr8gVK7P05XmTVp6GUJ_Ne27ERF9oT77NRe313w)
4. Download [PlantPAN_TF_annotation_filtered.tsv](https://mega.nz/file/eW5jEDzD#5y_PfsgiBfrVan8pgtdImu4P8byE0gH4ztkF7CNlXrE), put it into *'MAST_XML_parser'* folder.
5. Download [PlantPAN_meme_motifs (959 KB)](https://mega.nz/folder/zewBGZoZ#vbgjD8kxT81ah6q6YxV67A), put it into *'Run_MAST'* folder.
6. Download [PlantPAN_TF_annotation_filtered.tsv](https://mega.nz/file/eW5jEDzD#5y_PfsgiBfrVan8pgtdImu4P8byE0gH4ztkF7CNlXrE), put it into *'Run_MAST'*  folder.
7. Using bash shell, change current directory to *'Run_MAST'* ($cd full_path_to_the_folder_created_in_step_1/Run_MAST)
8. Run *'run_MAST_parallel.sh'* ($bash run_MAST_parallel.bash). The output folder (*'MAST_output'*) will appear in the current directory. Example output: [MAST_output (1.34 GB)](https://mega.nz/folder/OepnWDST#2Pw3pp1t0SdNH2ckBfbWtQ).
9. Open *'MAST_XML_parser.R'* (is located in *'MAST_XML_parser'*) in R Studio and run this script. The output file (*'mast_output_full.tsv'*) will appear in *'MAST_XML_parser'* directory. Example output: [mast_output_full.tsv (81.4 MB)](https://mega.nz/file/LLozGZyR#R0283KJ7J4s6_PmGbRPsPo0l_gDQlWrz5uv8Pi35ESI).
10. Open *'Annotate_MAST_output_full.R'* (is located in *'MAST_XML_parser'*) in R Studio and run this script. The output file (*'tf_analysis_input_annotated.tsv'*) will appear in *'MAST_XML_parser'* directory. Example output: [tf_analysis_input_annotated.tsv (104.1 MB)](https://mega.nz/file/LLozGZyR#R0283KJ7J4s6_PmGbRPsPo0l_gDQlWrz5uv8Pi35ESI).
11. **Master-regulators prediction:** put the table contains data on differential gene expression (DGE) into the folder you have created in the step 1. **NB:** the following columns must be in this table: GeneID (text or numeric), log2FC (numeric) (as shown below)

|   GeneID   | log2FC  |
| ---------- | ------- |
| 107809780  |  4.838  |
| 107760295  | -1.706  |

([example expression table (2 MB)](https://mega.nz/file/GewTWJbL#4mp5yTA-lLanGrGH247M_mLx-7wUEcAKslTrdxaO0u4)). Then open *'TF_regulons_enrichment_analysis.R '* (is located in *'TF_regulons_enrichment_analysis'*) in R Studio and run this script. The output file (*'DEG_enriched_regulons.tsv'*) will appear in *'TF_regulons_enrichment_analysis'* directory. Example output: [DEG_enriched_regulons.tsv (174 B)](https://mega.nz/file/HT4lDRgK#AfNMRrM9biKynge_6ymgact7Tmoik2s9j76ayxhLz7s).
12. **Prediction of TF families responsible for regulation of a certain group of genes:**

**NB: If you're interested in how *PlantPAN_TF_annotation_filtered.tsv* and chunked *PlantPAN_meme_motifs* were produced, you may perform the following steps:**

1. Create an empty folder on your machine (name it 'Promoter-analysis' or whatever you like).
2. Download all the folders from this repository into the folder you have created.
3. Download PWMs of TF binding sites (all plants) from [PlantPAN3.0](http://plantpan.itps.ncku.edu.tw/download/home.php). Put the file *'Transcription_factor_weight_matrix.txt'* into into the folder you have created.  [Alternative direct link to Transcription_factor_weight_matrix.txt (1.1 MB)](https://mega.nz/file/SD5HEbwR#0m7Buo6wWJPaxFHsU7qYlRj4UYI4iZCR5fVJzAH8TVk)
4. Download the ID mapping file (all plants) from [PlantPAN3.0](http://plantpan.itps.ncku.edu.tw/download/home.php). Put the file *'ID_mapping_all_plant.txt'* into into the folder you have created. [Alternative direct link to ID_mapping_all_plant.txt (2.8 MB)](https://mega.nz/file/zbhRiRKC#z9KUmrPrsJmxkAyZhvaZ2JDO5rMO-70mG0a8AotnGvk)
5. Use [RimGubaev's script](https://github.com/RimGubaev/extract_promoters) to extract promoters of your species' genes. Put the file *'Promoters.fa'* into *'Run_MAST'* directory.
6. Using bash shell, change current directory to *'Run_MAST'* ($cd full_path_to_folder_created_in_step_1/Run_MAST)
7. Run *'run_MAST_parallel.sh'* ($bash run_MAST_parallel.bash). The output folder (*'MAST_output'*) will appear in the current directory. Example output: [MAST_output (1.34 GB)](https://mega.nz/folder/OepnWDST#2Pw3pp1t0SdNH2ckBfbWtQ).
. To perform further analysis, go to step 7 of the previous section.

## Example datasets weblinks
 * [Promoters.fa (56.3 MB)](https://mega.nz/file/CS4RmbxA#eF2pFr8gVK7P05XmTVp6GUJ_Ne27ERF9oT77NRe313w)
 * [MAST_output (1.34 GB)](https://mega.nz/folder/OepnWDST#2Pw3pp1t0SdNH2ckBfbWtQ)
 * [PlantPAN_meme_motifs (959 KB)](https://mega.nz/folder/zewBGZoZ#vbgjD8kxT81ah6q6YxV67A)
 * [Example_expression_table.tsv (2 MB)](https://mega.nz/file/uOhnAbKY#4mp5yTA-lLanGrGH247M_mLx-7wUEcAKslTrdxaO0u4)
