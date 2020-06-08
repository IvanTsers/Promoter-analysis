library(data.table)

options(stringsAsFactors = F)

# Set working directory into the source file location automatically (NB: the source file must be opened in RStudio!)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Obtain the list of all PWMs
geneid_list <- read.csv('../ID_mapping_all_plant.txt', sep = '\t')[, 1:2]
colnames(geneid_list) <- c('Motif_PlantPAN_ID', 'TF_Gene_ID')

tf_annotation <- read.csv('../PlantPAN_annotation_download/Output/PlantPAN_TF_annotation_filtered.tsv', sep = '\t')
tf_annotation <- merge.data.frame(tf_annotation, geneid_list, by = 'TF_Gene_ID', all.x = T)
rm(geneid_list)

mast_output <- read.csv('mast_output_full.tsv', sep = '\t')
mast_output <- merge.data.table(mast_output, tf_annotation, by = 'Motif_PlantPAN_ID', all.x = T, allow.cartesian = T)
mast_output <- mast_output[, -c(2:4)]

# Rename the 'names' column to match GeneID column in the 'expression table' file
colnames(mast_output)[4] <- 'GeneID'
expression_table <- read.csv('../Example_expression_table.tsv', sep = '\t')

mast_output$GeneID <- as.character(mast_output$GeneID)
expression_table$GeneID <- as.character(expression_table$GeneID)
expression_plus_mast <- merge.data.table(expression_table, mast_output, by = 'GeneID', all.x = T, allow.cartesian = T)

# Now I highly recommend to replace the 'TF_name' and 'TF_family' values with ones given in the "tf_names_families_edited.tsv" file.
# The original names derived from PlantPAN site appeal to be not good enough for any kind of enrichment analysis.
# To make them suitable, TF names/families were revised (re-classified) manually:
#   1) more TF groups were itemized (e.g. AP2 family was subdivided into ERF and DREB)
#   2) synonims were deleted (e.g. 'ATMYB13;ATMYBLFGN;MYB13' entry was replaced with MYB13)

tf_names_edited <- read.csv('tf_names_families_edited.tsv', sep = '\t')
expression_plus_mast <- merge.data.table(expression_plus_mast, tf_names_edited, by = 'TF_name', all.x = T, allow.cartesian = T)
expression_plus_mast <- expression_plus_mast[, -c(1, 10, 14)]
colnames(expression_plus_mast) <- gsub('_new', '', colnames(expression_plus_mast))

# Only the part of the obtained table will be saved to file
write.table(expression_plus_mast[!is.na(expression_plus_mast$GeneID), c(1:3, 4, 12, 13)], 
            'tf_analysis_input_annotated.tsv', sep = '\t', quote = F, row.names = F)