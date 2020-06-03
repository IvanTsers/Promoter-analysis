# WARNING: this script is concieved to perform MULTI-THREADED downloads from the website.
# The number of threads (in our case it's 4) is declared in the line "cl <- makeCluster(4)". 
# One can reduce or increase it according to the number of avaliable CPU cores.

library(XML)
library(foreach)
library(doParallel)

options(stringsAsFactors = F)

# 'term' is a character string
search.term <- function(data, term) {
  ifelse(length(grep(term, data)) != 0, 
         data[which(grepl(term, data))+1], NA)
}

# Set working directory into the source file location automatically (NB: the source file must be opened in RStudio!)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

# Register the parallel backend
cl <- makeCluster(4)
registerDoParallel(cl)
rm(cl)

# Obtain the list of TF-encoding genes
geneid_list <- read.csv('./PlantPAN/ID_mapping_all_plant.txt', sep = '\t')
geneid_list <- unique(geneid_list$TF.Locus)

# Download all the annotations from the PlantPAN website
PlantPAN_annot <- foreach(g = 1:length(geneid_list), .combine = 'rbind') %do%
  {
  doc <- htmlTreeParse(paste('http://plantpan.itps.ncku.edu.tw/TF_info.php?GeneID=', geneid_list[g], sep = ''), useInternalNodes = T)
  container <- getNodeSet(doc, "//div[@class='container']")
  
  data <- capture.output(container[[3]])
  
  # Form a row for the output dataframe
  df <- data.frame(TF_Gene_ID = search.term(data, 'TF ID \\(Gene ID\\)'),
                   TF_family = search.term(data, 'TF Family'),
                   TF_name = search.term(data, 'Gene Name'),
                   Short_descr = search.term(data, 'Short Description'),
                   Species = search.term(data, 'Species'),
                   Invol_cond = search.term(data, 'Involved condition'))
  df <- lapply(df, function(x) gsub('.+">|</td.+', '', x))
  
  # Notify if the download number g is done
  print(paste('Done:', g))
  df
  }

PlantPAN_annot <- apply(PlantPAN_annot, 2, as.character)
PlantPAN_annot <- as.data.frame(PlantPAN_annot)
write.table(PlantPAN_annot, './PlantPAN/PlantPAN_annotation.tsv', sep = '\t', quote = F, row.names = F)

######################################################################
#                                                                    #
# THE FOLLOWING COMMANDS ARE TO EXTRACT ONLY THE WELL-ANNOTATED ROWS #
#                                                                    #
######################################################################

# Remove rows in which TF names are NA. It's considered to be a poor annotation
PlantPAN_annot_filtered <- PlantPAN_annot[!is.na(PlantPAN_annot$TF_name), ]

# Remove rows in which TF names are the same as GeneIDs. It's considered to be a poor annotation
PlantPAN_annot_filtered$TF_Gene_ID <- as.character(PlantPAN_annot_filtered$TF_Gene_ID)
PlantPAN_annot_filtered$TF_name <- as.character(PlantPAN_annot_filtered$TF_name)
PlantPAN_annot_filtered <- PlantPAN_annot_filtered[toupper(PlantPAN_annot_filtered$TF_name) != toupper(PlantPAN_annot_filtered$TF_Gene_ID), ]

# Remove the rest poorly annotated rows
PlantPAN_annot_filtered <- PlantPAN_annot_filtered[!grepl('XM_\\d+|<..|_POPTR|_ARALY|DRAFT_|Medtr|_BRADI', PlantPAN_annot_filtered$TF_name), ]

# Final clean-up
PlantPAN_annot_filtered$TF_name <- gsub('AT\\dG\\d.+; ', '', PlantPAN_annot_filtered$TF_name)
PlantPAN_annot_filtered <- data.frame(lapply(PlantPAN_annot_filtered, as.character), stringsAsFactors = F)

# Write the obtained well-annotated rows
write.table(PlantPAN_annot_filtered, '../PlantPAN_TF_annotation_filtered.tsv', sep = '\t', row.names = F, quote = F)
