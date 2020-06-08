# WARNING: this script is concieved to perform MULTI-THREADED parsing of multiple BIG .xml files.
# By default, the number of threads is in accordance with the number of avaliable CPU cores (see line 11).

library(doParallel)
library(XML)
library(rlist)

options(stringsAsFactors = F)

# Register the 'foreach %dopar%' backend
number_of_threads <- detectCores()-1
registerDoParallel(makeCluster(number_of_threads))

# Set working directory into the source file location automatically (NB: the source file must be opened in RStudio!)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

mast_output_dirs <- list.dirs('../Run_MAST/MAST_output/', recursive = F)
mast_output_dirs <- gsub('.*/', '', mast_output_dirs)

# All links to the certain values in lists (e.g. 'motifs', 'sequences') are given in accordance to the structure of .xml files
# First, we need to check whether some .xml files are 'empty' (i.e. contain no revealed matches).
# WARNING: parsing of multiple BIG .xml files may be slow.

empty_xmls <- foreach(i = 1:length(mast_output_dirs), .combine = 'rbind', .packages = 'XML') %dopar%
{
  dir.name <- mast_output_dirs[i]
  result_xml <- xmlParse(paste('../Run_MAST/MAST_output/', dir.name, '/mast.xml', sep =''))
  motifs <- xmlToList(result_xml)
  sequences <- motifs$sequences
  ifelse('sequence' %in% names(sequences), 
         {df <- data.frame(dir = dir.name, empty = F)}, 
         {df <- data.frame(dir = dir.name, empty = T)})
  df
}

empty_xmls <- empty_xmls[empty_xmls$empty == T, ]

mast_output_dirs <- subset(mast_output_dirs, !(mast_output_dirs %in% empty_xmls$dir))

# This loop extracts the most informative blocks of all .xml files and arranges it as single dataframe
# WARNING: parsing of multiple BIG .xml files may be slow.

dummy_row <- as.data.frame(t(rep(NA, 6)))
colnames(dummy_row) <- c('pos', 'gap', 'motif', 'pvalue', 'strand', 'match')

mast_output <- foreach(y = 1:length(mast_output_dirs), .combine = 'rbind', .packages = c('XML', 'doParallel', 'rlist'), .verbose = T) %dopar%
{
  dir.name <- mast_output_dirs[y]
  result_xml <- xmlParse(paste('../Run_MAST/MAST_output/', dir.name, '/mast.xml', sep =''))
  motifs <- xmlToList(result_xml) # Relatively slow process
  sequences <- motifs$sequences
  motifs <- motifs$motifs
  rm(result_xml)
  
  seq_table <- foreach(k = 2:length(sequences), .combine = 'rbind') %do%
  {
    sequence <- sequences[[k]]
    sequence <- list.flatten(sequence)
    ifelse('seg.hit' %in% names(sequence),
           {tab <- as.data.frame(rbind(t(as.data.frame(sequence[names(sequence) == 'seg.hit']))))
           tab <- cbind(tab, names = sequence$.attrs[4])},
           {tab <- cbind(dummy_row, names = sequence$.attrs[4])} # build a row with NAs instead of 'hit'
    )
    tab
  }
  
  ms <- motifs[names(motifs) == 'motif']
  rm(motifs)
  ms <- foreach(z = 1:length(ms), .combine = 'rbind') %do%
    
  {
    as.data.frame(rbind(ms[[z]][c(1, 3)]))
  }
  colnames(ms) <- c('motif', 'Motif_PlantPAN_ID')
  merge(seq_table, ms, by = 'motif', all.x = T)[, -6]
}

write.table(mast_output, 'mast_output_full.tsv', sep = '\t', quote = F, row.names = F)

# Please use "Annotate_mast_output_full.R" to fuse mast_output_full.tsv and PlantPAN_TF_annotation_filtered.tsv 