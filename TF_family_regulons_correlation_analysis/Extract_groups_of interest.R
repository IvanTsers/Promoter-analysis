options(stringsAsFactors = F)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

input_table <- read.csv('../Example_expression_table.tsv', sep = '\t')

annotation_of_interest <- c('Chitinases', 'Expansins', 'Extensins', 'XTH', 'Polygalacturonan lyases', 'Polygalacturonases', 
                            'Rhamnogalacturonate lyases', 'Ubiquitin ligases E3', 'Calmodulin-like', 'Apoptic-like cell death') 

genes_of_interest <- list()
genes_of_interest <- lapply(annotation_of_interest, function(x) input_table[input_table$Gene_group == x, ])
genes_of_interest <- do.call('rbind', genes_of_interest)

write.table(genes_of_interest, 'Genes_of_interest.tsv', sep = '\t', quote = F, row.names = F)