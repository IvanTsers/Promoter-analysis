# WARNING: this script is concieved to perform MULTI-THREADED statistical analysis.
# By default, the number of threads is in accordance with the number of avaliable CPU cores (see line 9).

library(doParallel)

options(stringsAsFactors = F)

# Register the 'foreach %dopar%' backend
number_of_threads <- detectCores()-1
registerDoParallel(makeCluster(number_of_threads))

# Automatically set WD to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

input_table <- read.csv('../MAST_XML_parser/tf_analysis_input_annotated.tsv', sep = '\t')

# Processing the input dataframe: remove duplicates to count genes properly
no_duplications <- rownames(unique(input_table[, c("GeneID", "TF_name")]))
input_table <- input_table[rownames(input_table) %in% no_duplications, ]
rm(no_duplications)

input_table$nonDEG <- ifelse(is.na(input_table$log2FC), 1, 0)
input_table$DEG <- ifelse(!is.na(input_table$log2FC), 1, 0)

TF_list <- unique(input_table$TF_name)
TF_list <- TF_list[!is.na(TF_list)]

k <- sum(input_table$DEG)

enriched_regulons <- foreach(i = 1:length(TF_list), .combine = 'rbind') %dopar%
{
  predicted_regulon <- subset(input_table, input_table$TF_name == TF_list[i])
  
  x <- sum(predicted_regulon$DEG) # the number of DEGs in the given regulon
  m <- sum(predicted_regulon$DEG) + sum(predicted_regulon$nonDEG) # total genes in the given regulon
  n <- sum(subset(input_table, input_table$TF_name != TF_list[i])[, grepl('^DEG$|^nonDEG$', colnames(input_table))]) # total genes NOT in the given regulon
  pval <- round(phyper(x-1, m, n, k, lower.tail = F), digits = 3)
  df = data.frame(TF_name = TF_list[i], Total_DEG = x,
                  Total_nonDEG = sum(predicted_regulon$nonDEG), P_value = pval)
  df
}

enriched_regulons$FDR <- p.adjust(enriched_regulons$P_value, method = 'BH')
enriched_regulons <- enriched_regulons[enriched_regulons$FDR < 0.05, ]

write.table(enriched_regulons, 'DEG_enriched_regulons.tsv', sep = '\t', quote = F, row.names = F)
