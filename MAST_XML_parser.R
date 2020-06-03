# WARNING: this script is concieved to perform MULTI-THREADED parsing of multiple BIG XML files.
# The number of threads (in our case it's 6) is declared in the line "cl <- makeCluster(6)". 
# One can reduce or increase it according to the number of avaliable CPU cores.

library(doParallel)
library(XML)
library(rlist)

options(stringsAsFactors = F)

# This function groups matrices from individual files into j new files by d matrices maximum in each
file_packer <- function(X, d){
  dir.create('../mast_merged')
  l <- length(X)
  j <- ceiling(l/d) # j - the main iterator, which defines i
  for (i in 1:j){
    a = (i-1)*d + i # the left border of closed interval of file names
    ifelse(i!=j, (b = d*i+i), (b = l)) # the right border of closed interval of file names
    lines <- lapply(files[a:b], readLines) # read d files (numbered from a to b in list X) line by line
    write(c(header, unlist(lines)), file = paste('../mast_merged/part_', i, sep = ''))
  }
}

setwd('/home/rim/Arbeit/Ntab_promoters_fall_2019/Prom_motifs/jaspar_splitted/')

files <- list.files() #get vector with all motifs files names
header <- readLines(files[1], n=9) #these strings will be added at all new files
system("sed -i '1,9d' MA*") #invokes bash command to cut out first 9 lines from every file
file_packer(files, 50)

setwd('../')
files <- list.files('./mast_merged/')

cl <- makeCluster(6)
registerDoParallel(cl)
rm(cl)

foreach(i = 1:length(files)) %dopar% 
  {
    dir.create(paste('./MAST_out/', files[i], sep = ''))
    system(
      paste('mast ',
            getwd(), '/mast_merged/', files[i],
            ' /home/rim/References/Nicotiana_tabacum/promoters_500/promoters.fa -oc ',
            getwd(),'/MAST_out/', files[i], ' -mf ', files[i],
            sep = '')
           )
  }

#Parsing of XML files will be done by XML library for R. We just need to re-organize the data obtained by it
setwd('/home/rim/Arbeit/Ntab_promoters_fall_2019/Prom_motifs/')

#Warning! There're time-consuming operations below
#All links to the certain values in lists are made according to the structure of .xml file
#The general loop (with 'y' iterator) 

m <- motifs[names(motifs) == 'motif']
m <- foreach(z = 1:length(m), .combine = 'rbind') %do%
{
  as.data.frame(rbind(m[[z]][c(1, 3)]))
}

colnames(m) <- c('motif', 'Motif_JASPAR_ID') #table with motifs

mast_output <- foreach(y = 1:length(files), .combine = 'rbind') %dopar%
{
    result_xml <- xmlParse(paste(getwd(), '/MAST_out/1000_bp/part_', y, '/mast.xml', sep =''))
    motifs <- xmlToList(result_xml) #Relatively slow process
    sequences <- motifs$sequences #get sequences block
    motifs <- motifs$motifs #get motifs block
    rm(result_xml)

    dummy_row <- as.data.frame(t(rep(NA, 6)))
    colnames(dummy_row) <- c('pos', 'gap', 'motif', 'pvalue', 'strand', 'match')
    
    seq_table <- foreach(k = 2:length(sequences), .combine = 'rbind') %dopar% #Relatively slow process
      {
        sequence <- sequences[[k]]
        sequence <- list.flatten(sequence)
        ifelse('seg.hit' %in% names(sequence),
          {ts <- as.data.frame(rbind(t(as.data.frame(sequence[names(sequence) == 'seg.hit']))))
          ts <- cbind(ts, names = sequence$.attrs[4])},
          {ts <- cbind(dummy_row, names = sequence$.attrs[4])} # this builds a row with NAs instead 'hit'
        )
        ts
      }
    
    ms <- motifs[names(motifs) == 'motif']
    rm(motifs)
    ms <- foreach(z = 1:length(ms), .combine = 'rbind') %do%

      {
        as.data.frame(rbind(ms[[z]][c(1, 3)]))
      }
    colnames(ms) <- c('motif', 'Motif_JASPAR_ID')
    merge(seq_table, ms, by = 'motif', all.x = T)[, -6]
}


desc <- read.csv('JASPAR_core_desc.csv')
colnames(desc)[1] <- 'Motif_JASPAR_ID'
mast_output <- merge(mast_output, desc, by = 'Motif_JASPAR_ID', all.x = T)

mast_output <- mast_output[, -c(2, 4)]

#make strand' values shorter (to one letter)
mast_output[, 4] <- ifelse(mast_output[, 4] == 'forward', 'F', 'R')

#describe motif start position in relation to TSS
prom_len <- read.csv('/home/rim/References/Nicotiana_tabacum/promoters_1000/promoters.bed', sep = '\t', header = F)[, 2:4]
prom_len$Prom.length <- prom_len[, 2] - prom_len[, 1]
colnames(prom_len)[3] <- 'names'

mast_output <- merge(mast_output, prom_len[, 3:4], by = 'names', all.x = T)
mast_output$pos <- as.numeric(as.character(mast_output$pos))
mast_output$Start.pos <-  mast_output[, 3] - mast_output[, 10]
mast_output <- mast_output[, c(2, 1, 4, 11, 5, 6:9)]

colnames(mast_output)[2:9] <- c('GeneID', 'P-value', 'Start pos', 'Strand', 
                                'TFs name', 'Motif of species', 'TFs class', 'TFs family')


write.table(mast_output, 'TFs_binding_sites_NT_1000.tsv', sep = '\t', quote = F, row.names = F)


####################################################################
# WARNING: this script is concieved to perform MULTI-THREADED parsing of multiple BIG .xml files.
# The number of threads (in our case it's 6) is declared in the line "cl <- makeCluster(6)". 
# One can reduce or increase it according to the number of avaliable CPU cores.

library(doParallel)
library(XML)
library(rlist)

options(stringsAsFactors = F)

# Set working directory into the source file location automatically (NB: the source file must be opened in RStudio!)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

mast_output_dirs <- list.dirs('./MAST_output', recursive = F)
mast_output_dirs <- gsub('.*/', '', mast_output_dirs)

# All links to the certain values in lists (e.g. 'motifs', 'sequences') are made in accordance to the structure of .xml file
# First, we need to check whether some .xml files are 'empty' (i.e. contain no revealed matches).

empty_xmls <- foreach(i = 1:length(mast_output_dirs), .combine = 'rbind') %do%
{
  dir.name <- mast_output_dirs[i]
  result_xml <- xmlParse(paste('MAST_output/', dir.name, '/mast.xml', sep =''))
  motifs <- xmlToList(result_xml)
  sequences <- motifs$sequences
  ifelse('sequence' %in% names(sequences), 
         {df <- data.frame(dir = dir.name, empty = F)}, 
         {df <- data.frame(dir = dir.name, empty = T)})
  df
}

empty_xmls <- empty_xmls[empty_xmls$empty == T, ]

mast_output_dirs <- subset(mast_output_dirs, !(mast_output_dirs %in% empty_xmls$dir))

cl <- makeCluster(6)
registerDoParallel(cl)
rm(cl)

# This loop extracts the most informative blocks of all .xml filesa and arranges it as single dataframe
# WARNING: parsing of multiple BIG .xml files may be slow.
mast_output <- foreach(y = 1:length(mast_output_dirs), .combine = 'rbind') %do%
{
  dir.name <- mast_output_dirs[y]
  result_xml <- xmlParse(paste('MAST_output/', dir.name, '/mast.xml', sep =''))
  motifs <- xmlToList(result_xml) # Relatively slow process
  sequences <- motifs$sequences
  motifs <- motifs$motifs
  rm(result_xml)
  
  dummy_row <- as.data.frame(t(rep(NA, 6)))
  colnames(dummy_row) <- c('pos', 'gap', 'motif', 'pvalue', 'strand', 'match')
  
  seq_table <- foreach(k = 2:length(sequences), .combine = 'rbind') %dopar% # Relatively slow process
  {
    sequence <- sequences[[k]]
    sequence <- list.flatten(sequence)
    ifelse('seg.hit' %in% names(sequence),
           {tab <- as.data.frame(rbind(t(as.data.frame(sequence[names(sequence) == 'seg.hit']))))
           tab <- cbind(tab, names = sequence$.attrs[4])},
           {tab <- cbind(dummy_row, names = sequence$.attrs[4])} # build a row with NAs instead 'hit'
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

desc <- IDs[, 1:2]
colnames(desc)[1] <- 'Motif_PlantPAN_ID'
mast_output <- merge(mast_output, desc, by = 'Motif_PlantPAN_ID', all.x = T)

mast_output <- mast_output[, -c(2, 4)]

# Make strand values shorter (to one letter)
mast_output[, 4] <- ifelse(mast_output[, 4] == 'forward', 'F', 'R')

# Describe motif start position in relation to TSS
prom_len <- read.csv('/home/rim/References/Nicotiana_tabacum/promoters_1000/promoters.bed', sep = '\t', header = F)[, 2:4]
prom_len$Prom.length <- prom_len[, 2] - prom_len[, 1]
colnames(prom_len)[3] <- 'names'

mast_output <- merge(mast_output, prom_len[, 3:4], by = 'names', all.x = T)
mast_output$pos <- as.numeric(as.character(mast_output$pos))
mast_output$Start.pos <-  mast_output[, 3] - mast_output[, 7]

mast_output <- merge(mast_output, PlantPAN_annot_filtered, by = 'TF_Gene_ID', all.x = T)

mast_output <- mast_output[!is.na(mast_output$Motif_PlantPAN_ID), ]
mast_output <- mast_output[, c(2, 10, 5, 8, 6, 1, 3, 9, 11:13)]

colnames(mast_output)[1] <- 'GeneID'

mast_output[, 3] <- as.numeric(mast_output[, 3])
mast_output[, 4] <- as.numeric(mast_output[, 4])

write.table(mast_output, 'TFs_binding_sites_NT_1000.tsv', sep = '\t', quote = F, row.names = F)
