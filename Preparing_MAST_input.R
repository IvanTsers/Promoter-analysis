# WARNING: you need jaspar2meme (the tool from MEME suite) to be installed on your system.
# MEME Suite installation: http://meme-suite.org/doc/install.html?man_type=web

# The downloadable PlantPAN 3.0 position weight matrices (PWMs, hereinafter referred to as matrices) apparently lacks of uniform format. 
# Instead, they are medley of position weights and position probabilities in JASPAR-like format.
# Moreover, there are diffrent fractional number formats:
#   1) no fixed number of digits after point;
#   2) exponential format in some positions;
#   3) matrices of integers are mixed with fractions.
# Transcription_factor_weight_matrix.txt is to be formatted into .MEME via jaspar2meme tool,
# and we are to edit PlantPAN 3.0 matrices just to create an input for this tool.
# Only after such formatting one can load the matrices in MAST.

options(stringsAsFactors = F)

# This function groups matrices from individual files into j new files by d matrices maximum in each
file_sorter <- function(X, d){
  l <- length(X)
  j <- ceiling(l/d) #j - the main iterator, which defines i
  for (i in 1:j){
    folder.name <- paste('./part_', i, sep = '')
    dir.create(folder.name)
    a = (i-1)*d + i # the left border of closed interval of file names
    ifelse(i!=j, (b = d*i+i), (b = l)) # the right border of closed interval of file names
    file.copy(files[a:b], folder.name) # read d files (numbered from a to b in list X) line by line
  }
}

# Set working directory into the source file location automatically (NB: the source file must be opened in RStudio!)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

# Get downloaded PlantPAN3.0 database (PP stands for "PlantPan")
PP_matrix <- readLines('../Transcription_factor_weight_matrix.txt')

########################
#                      #
# EDITING THE MATRICES #
#                      #
########################

# Delete all numbers after 4 digits past 0. (i. e. left only 0.dddd)
PP_matrix <-gsub('(?<=\t\\[0\\.\\d{4})\\d+', '', PP_matrix, perl = T)
PP_matrix <-gsub('(?<=\t0\\.\\d{4})\\d+', '', PP_matrix, perl = T)

# Delete "0."
PP_matrix <- gsub('0\\.', '', PP_matrix)

# Delete all zeros before non-zero digits
PP_matrix <- gsub('\t0+', '\t', PP_matrix)
PP_matrix <- gsub('\t\\[0+', '\t\\[', PP_matrix)

# Restore lost "native" zeros.
PP_matrix <- gsub('\t\t', '\t0\t', PP_matrix)
PP_matrix <- gsub('\t\\[\t', '\t\\[0\t', PP_matrix)
PP_matrix <- gsub('\t\t', '\t0\t', PP_matrix)
PP_matrix <- gsub('\t\\[\t', '\t\\[0\t', PP_matrix)

# Replace all "e-NN" values with zero (it's assumed that these numbers are tend to zero)
PP_Es <- grep('e-\\d+', PP_matrix)
PP_matrix[PP_Es] <- gsub('e-\\d+', '', PP_matrix[PP_Es])
PP_matrix[PP_Es] <- gsub('\\d\\.\\d+', '0', PP_matrix[PP_Es])
PP_matrix <- gsub('\\.\\d+', '', PP_matrix)
rm(PP_Es)

###########################
#                         #
#  CONVERT THE MATRICES   #
#                         #
###########################

dir.create('Separated_matrices')
setwd('./Separated_matrices/')
           
l <- (length(PP_matrix)/5)-1

# Write all matrices as separate files

for (i in 0:l) {
  header <- PP_matrix[1+5*i]
  header <- gsub('>', '', header)
  
  matr <- PP_matrix[(1+5*i):(1+5*i+4)]
  matr <- matr[-1]
  matr <- gsub('A\t\\[|C\t\\[|G\t\\[|T\t\\[|\t\\]', '', matr)
  
  fileConn <-file(paste(header, '.pfm', sep = ''))
  writeLines(matr, fileConn)
  close(fileConn)
}
rm(fileConn, header, i, matr)

# Delete all the files containing PWMs of poorly annotated TFs

PlantPAN_annot_filtered <- read.csv('../../PlantPAN_annotation_download/Output/PlantPAN_TF_annotation_filtered.tsv', sep = '\t')
ids <- read.csv('../../ID_mapping_all_plant.txt', sep = '\t')

colnames(ids)[2] <- 'TF_Gene_ID'
ids <- ids[ids$TF_Gene_ID %in% PlantPAN_annot_filtered$TF_Gene_ID, ]
files_to_keep <- unique(ids$Matrix.ID)
files_to_keep <- paste(files_to_keep, '.pfm', sep = '')

files <- list.files()
files <- subset(files, !(files %in% files_to_keep))
file.remove(files)
rm(files_to_keep)

# Put all the matrices from the individual files in 'grouped' files containing 50 matrices.
# 50 is the optimal number of position weight matrices to be loaded in MAST.

files <- list.files()
file_sorter(files, 50)

dir.create('../PlantPAN_meme_motifs')

# You may start jaspar2meme from the shell, e.g.:
#    $ cd PATH_TO_SEPARATED_MATRICES_FOLDER
#    $ jaspar2meme -pfm ./part_1 > ../PlantPAN_meme_motifs/PlantPAN_part_1.meme
# But I suggest running the necessary commands right from here using the 'system' function.

# NB: there are troubles with TFmatrixID_0251.pfm and TFmatrixID_0434.pfm since there are positions summed to zero.
# These matrices ARE NOT represented in the further analysis.

system(paste('cd', getwd(), sep = ' '))
i <- ceiling(length(files)/50)

for (k in 1:i){
system(paste('jaspar2meme  -pfm ',
             './part_', k, 
             ' > ../PlantPAN_meme_motifs/PlantPAN_part_', k, '.meme',
             sep = ''))
}

system(paste("sed -i 's/,/./g' ./PlantPAN_meme_motifs/*.meme", sep = ''))

# Delete the big folder 'Separated_matrices'
setwd('../')
unlink('Separated_matrices/', recursive = T)

# Now one should move 'PlantPAN_meme_motifs' folder into the same directory as 'run_MAST_parallel.sh'