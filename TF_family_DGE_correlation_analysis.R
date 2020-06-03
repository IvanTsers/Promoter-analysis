library(ggpubr)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(grid)

options(stringsAsFactors = F)

# Set working directory into the source file location automatically (NB: the source file must be opened in RStudio!)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

# Load the input dataframe containing gene groups of interest. Groups are denoted in the "Category" column
input_table <- read.csv('Groups_of_interest.tsv', sep = '\t')

# Processing the input dataframe: this need to count the number of genes in the groups
groups <- unique(input_table$Category)
input_table$Up_regulated <- ifelse(input_table$log2FC > 0, 1, 0)
input_table$Down_regulated <- ifelse(input_table$log2FC < 0, 1, 0)

# Processing the input dataframe: remove duplicates to count genes properly
no_duplications <- rownames(unique(input_table[, c("GeneID", "TF_family", 'Category')]))
input_table <- input_table[rownames(input_table) %in% no_duplications, ]
rm(no_duplications)

# Initialize the list of named dataframes (each dataframe is for certain gene group)
counts_genes <- lapply(groups, function(x)
  data.frame(TF_family = unique(input_table[input_table$Category == x, 'TF_family'])))
names(counts_genes) <- groups

# Count up- and down-regulated genes belonging to each group, then count total regulated genes in the groups
for (i in 1:length(groups)) {
  t <- counts_genes[[i]]$TF_family
  counts_genes[[i]]$Up_regulated <- sapply(t, function(x) sum(input_table[input_table$TF_family == x & input_table$Category == groups[i], 'Up_regulated']))
  counts_genes[[i]]$Down_regulated <- sapply(t, function(x) sum(input_table[input_table$TF_family == x & input_table$Category == groups[i], 'Down_regulated']))
  counts_genes[[i]]$Total <- counts_genes[[i]]$Up_regulated + counts_genes[[i]]$Down_regulated
}

# Create two folders for output plots
dir.create('Significant')
dir.create('Non-significant')

# The code below is looped. Group names serve as iterators
for (z in 1:length(groups)){
chosen_group_name <- groups[[z]]

{
# Get top 10 TF families by total counts of genes related to them  
top10 <- counts_genes[[paste(chosen_group_name)]]
top10 <- top10[order(top10$Total, decreasing = T), ][1:10, ]

# Generate output plot for DEG distribution (one with grey bars)
barplot_total_degs <- ggplot(top10[, c(1, 4)], aes(x = reorder(TF_family, -Total), y = Total)) + 
  geom_bar(stat = 'identity', width = 0.8, fill = '#A9A9A9') +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0.01, 0)) +
  xlab("\nTranscription factor family")+ylab('Number of DEGs') +
  theme_classic() +
  theme(
    axis.text.x = element_text(color = 'black', size = 11, angle = 45, hjust = 1, vjust = 0.9),
    axis.text.y = element_text(color = 'black', size = 11),
    axis.title = element_text(size = 11),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5, size = 11))

# Extract the table of certain gene group from the processed general input data frame
chosen_group <- input_table[input_table$Category == chosen_group_name, ]
chosen_group <- chosen_group[, c('GeneID', 'TF_family', 'Up_regulated', 'Down_regulated')]

# Forming the list of "+TF/-TF" dataframes for each TF family from top 10
candidates <- list()
tf_families_names <- top10$TF_family

for (i in 1:length(tf_families_names)){
  sites_presence <- list(
    chosen_group[chosen_group$TF_family == tf_families_names[i], c("GeneID", 'Up_regulated', 'Down_regulated')],
    chosen_group[chosen_group$TF_family != tf_families_names[i], c("GeneID", 'Up_regulated', 'Down_regulated')]
  )

  sites_presence[[2]] <- sites_presence[[2]][!duplicated(sites_presence[[2]]$GeneID), ]
  sites_presence[[2]] <- sites_presence[[2]][!(sites_presence[[2]]$GeneID %in% sites_presence[[1]]$GeneID), ]
  
  sites_presence <- lapply(c(1, 2), function(x) sites_presence[[x]][, 2:3])
  sites_presence <- lapply(c(1, 2), function(x) data.frame(t(colSums(sites_presence[[x]]))))

  names(sites_presence) <- c(
    paste('+', tf_families_names[i], sep = ''),
    paste('-', tf_families_names[i], sep = '')
  )

# Ranking the "+TF/-TF" groups
  for (y in 1:2){
    sites_presence[[y]]$Total <- rowSums(sites_presence[[y]])
    sites_presence[[y]]$Percent_up <- ifelse(
    sites_presence[[y]]$Total > 1,
    sites_presence[[y]]$Up_regulated/sites_presence[[y]]$Total*100, 0)
    sites_presence[[y]]$Percent_down <- ifelse(
    sites_presence[[y]]$Total > 1,
    sites_presence[[y]]$Down_regulated/sites_presence[[y]]$Total*100, 0)
  }
  ratio_1 <- 2^abs(log(sites_presence[[1]]$Percent_up/sites_presence[[1]]$Percent_down, 2))
  ratio_2 <- 2^abs(log(sites_presence[[2]]$Percent_up/sites_presence[[2]]$Percent_down, 2))
  sites_presence$Rank <- abs(ratio_1 - ratio_2)
  rm(ratio_1, ratio_2)
  
  candidates[[i]] <- sites_presence
  names(candidates)[i] <- tf_families_names[i]
}

# Prepare the data for creating the "+TF/-TF" plot (one with colored bars)
top6_candidates <- list()

# The first three groups in this list are the first three groups from the top 10 families
top6_candidates <- candidates[1:3]
candidates <- candidates[which(!(names(candidates) %in% names(top6_candidates)))]
candidates <- candidates[order(sapply(candidates, function(x) x[[3]], simplify = TRUE), decreasing = TRUE)]

# The second three groups in this list are the three top-ranked groups
top6_candidates[4:6] <- candidates[1:3]
rm(candidates)

top6_candidates_counts <- top6_candidates

for(i in 1: length(top6_candidates)){
  top6_candidates[[i]][[1]] <- top6_candidates[[i]][[1]][, -c(1:3)]
  top6_candidates[[i]][[2]] <- top6_candidates[[i]][[2]][, -c(1:3)]
  top6_candidates[[i]][[3]] <- NULL
}

for(i in 1: length(top6_candidates_counts)){
  top6_candidates_counts[[i]][[1]] <- top6_candidates_counts[[i]][[1]][, -c(3:5)]
  top6_candidates_counts[[i]][[2]] <- top6_candidates_counts[[i]][[2]][, -c(3:5)]
  top6_candidates_counts[[i]][[3]] <- NULL
}

top6_candidates <- lapply(top6_candidates, function(x) melt(x, id.vars = NULL))
top6_candidates <- do.call('rbind', top6_candidates)
colnames(top6_candidates)[3] <- 'TF_family'


top6_candidates_counts <- lapply(top6_candidates_counts, function(x) melt(x, id.vars = NULL))
top6_candidates_counts <- do.call('rbind', top6_candidates_counts)
colnames(top6_candidates_counts)[c(2, 3)] <- c('Counts', 'TF_family')

top6_candidates <- cbind(top6_candidates, Counts = top6_candidates_counts$Counts)
rm(top6_candidates_counts)

top6_candidates$TF_group <- gsub('\\+|\\-', '', top6_candidates$TF_family)

tf_families_names <- unique(top6_candidates$TF_group)

# Implement Fisher's exact test
fisher_significance <- list()

for (i in 1:6){
  M <- top6_candidates[top6_candidates$TF_group == tf_families_names[i], ]
  M <- matrix(M[1:4, 4], nrow = 2, byrow = TRUE)
  fisher_significance[[i]] <- fisher.test(M)$p.value
  names(fisher_significance)[i] <- tf_families_names[i]
}
}
rm(M)

# Arrange the obtained dataframes to create the plot
top6_candidates$TF_family <- factor(top6_candidates$TF_family, levels = unique(top6_candidates$TF_family))

# Create the "+TF/-TF" plot (one with colored bars)
barplot_top_ratios <- 
  ggplot(top6_candidates, aes(x = TF_family, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = 'stack', width = 0.9) +
  scale_fill_manual(values=c("#ff6d6d", '#82b5f0')) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 100, by = 20)) +
  scale_x_discrete(expand = c(0, 0)) +
  xlab("\nTranscription factor family")+ylab('% of DEGs') +
  geom_text(data = subset(top6_candidates, Counts != 0), aes(label = Counts), 
            position = position_stack(vjust = 0.5), color = 'black', size = 3.7) +
  theme_classic() +
  theme(
    axis.text.y = element_text(color = 'black', size = 11),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 0.9, color = 'black'),
    axis.title.x = element_text(size = 11, vjust = 4),
    axis.title.y = element_text(size = 11),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(size = 11, hjust = 0.5),
    legend.position = 'none')

# Save the .png file containing  barplot_total_degs and barplot_top_ratios to the corresponding folder:
label_significance <- paste(names(which(fisher_significance < 0.05)), collapse = ', ')
label_significance <- ifelse(label_significance == '', 'none', label_significance)

output <- ggarrange(barplot_total_degs, barplot_top_ratios, widths = c(2/5, 3/5))
output <- annotate_figure(output, top = '', 
                          fig.lab = paste('Difference is significant (p<0.05) for: ', label_significance, '\n', sep = ''),
                          fig.lab.pos = 'top.right', fig.lab.size = 8)

path_name <- ifelse(label_significance == 'none', './Non-significant', './Significant')
png_name <- paste(chosen_group_name, '.png', sep = '')
png_path <- paste(path_name, '/', png_name, sep = '')

ggsave(png_path,  output, width = 18, height = 6, units = 'cm', dpi = 350, limitsize = F)

# Notify if the barplots done
print(paste('The', chosen_group_name, 'group is processed!', sep = ' '))
}
rm(y, z, t, i)
