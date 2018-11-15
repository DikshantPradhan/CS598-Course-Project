library(dplyr)

# ExE = read.csv('Documents/jensn lab/yeast/Data File S1. Raw genetic interaction datasets: Pair-wise interaction format/SGA_ExE.txt', header = TRUE, sep = "\t")
# NxN = read.csv('Documents/jensn lab/yeast/Data File S1. Raw genetic interaction datasets: Pair-wise interaction format/SGA_NxN.txt', header = TRUE, sep = "\t")
# NxE = read.csv('Documents/jensn lab/yeast/Data File S1. Raw genetic interaction datasets: Pair-wise interaction format/SGA_ExN_NxE.txt', header = TRUE, sep = "\t")

full_interaction_data <- read.csv('GitHub/CS598-Course-Project/full_interaction_data.csv', header = TRUE, sep = ",")

test_genes <- read.csv('/home/dikshant/GitHub/CS598-Course-Project/GOTermGeneListForClassificationYeastTest.txt', header = TRUE, sep = '\t') 
training_genes <- read.csv('/home/dikshant/GitHub/CS598-Course-Project/GOTermGeneListForClassificationYeastTraining.txt', header = TRUE, sep = '\t') 
all_genes <- rbind(test_genes, training_genes)

genes <- unique(all_genes$ESTID)

reduce_dataset <- function(df, genes){
  query_idxs <- c()
  array_idxs <- c()
  for (gene in genes){
    print(gene)
    query_idxs <- c(query_idxs, grep(gene, df$Query.allele.name))
    array_idxs <- c(array_idxs, grep(gene, df$Array.allele.name))
  }
  
  idxs <- intersect(query_idxs, array_idxs)
  print(idxs)
}

data <- filter(data, grepl(paste(genes, collapse="|"), Query.Strain.ID))
data <- filter(data, grepl(paste(genes, collapse="|"), Array.Strain.ID))

# reduce_dataset(ExE, genes)
