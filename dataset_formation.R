library(dplyr)
library(readr)

# ExE = read.csv('Documents/jensn lab/yeast/Data File S1. Raw genetic interaction datasets: Pair-wise interaction format/SGA_ExE.txt', header = TRUE, sep = "\t")
# NxN = read.csv('Documents/jensn lab/yeast/Data File S1. Raw genetic interaction datasets: Pair-wise interaction format/SGA_NxN.txt', header = TRUE, sep = "\t")
# NxE = read.csv('Documents/jensn lab/yeast/Data File S1. Raw genetic interaction datasets: Pair-wise interaction format/SGA_ExN_NxE.txt', header = TRUE, sep = "\t")

# full_interaction_data <- read.csv('GitHub/CS598-Course-Project/full_interaction_data.csv', header = TRUE, sep = ",")

# test_genes <- read.csv('/home/dikshant/GitHub/CS598-Course-Project/GOTermGeneListForClassificationYeastTest.txt', header = TRUE, sep = '\t') 
# training_genes <- read.csv('/home/dikshant/GitHub/CS598-Course-Project/GOTermGeneListForClassificationYeastTraining.txt', header = TRUE, sep = '\t') 
# all_genes <- rbind(test_genes, training_genes)
# 
# query_genes <- test_genes$ESTID
# reference_genes <- training_genes$ESTID
# genes <- unique(all_genes$ESTID)
# 
# # data <- filter(data, grepl(paste(genes, collapse="|"), Query.Strain.ID))
# # data <- filter(data, grepl(paste(genes, collapse="|"), Array.Strain.ID))
# 
# # reduce_dataset(ExE, genes)
# single_mutant_data <- read_delim("~/Documents/jensn lab/yeast/Data File S1. Raw genetic interaction datasets: Pair-wise interaction format/strain_ids_and_single_mutant_fitness.csv", 
#                                     "\t", escape_double = FALSE, trim_ws = TRUE)

parse_fitness_helper <- function(i, j, ij, pred_tol = 0.1, class_tol = 0.1){
  pred_ij <- i*j
  pred_match <- near(pred_ij, ij, pred_ij*pred_tol)
  if (pred_match){return('independent')}
  
  i_j_match = near(i,j, i*class_tol)
  i_ij_match = near(i,ij,i*class_tol)
  j_ij_match = near(j,ij,j*class_tol)

  if (i_j_match){
    if (j_ij_match & i_ij_match){return('coequal')}
    if (j < ij){return('synergistic')}
  }
  if (i < j){
    if (j <= ij){return('masking')}
    if ((i <= ij) & (ij <= j)){return('suppressive')}
  }
  if ((ij < i) & (ij < j)){return('antagonistic')}

  return('unknown')
}

parse_fitness <- function(i, j, ij, pred_tol = 0.1, class_tol = 0.1){
  classification = parse_fitness_helper(i, j, ij, pred_tol, class_tol)
  if (classification == 'unknown'){
    classification = parse_fitness_helper(j, i, ij, pred_tol, class_tol)
  }
  return(classification)
}

double_mutant_fitness <- function(reference_gene, query_gene){
  ref_idxs <- grep(pattern = reference_gene, double_mutant_data$Query.Strain.ID)
  qer_idxs <- grep(pattern = query_gene, double_mutant_data$Array.Strain.ID)
  
  idxs <- intersect(ref_idxs, qer_idxs)
  
  return(mean(double_mutant_data$Double.mutant.fitness))
}

single_mutant_fitness <- function(gene){
  idx <- grep(gene, single_mutant_data$`Strain ID`)
  if (length(idx) > 1){
    idx <- idx[1]
  }
  return(single_mutant_data$`Single mutant fitness (26°)`[idx])
}

interaction_data_frame <- function(reference_genes, query_genes){
  reference_genes <- unique(reference_genes)
  query_genes <- unique(query_genes)
  df <- data.frame(matrix(nrow = length(query_genes), ncol = length(reference_genes)))
  rownames(df) <- query_genes
  colnames(df) <- reference_genes
  
  for (ref in reference_genes){
    ref_fitness <- single_mutant_fitness(ref)
    if (length(ref_fitness) < 1){next}
    if (is.nan(ref_fitness)){next}
    print(ref)
    for (qer in query_genes){
      qer_fitness <- single_mutant_fitness(qer)
      #print(qer_fitness)
      if (length(qer_fitness) < 1){next}
      if (is.nan(qer_fitness)){next}
      double_fitness <- double_mutant_fitness(ref, qer)
      interaction <- parse_fitness(ref_fitness, qer_fitness, double_fitness)
      df[qer,ref] <- interaction
    }
  }
  
  return(df)
}

df <- interaction_data_frame(reference_genes, query_genes)
write.csv(df, 'interaction_df.csv')