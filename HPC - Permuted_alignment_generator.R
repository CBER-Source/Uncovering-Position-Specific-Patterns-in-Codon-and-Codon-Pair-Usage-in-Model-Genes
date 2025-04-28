library(foreach)
library(doParallel)
library(doRNG)

#Set up a parallel backend for the codon pair count function.
#Choose an integer for the number of cores to use during parallelization.
clustnum <- 
co <- clustnum
cl <- makeForkCluster(nnodes = co)
registerDoParallel(cl, clustnum)
registerDoRNG(seed = 2021, once = TRUE)

#Character vector of gene names to analyze.
genes <- c()

for (g in genes) {
  
  #Set the working directory.
  pwd <- paste0()    
  setwd(pwd)
  
  #Load the saved nucleotide matrix and convert it into a codon matrix.
  load(paste0(g, "_c.mat"))
  cod_mat <- foreach(n = 1:(ncol(cod_mat)/3), .combine = 'cbind') %do% {
    
    paste0(cod_mat[,1+3*(n-1)], cod_mat[,2+3*(n-1)], cod_mat[,3+3*(n-1)])
    
  }
  
  cod_mat <- cod_mat[,1:(ncol(cod_mat)-1)]
  
  #aln/pos matrix permutation function
  codon_permutation <- function(codon_matrix, model_type){
    
    if (model_type == "pos") {
      
      amino_acid_matrix <- foreach(r = 1:nrow(codon_matrix), .combine = 'rbind') %do% {
        
        seqinr::translate(unlist(strsplit(codon_matrix[r,], split = character(0))))
        
      }
      
      #loop through the entire matrix and create a permutation of each column of codons
      #the results are then row bound to generate a matrix
      permuted_codon_matrix <- matrix(data = NA, nrow = nrow(codon_matrix), ncol = ncol(codon_matrix))
      
      foreach(k = 1:ncol(codon_matrix), .multicombine = T, .inorder = TRUE) %do% {
        
        foreach(u = unique(amino_acid_matrix[,k])) %do% {
          codon_seq <- codon_matrix[,k]
          aa_seq <- amino_acid_matrix[,k]
          synonymous_codon_positions <- which(aa_seq==u)
          shuffled_codons <- sample(x = codon_seq[synonymous_codon_positions], size = length(synonymous_codon_positions), replace = FALSE)
          permuted_codon_matrix[synonymous_codon_positions, k] <- shuffled_codons
          
        }
        
      }
      
    } else if (model_type == "aln") {
      
      amino_acid_matrix <- foreach(r = 1:nrow(codon_matrix), .combine = 'rbind') %do% {
        
        seqinr::translate(unlist(strsplit(codon_matrix[r,], split = character(0))))
        
      }
      
      #loop through the entire matrix and create a permutation of each column of codons
      #the results are then row bound to generate a matrix
      permuted_codon_matrix <- matrix(data = NA, nrow = nrow(codon_matrix), ncol = ncol(codon_matrix))
        
      foreach(k = 1:nrow(codon_matrix), .multicombine = T, .inorder = TRUE) %do% {
        
        foreach(u = unique(amino_acid_matrix[k,])) %do% {
          codon_seq <- codon_matrix[k,]
          aa_seq <- amino_acid_matrix[k,]
          synonymous_codon_positions <- which(aa_seq==u)
          shuffled_codons <- sample(x = codon_seq[synonymous_codon_positions], size = length(synonymous_codon_positions), replace = FALSE)
          permuted_codon_matrix[k,synonymous_codon_positions] <- shuffled_codons
          
        }
        
      }
      
    } else {
      #error detection
      stop("Unsupported Model Type")
    }
    
    return(permuted_codon_matrix)
  }
  
  #Count codon pairs in a matrix
  #The matrix elements must be codons

  #Only count alignments where the stop codon has been removed with this function.
  sequential_codon_pair_count <- function(codon_matrix) {
    
    #Index over each codon pair and count the codons at the paired location
    cp_counts <- foreach(h = 1:(ncol(codon_matrix)-1), .combine = 'rbind', .multicombine = T, .inorder = TRUE) %do% {
      
      plyr::count(cbind(codon_matrix[,h], codon_matrix[,h+1]))
      
    }
    
    #Some codon pairs will be repeated in the counts performed in the last step.
    #The names of all codon pairs are extracted here so that the results for codon pairs with the same name can be combined
    cp_names <- foreach(y = 1:nrow(cp_counts), .combine = 'c', .multicombine = T, .inorder = TRUE) %do% {
      
      paste0(cp_counts[y, 1:2], collapse = "")
      
    }
    
    #Reformat the codon pair count matrix into a vector.
    cp_counts <- cp_counts[,3]
    names(cp_counts) <- cp_names
    
    #Combine the codon pair counts by name.
    cp_counts <- tapply(cp_counts, names(cp_counts), sum)
    
    return(cp_counts)
  }
  
  model_vector <- c("pos", "aln")

  for (loop_model in model_vector) {

    total_permuations <- 10000
    
    foreach(p = 1:total_permuations) %dopar% {
      
      #Create a permuted matrix
      permutation_matrix <- codon_permutation(cod_mat, loop_model)
      
      #Count the permuted matrix's codon pair usage
      permutation_cp_count <- sequential_codon_pair_count(permutation_matrix)
      
      save(permutation_matrix, file = paste0(pwd, "permutation_test/", loop_model, "/", g, "_", loop_model, "_", p, "_p.mat"))
      save(permutation_cp_count, file = paste0(pwd, "permutation_test/", loop_model, "/", g, "_", loop_model, "_", p, "_p_cp.counts"))
      
    }

  }
  
}

#Close parallel backend.
stopCluster(cl)
