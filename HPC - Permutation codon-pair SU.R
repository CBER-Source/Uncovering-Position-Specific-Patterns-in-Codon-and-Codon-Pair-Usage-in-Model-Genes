library(foreach)
library(doParallel)

#Character vector of gene names to analyze.
genes <- c()

#Set up a parallel backend for the codon pair count function.
#Choose an integer for the number of cores to use during parallelization.
clustnum <- 
co <- clustnum
cl <- makeForkCluster(nnodes = co)
registerDoParallel(cl, clustnum)

for (g in genes) {
  
  #Set the working directory.
  pwd <- paste0()  
  setwd(pwd)
  
  #Load the saved nucleotide matrix and convert it into a codon matrix
  load(file = paste0(g, "_codon_alignment.mat"))
  cod_mat <- cod_mat[,1:(ncol(cod_mat)-1)]
  
  #Vector of random codon sampling method names.
  model_vector <- c("pos", "aln")
  
  for (loop_model in model_vector) {
    
    total_permuations <- 10000
    
    #Calculate codon-pair symmetric uncertainty for each permutation.
    permutation_SU <- foreach(p = 1:total_permuations, .combine = 'rbind', .multicombine = T, .inorder = T) %dopar% {
      
      load(file = paste0(pwd, "permutation_test/", loop_model, "/", g, "_", loop_model, "_", p, "_p.mat"))
      
        foreach(u = 1:(ncol(permutation_matrix)-1), .combine = 'c') %do% {
        
        count <-  cbind(permutation_matrix[,u], permutation_matrix[,u+1])
        
        #Remove irrelevant codon-pairs.
        bad_inds <- c(grep("-", count[,1]),
        grep("-", count[,2]),
        grep("N", count[,1]),
        grep("N", count[,2]))
  
        bad_inds <- unique(bad_inds)
         
        if (length(bad_inds) == nrow(count)) {
          NA
        } else if (length(bad_inds) == (nrow(count)-1)) {
          NA
        } else {
          
          if (length(bad_inds) > 0) {
            count <- count[-bad_inds,]
          }
          
          e1 <- infotheo::entropy(count[,1])
          e2 <- infotheo::entropy(count[,2])
          mut <- infotheo::mutinformation(count[,1], count[,2])
          
          normalize <- function(calculated_mutual_information, Entropy_X, Entropy_Y)
          {
            if (Entropy_X==0 & Entropy_Y==0) {
              1
            } else if (Entropy_X==0 | Entropy_Y==0) {
              0
            } else {
              2*calculated_mutual_information/(Entropy_X+Entropy_Y)
            }
          }
          
          normalize(mut, e1, e2)
        }
        
      }
    }
    
    write.csv(permutation_SU, file = paste0(g, "_", loop_model, "_SU.csv"))
    save(permutation_SU, file = paste0(g, "_", loop_model, ".SU"))
    
    #Calculate amino acid pair symmetric uncertainty for each permutation.
    permutation_aap_SU <- foreach(p = 1:total_permuations, .combine = 'rbind', .multicombine = T, .inorder = T) %dopar% {

      load(file = paste0(pwd, "permutation_test/", loop_model, "/", g, "_", loop_model, "_", p, "_p.mat"))
      
      permutation_matrix <- foreach(n = 1:nrow(permutation_matrix), .combine = 'rbind') %do% {
        
        seqinr::translate(unlist(strsplit(paste0(permutation_matrix[n,], collapse = ""), split = character(0))))
        
      }
      
        foreach(u = 1:(ncol(permutation_matrix)-1), .combine = 'c') %do% {
        
        count <-  cbind(permutation_matrix[,u], permutation_matrix[,u+1])
        
        #Remove irrelevant amino acid pairs.
        bad_inds <- c(grep("X", count[,1]),
                      grep("X", count[,2]),
                      grep("\\*", count[,1]),
                      grep("\\*", count[,2]))
        
        bad_inds <- unique(bad_inds)
        
        if (length(bad_inds) == nrow(count)) {
          NA
        } else if (length(bad_inds) == (nrow(count)-1)) {
          NA
        } else {
          
          if (length(bad_inds) > 0) {
            count <- count[-bad_inds,]
          }
          
          e1 <- infotheo::entropy(count[,1])
          e2 <- infotheo::entropy(count[,2])
          mut <- infotheo::mutinformation(count[,1], count[,2])
          
          normalize <- function(calculated_mutual_information, Entropy_X, Entropy_Y)
          {
            if (Entropy_X==0 & Entropy_Y==0) {
              1
            } else if (Entropy_X==0 | Entropy_Y==0) {
              0
            } else {
              2*calculated_mutual_information/(Entropy_X+Entropy_Y)
            }
          }
          
          normalize(mut, e1, e2)
        }
        
      }
    }
    
    save(permutation_aap_SU, file = paste0(g, "_", loop_model, "_aap.SU"))
    
    
  }
  
}

#Close parallel backend.
stopCluster(cl)
