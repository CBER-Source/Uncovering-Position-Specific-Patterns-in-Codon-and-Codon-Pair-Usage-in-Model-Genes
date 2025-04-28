library(foreach)
library(doParallel)
library(seqinr)

#Convert a nucleotide alignment into a codon alignment.
nt_2_codon_mat <- function(nt_matrix) {
  
  codon_mat <- foreach(n = 1:(ncol(nt_matrix)/3), .combine = 'cbind') %do% {
    
    paste0(nt_matrix[,1+3*(n-1)], nt_matrix[,2+3*(n-1)], nt_matrix[,3+3*(n-1)])
    
  }
  
  return(codon_mat)
  
}

#Convert a nucleotide alignment into an amino acid alignment.
nt_2_aa_mat <- function(nt_matrix) {
  
  aa_mat <- foreach(n = 1:nrow(nt_matrix), .combine = 'rbind') %do% {
    
    seqinr::translate(nt_matrix[n,])
    
  }
  
  return(aa_mat)
  
}

#alignment-specific amino acid pair observed to expected ratio
aln_aap_OE <- function(amino_acid_matrix) {
  
  all_paired_counts <- foreach(u = 1:(ncol(amino_acid_matrix)-1), .combine = 'rbind') %do% {
    
    AAcount <-  cbind(amino_acid_matrix[,u], amino_acid_matrix[,u+1])
    
    bad_inds <- c(grep("X", AAcount[,1]),
                  grep("X", AAcount[,2]),
                  grep("\\*", AAcount[,1]),
                  grep("\\*", AAcount[,2]))
    
    bad_inds <- unique(bad_inds)
    
    if (length(bad_inds) == nrow(AAcount)) {
      NA
    } else if (length(bad_inds) == (nrow(AAcount)-1)) {
      NA
    } else {
      
      if (length(bad_inds) > 0) {
        AAcount <- AAcount[-bad_inds,]
      }

      return(AAcount)
      
    }

  }
  
  aaf <- plyr::count(all_paired_counts[,1])
  aaf$freqn <- aaf[,2]/sum(aaf[,2])
  
  aapf <- plyr::count(all_paired_counts)

  aapf$freqn <- aapf[,3]/sum(aapf[,3])
  
  f12_f1f2 <- foreach(r = 1:nrow(aapf), .combine = 'c') %do% {
    
    aa1 <- which(aapf[r,1]==aaf[,1])
    aa2 <- which(aapf[r,2]==aaf[,1])
    
    result_val <- aapf$freqn[r]/(aaf$freqn[aa1]*aaf$freqn[aa2])
    return(result_val)
    
  }
  
  aapf$freqmn <- f12_f1f2
  aapf$g_name <- g
  
  return_dat <- aapf[,c(6,1,2,5)]
  
  return(return_dat)
  
}


#position-specific amino acid pair observed to expected ratio
aap_OE <- function(amino_acid_matrix) {
  
  OE <- foreach(u = 1:(ncol(amino_acid_matrix)-1), .combine = 'rbind') %do% {
    
    AAcount <-  cbind(amino_acid_matrix[,u], amino_acid_matrix[,u+1])
    
    bad_inds <- c(grep("X", AAcount[,1]),
                  grep("X", AAcount[,2]),
                  grep("\\*", AAcount[,1]),
                  grep("\\*", AAcount[,2]))
    
    bad_inds <- unique(bad_inds)
    
    if (length(bad_inds) == nrow(AAcount)) {
      NA
    } else if (length(bad_inds) == (nrow(AAcount)-1)) {
      NA
    } else {
      
      if (length(bad_inds) > 0) {
        AAcount <- AAcount[-bad_inds,]
      }
      
      f1 <- plyr::count(AAcount[,1])
      f1$freqn <- f1[,2]/sum(f1[,2])
      
      f2 <- plyr::count(AAcount[,2])
      f2$freqn <- f2[,2]/sum(f2[,2])
      
      f12 <- plyr::count(AAcount)
      f12$freqn <- f12[,3]/sum(f12[,3])
      
      f12_f1f2 <- foreach(r = 1:nrow(f12), .combine = 'c') %do% {
        
        aa1 <- which(f12[r,1]==f1[,1])
        aa2 <- which(f12[r,2]==f2[,1])
        
        result_val <- f12$freqn[r]/(f1$freqn[aa1]*f2$freqn[aa2])
        return(result_val)
        
      }
      
      f12$freqmn <- f12_f1f2
      f12$g_name <- g
      
      return_dat <- f12[,c(6,1,2,5)]
      
      return_dat <- cbind(u, return_dat)
      
      return(return_dat)
      
    }
  }
  
  return(OE)
  
}

#convert mutual information to SU
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

#amino acid pair SU
aap_SU <- function(amino_acid_matrix) {
  
  SU <- foreach(u = 1:(ncol(amino_acid_matrix)-1), .combine = 'c') %do% {
    
    count <-  cbind(amino_acid_matrix[,u], amino_acid_matrix[,u+1])
    
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
      
      normalize(mut, e1, e2)
      
    }
  }
  
  return(SU)
  
}

#codon-pair SU
cp_SU <- function(codon_matrix) {
  
  foreach(u = 1:(ncol(codon_matrix)-1), .combine = 'c') %do% {
    
    count <-  cbind(codon_matrix[,u], codon_matrix[,u+1])
    
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
      
      normalize(mut, e1, e2)
      
    }
  }
  
}
  