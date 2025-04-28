library(foreach)
library(doParallel)
library(cluster)
library(ggplot2)
library(seqinr)

#Character vector of gene names to analyze.
genes <- c()

for (g in genes) {
  
  #Set the working directory.
  pwd <- paste0()  
  setwd(pwd)
  
  load(file = paste0(g, "_c.mat"))

  #Count each position's codon usages.
  c_count <- foreach(c = 1:(ncol(cod_mat)/3)) %do% {
    
    plyr::count(cod_mat[,(1:3)+3*(c-1)])
    
  }
  
  #Convert the codon alignment to a codon-pair alignment.
  cod_mat <- foreach(c = 1:(ncol(cod_mat)/3), .combine = 'cbind') %:% foreach(r = 1:nrow(cod_mat), .combine = 'c') %do% {
    
    paste0(cod_mat[r,(1:3)+3*(c-1)], collapse = "")
    
  }
  
  #Sum the codon usages for the entire alignment.
  c_count_total <- foreach(n = 1:length(c_count), .combine = "rbind") %do% {
    
    c_count[[n]]
    
  }
  
  c_count_names <- foreach(n = 1:nrow(c_count_total), .combine = "c") %do% {
    
    paste0(c_count_total[n,1:3], collapse = "")
    
  }
  
  c_count_total <- c_count_total[,4]
  names(c_count_total) <- c_count_names
  
  c_count_total <-  tapply(c_count_total, names(c_count_total), sum) 
  real_codons <- toupper(seqinr::words())
  
  #Remove stop codons.
  real_codons <- real_codons[-c(which(real_codons=="TAG"), 
                                which(real_codons=="TGA"), 
                                which(real_codons=="TAA"))]
  
  c_count_real_names <- foreach(n = 1:length(c_count_total), .combine = "c") %do% {
    
    names(c_count_total)[n] %in% real_codons
    
  }
  
  c_count_total <- c_count_total[which(c_count_real_names==T)]
  c_count_total <- c_count_total[order(names(c_count_total))]

  #Translate codon counts to amino acid counts.
  aa_names <- foreach(r = real_codons, .combine = 'c') %do% {
    
    seqinr::translate(unlist(strsplit(r, split = character(0))))
    
  }
  
  aa_order <- order(aa_names)
  aa_names <- aa_names[aa_order]
  c_count_total <- c_count_total[aa_order]
  real_codons <- real_codons[aa_order]

  #Calculate alignment-specific rscu.
  rscu <- foreach(a = unique(aa_names), .combine = 'c') %do% {
    
    degen <- which(aa_names==a)
    nac_all <- sum(c_count_total[degen])
    
    RFac_all <- foreach(d = degen, .combine = 'c') %do% {
      
      c_count_total[d]/nac_all
      
    }
    
    RFac_all*length(degen)
    
  }
  
  #Count each position's codon-pair usages.
  cp_counts <- foreach(h = 1:(ncol(cod_mat)-1), .combine = 'rbind') %do% {
    
    plyr::count(cbind(cod_mat[,h], cod_mat[,h+1]))
    
  }

  #Sum the codon-pair counts for the alignment.
  cp_names <- foreach(y = 1:nrow(cp_counts), .combine = 'c') %do% {
    
    paste0(cp_counts[y, 1:2], collapse = "")
    
  }
  
  cp_counts <- cp_counts[,3]
  names(cp_counts) <- cp_names
  cp_counts <- tapply(cp_counts, names(cp_counts), sum)  
  
  cp_count_names <- foreach(c = 1:length(cp_counts), .combine = 'rbind') %do% {
    
    cp_string <- unlist(strsplit(names(cp_counts)[c], split = character()))
    c(paste0(cp_string[1:3], collapse = ""), paste0(cp_string[4:6], collapse = ""))
    
  }
  
  #Remove irrelevant codon-pairs.
  c1_good_names <- foreach(c = 1:nrow(cp_count_names), .combine = "c") %do% {
    
    cp_count_names[c,1] %in% real_codons
    
  }
  
  c2_good_names <- foreach(c = 1:nrow(cp_count_names), .combine = "c") %do% {
    
    cp_count_names[c,2] %in% real_codons
    
  }
  
  
  cp_good_names <- foreach(c = 1:nrow(cp_count_names), .combine = 'c') %do% {
    
    (c1_good_names[c] == T) && (c2_good_names[c] == T)
    
  }
  
  cp_counts <- cp_counts[which(cp_good_names==T)]
  
  #Codon-pair counts labelled by their amino acid.
  cp_info <- foreach(cp = 1:length(cp_counts), .combine = 'rbind') %do% {
    
    cp_id <- names(cp_counts)[cp]
    aa_id <- paste0(seqinr::translate(unlist(strsplit(cp_id, split = character(0)))), collapse = "")
    loop_count <- cp_counts[cp]
    
    c(cp_id, aa_id, loop_count)
    
  }
  
  cp_freqs <- as.numeric(cp_info[,3])/sum(as.numeric(cp_info[,3]))
  cp_info <- cbind(cp_info, cp_freqs)
  
  #Calculate alignment-specific rscpf
  rscpf <- foreach(a = unique(cp_info[,2]), .combine = 'c') %do% {

    degen <- which(cp_info[,2]==a)
    nac_all <- sum(as.numeric(cp_info[degen,3]))

    RFac_all <- foreach(d = degen, .combine = 'c') %do% {

      as.numeric(cp_info[d,3])/nac_all

    }

    no_name_rspcu <- RFac_all
    
    names(no_name_rspcu) <- cp_info[degen, 1]
    no_name_rspcu

  }

  #Make an empty dataframe then fill it with RSCPF values.
  cp_counts_df <- expand.grid(Codon1 = real_codons, Codon2 = real_codons)
  cp_vec <- rep(NA, times = nrow(cp_counts_df))


  foreach(e = 1:length(rscpf)) %do% {

    cp <- unlist(strsplit(names(rscpf)[e], split = character(0)))

    c1 <- which(cp_counts_df$Codon1 == paste0(cp[1:3], collapse = ""))
    c2 <- which(cp_counts_df$Codon2 == paste0(cp[4:6], collapse = ""))
    cp_vec[intersect(c1, c2)] <- rscpf[e]

  }

  cp_counts_df$Counts <- cp_vec
  ic_rscpf_df <- cp_counts_df
  save(ic_rscpf_df, file = "ic_rscpf.df")
  
  #Calculate codon, amino acid, and amino acid pair usage frequencies.
  c_freqs <- c_count_total/sum(c_count_total, na.rm = T)
  
  aa_counts <- c_count_total
  names(aa_counts) <- aa_names
  aa_counts <- tapply(aa_counts, names(aa_counts), sum)
  aa_freqs <- aa_counts/sum(aa_counts, na.rm = T)
  
  aa_pair_counts <- as.numeric(cp_info[,3])
  names(aa_pair_counts) <- cp_info[,2]
  aa_pair_counts <- tapply(aa_pair_counts, names(aa_pair_counts), sum)
  aa_pair_freqs <- aa_pair_counts/sum(aa_pair_counts)
  
  #Calculate alignment-specific cps.
  cps <- foreach(b = 1:nrow(cp_info), .combine = 'c') %do% {
    
    c_pair <- unlist(strsplit(cp_info[b, 1], split = character(0)))
    aa_pair<- unlist(strsplit(cp_info[b, 2], split = character(0)))
    
    fab <- aa_pair_freqs[which(names(aa_pair_freqs)==cp_info[b, 2])]
    fa <- aa_freqs[which(names(aa_freqs)==aa_pair[1])]
    fb <- aa_freqs[which(names(aa_freqs)==aa_pair[2])]
    fafb <- fa*fb
    
    c12 <- as.numeric(cp_info[b,4])
    c1 <- c_freqs[which(names(c_freqs) == paste0(c_pair[1:3], collapse = ""))]
    c2 <- c_freqs[which(names(c_freqs) == paste0(c_pair[4:6], collapse = ""))]
    c1c2 <- c1*c2
    
    loop_cps <- log10(c12/((c1c2/fafb)*fab))
    names(loop_cps) <- cp_info[b,1]
    loop_cps
    
  }

  #Make an empty dataframe then fill it with CPS values.
  cp_counts_df <- expand.grid(Codon1 = real_codons, Codon2 = real_codons)
  cp_vec <- rep(NA, times = nrow(cp_counts_df))
  
  foreach(e = 1:length(cps)) %do% {
    
    cp <- unlist(strsplit(names(cps)[e], split = character(0)))
    
    c1 <- which(cp_counts_df$Codon1 == paste0(cp[1:3], collapse = ""))
    c2 <- which(cp_counts_df$Codon2 == paste0(cp[4:6], collapse = ""))
    cp_vec[intersect(c1, c2)] <- cps[e]
    
  }
  
  cp_counts_df$Counts <- cp_vec
  
  ic_cps_df <- cp_counts_df
  save(ic_cps_df, file = "ic_cps.df")
  
}


