library(seqinr)
library(foreach)
library(doParallel)
library(ggplot2)
library(ggbreak)
library(stringr)
library(gridExtra)
library(cowplot)
library(ggh4x)
library(scales)
library(stringr)

#Set the parent working directory.
home_wd <- paste0()
setwd(home_wd)
source("SU functions.R")

#Character vector of gene names to analyze.
genes <- c()

#Set the results working directory.
genewd <- paste0()

for (g in genes) {

  #Set the working directory.
  pwd <- paste0()
  setwd(pwd)
  
  #Load RSCPF data.
  gene_df <- read.csv(file = paste(g, "Supplemental File 6.csv"), as.is = T, na.strings = NULL)
  
  #Extract human sequence for reference.
  load(paste0(g, "_c.mat"))
  
  if (g == "vwf_75") {
    human_ind <- 150
  } else if (g == "f8_75") {
    human_ind <- 111
  } else if (g == "f9_75") {
    human_ind <- 75
  } else if (g == "ADAMTS13_75") {
    human_ind <- 106
  } 
  
  #Convert the nucleotide sequence to codon, codon-pair, amino acid, and amino acid pair sequences.
  gene_wt <- cod_mat[human_ind,]

  gene_wt_codseq <- foreach(n = 1:(length(gene_wt)/3), .combine = 'c') %do% {
    
    paste0(gene_wt[1+3*(n-1)], gene_wt[2+3*(n-1)], gene_wt[3+3*(n-1)])
    
  }
  
  gene_wt_codseq <- gene_wt_codseq[-length(gene_wt_codseq)]
  
  gene_wt_cpseq <- foreach(n = 1:(length(gene_wt_codseq)-1), .combine = 'c') %do% {
    
    paste0(gene_wt_codseq[n], gene_wt_codseq[n+1])
    
  }
  
  save(gene_wt_cpseq, file = "gene_wt.cpseq")
  
  gene_wt_aaseq <- seqinr::translate(gene_wt)
  
  gene_wt_aapseq <- foreach(n = 1:(length(gene_wt_aaseq)-1), .combine = 'c') %do% {
    
    paste0(gene_wt_aaseq[n], gene_wt_aaseq[n+1])
    
  }
  
  #Generate codon-pair optimality data for the plotting script.
  supp6_update <- foreach(u = 1:length(gene_df$Position), .combine = 'rbind') %do% {
    
    position <- gene_df$Position[u]
    pos_dat <- gene_df[which(gene_df$Position==position),]
    
    #which amino acid pair is most common at this position? 
    pos_aa <- pos_dat$Position.specific.usage
    names(pos_aa) <- pos_dat$Amino.acid.pair
    
    pos_aa <- tapply(pos_aa, names(pos_aa), sum)
    pos_aa <- names(pos_aa)[which(pos_aa == max(pos_aa))]
    
    #which amino acid pair is u?
    ucp <- gene_df$Codon.pair[u]
    uaa <- gene_df$Amino.acid.pair[u]
    uaa_rscpf <- gene_df$Position.specific.RSCPF[u]
    
    #which amino acid pair is in the human sequence. 
    human_cp <- gene_wt_cpseq[position]
    human_aa <- gene_wt_aapseq[position]
    
    pos_aa_dat <- pos_dat[which(pos_dat$Amino.acid.pair==uaa),]
    
    #if u is synonymous to the human codon-pair.
    if (uaa == human_aa) {
      
      aln_max_dat <- gene_df[which(gene_df$Amino.acid.pair==uaa),]
      cp <- aln_max_dat$Codon.pair
      aln_rscpf <- aln_max_dat$Alignment.specific.RSCPF
      aln_max_dat <- plyr::count(cbind(cp, aln_rscpf))
      cp <- aln_max_dat[,1]
      aln_max_dat <- aln_max_dat[,2]
      names(aln_max_dat) <- cp
      
      aln_max <- aln_max_dat[which(aln_max_dat==max(aln_max_dat))]
      
      
      #only show results for most common amino acid pair at this position
      pos_max_ind <- which(pos_aa_dat$Position.specific.RSCPF==max(pos_aa_dat$Position.specific.RSCPF))
      pos_max <- pos_aa_dat$Position.specific.RSCPF[pos_max_ind]
      names(pos_max) <- pos_aa_dat$Codon.pair[pos_max_ind]
      
      #If u has maximum position-specific RSCPF.
      if (ucp %in% names(pos_max)) {
        
        #If u has maximum alignment-specific RSCPF.
        if (ucp %in% names(aln_max)) {
          
          #If u is the human codon-pair.
          if (ucp == human_cp) {
            pos_aln_compare <- "Human codon-pair was alignment- and position-specific optimal."
          } else {
            pos_aln_compare <- "Non-human codon-pair was alignment- and position-specific optimal."
          }
          #If u does not have maximum alignment-specific RSCPF.
        } else {
          #If u is the human codon-pair.
          if (ucp == human_cp) {
            pos_aln_compare <- "Human codon-pair was only position-specific optimal."
          } else {
            pos_aln_compare <- "Non-human codon-pair was position-specific optimal."
          }
          
        }
        #If u does note have maximum position-specific RSCPF.
      } else {
        #If u is the human codon-pair.
        if (ucp == human_cp) {
          if (ucp %in% names(aln_max)) {
            pos_aln_compare <- "Human was optimal alignment-specific codon-pair."
            
          } else {
            pos_aln_compare <- "Non-optimal codon-pair."
            
          }
          #If u is not the human codon-pair.
        } else {
          pos_aln_compare <- "Non-optimal codon-pair."
          
        }
        
      } 
      
      #if u is not synonymous to the human codon-pair.
    } else {
      pos_aln_compare <- "Non-synonymous with human codon-pair."
    } 
    
    
    c("g_name" = unique(pos_aa_dat$g_name), 
      "Position" = gene_df$Position[u],
      "Codon.pair" = gene_df$Codon.pair[u],
      "Optimization.Candidate" = pos_aln_compare)
    
  }
  
  save(supp6_update, file = paste0(g, " supplemental file 6 - aln-pos.update"))
  
  #Generate data that is summarized in table 5.
  table5 <- foreach(u = unique(gene_df$Position), .combine = 'rbind') %do% {
    
    position <- as.numeric(u)
    
    #which amino acid pair is in the human sequence?
    human_cp <- gene_wt_cpseq[position]
    human_aa <- gene_wt_aapseq[position]
    
    pos_dat <- gene_df[which(gene_df$Position==position),]
    pos_dat <- pos_dat[which(pos_dat$Amino.acid.pair==human_aa),]
    
    pos_cp_opt <- pos_dat$Codon.pair[which(pos_dat$Position.specific.RSCPF==max(pos_dat$Position.specific.RSCPF))]

    
    aln_max_dat <- gene_df[which(gene_df$Amino.acid.pair==human_aa),]
    cp <- aln_max_dat$Codon.pair
    aln_rscpf <- aln_max_dat$Alignment.specific.RSCPF
    aln_max_dat <- plyr::count(cbind(cp, aln_rscpf))
    cp <- aln_max_dat[,1]
    aln_max_dat <- aln_max_dat[,2]
    names(aln_max_dat) <- cp
    
    aln_cp_opt <- aln_max_dat[which(aln_max_dat==max(aln_max_dat))]
    
    
    c("g_name" = unique(pos_dat$g_name), 
      "Position" = u,
      "Human amino acid pair" = human_aa,
      "Human codon-pair" = human_cp,
      "Alignment-specific optimal codon-pair" = paste(names(aln_cp_opt), collapse = ", "),
      "Alignment-specific change?" = !(human_cp %in% names(aln_cp_opt)),
      "Position-specific codon-pair" = paste(pos_cp_opt, collapse = ", "),
      "Position-specific change?" = !(human_cp %in% pos_cp_opt),
      "Alignment-specific optimal codon-pair = position-specific optimal codon-pair" = (names(aln_cp_opt) %in% pos_cp_opt))
    
  }
  
  #Count the optimization results for each position.
  aln_vs_pos <- plyr::count(table5[,c(6,8)])
  aln_equal_pos <- plyr::count(table5[,9])
  aln_equal_pos <- cbind("Aln == Pos opt cp?", aln_equal_pos)
  colnames(aln_equal_pos) <- colnames(aln_vs_pos)
  
  aln_vs_pos <- rbind(aln_vs_pos, aln_equal_pos)
  write.csv(aln_vs_pos, file = paste0(g, "_table5-cp-position.csv"))
  
  pos_opt_type <- c("WT codon-pair is as common as other codon-pair(s) at this position.",
                    "WT codon-pair is less common than other, equally common codon-pairs at this position.",
                    "WT codon-pair is the most common codon-pair at this position.",
                    "WT codon-pair is less common than another codon-pair at this position.",
                    "WT codon-pair was not observed at this position.")
  
}
