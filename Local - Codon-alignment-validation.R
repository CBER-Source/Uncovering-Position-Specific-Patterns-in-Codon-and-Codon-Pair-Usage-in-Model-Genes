library(Biostrings)
library(foreach)
library(doParallel)

#Character vector of gene names to analyze.
genes <- c()

#Set the results working directory.
genewd <- paste0()
setwd(genewd)

for (g in genes) {
 
  #Set the working directory.
  pwd <- paste0()
  setwd(pwd)
  
  #Read coding sequences and extract their annotations.
  refseq_set <- seqinr::read.fasta(file = paste0(g, "_entrez.fasta"))
  
  refseq_annotations <- foreach(r = 1:length(refseq_set), .combine = 'c') %do% {
    
    annot <- unlist(strsplit(attributes(refseq_set[[r]])$Annot, split = " "))
    paste(annot[2:length(annot)], collapse = " ")
    
  }
  
  #Read aligned sequences.
  macseq_set <- seqinr::read.fasta(file = paste0(g, "_entrez_NT_NT_NT.fasta"))
  
  nm <- names(macseq_set)
  rm <- names(refseq_set)
  
  #Ensure the aligned and unaligned sequences are in the same order.
  ncheck <- foreach(n = nm, .combine = 'c') %do% {
    
    which(rm == n)
    
  }
  
  refseq_set <- seqinr::read.fasta(file = paste0(g, "_entrez.fasta"), seqonly = T)
  macseq_set <- seqinr::read.fasta(file = paste0(g, "_entrez_NT_NT_NT.fasta"), seqonly = T)
  
  #Read in metadata
  acc <- read.csv("reference_accessions.csv")
  des <- read.csv(paste0(g, "_des.csv"))
  
  rm <- rm[ncheck]
  refseq_set <- refseq_set[ncheck]
  

  stat_sum <- foreach(s = 1:length(macseq_set), .combine = 'rbind') %do% {
    
    #Index the aligned nucleotide sequences by codon and remove gap codons.
    mac_align <- gsub("(.{3})", "\\1 ", macseq_set[[s]])
    mac_align <- unlist(strsplit(mac_align, split = " "))
    gap_inds <- which(mac_align == "---")

    mac_align <- mac_align[-gap_inds]
    mac_align <- paste0(mac_align, collapse = "")

    #nt pairwise alignment
    ref_align <- refseq_set[[s]]
    
    palign <- pairwiseAlignment(ref_align, mac_align)
    
    fname <-  acc[s,]
    
    #Calculate percent identity for the alignment.
    percent_identity <- pid(palign)
    
    #If the sequences are not identical, identify which codons have changed.
    if (percent_identity!=100) {
      writePairwiseAlignments(palign, file = paste0(fname, "_orf_pairwise_alignment.txt"))
      
      ref_codon <- pattern(palign)
      ref_codon <- gsub("(.{3})", "\\1 ", ref_codon)
      ref_codon <- unlist(strsplit(ref_codon, split = " "))
      
      mac_codon <- subject(palign)
      mac_codon <- gsub("(.{3})", "\\1 ", mac_codon)
      mac_codon <- unlist(strsplit(mac_codon, split = " "))
      
      unmatched_codons <- ref_codon==mac_codon
      ref_codon <- ref_codon[which(unmatched_codons==F)]
      mac_codon <- mac_codon[which(unmatched_codons==F)]
      
      ref_aa <- seqinr::translate(unlist(strsplit(paste0(ref_codon, collapse = ""), split = character(0))))
      mac_aa <- seqinr::translate(unlist(strsplit(paste0(mac_codon, collapse = ""), split = character(0))))
      
      diff_table <- cbind(codon_ind = which(unmatched_codons==F), 
                          ref_codon, 
                          mac_codon,
                          ref_aa,
                          mac_aa)
      
      write.csv(diff_table, file = paste0(fname, "_codon_differences.csv"))
      
    }
    
    #Output the accession number and its percent identity.
    c(fname, percent_identity)
    
  }
  
  stat_numeric <- as.numeric(stat_sum[,2])
  write.csv(stat_sum, file = paste0(g, "_pid_summary.csv"))
  
}
