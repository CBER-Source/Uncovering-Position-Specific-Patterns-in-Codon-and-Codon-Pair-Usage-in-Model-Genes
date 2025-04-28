library(doParallel)
library(foreach)
library(read.gb)
library(stringr)

#Character vector of gene names to analyze.
genes <- c()

for (g in genes) {
  
  #Set the working directory.
  pwd <- paste0()
  setwd(pwd)
  
  #Read in metadata for sequence filtration.
  hit_table <- read.csv(file = paste0("hit.csv"))
  des_table <- read.csv(file = paste0("des.csv"))
  blast_results <- read.gb("sequence.gb")
  
  #Make a vector containing each sequence's transcript variant.
  transcript_variants <- foreach(l = 1:length(blast_results), .combine = 'c') %do% {
    
    string <- blast_results[[l]]$DEFINITION
    
    str_extract(string, "(?<=transcript variant )\\w+")
    
  }
  
  des_table <- cbind(des_table, transcript_variants)
  
  #Match species with accession numbers.
  species <- unique(des_table[,2])
  
  species_acc <- foreach(s = species) %do% {
    
    des_table[which(des_table[,2]==s),9]
    
  }
  
  #Match accession numbers with hits by species.
  species_hits <- foreach(s = 1:length(species_acc)) %:% foreach(a = 1:length(species_acc[[s]]), .combine = 'rbind') %do% {
    
    hit_table[which(hit_table[,2]==species_acc[[s]][a]),]
    
  }
  
  species_dref <- foreach(s = 1:length(species_hits)) %do% {
    
    species_hits[[s]][which(species_hits[[s]][,12]==max(species_hits[[s]][,12])),2]
    
  }

  #Which transcript has the maximum hit for each species.
  reference_accessions <- foreach(s = 1:length(species_dref), .combine = 'c') %do% {
    
    if (length(species_dref[[s]])>1) {
      
      transcript_variants <- foreach(n = species_dref[[s]], .combine = 'c') %do% {
        
        des_table[which(des_table[,9]==n),10]
        
      }
      
      transcript_ind <- which(transcript_variants==min(transcript_variants, na.rm = T))
      
      species_dref[[s]][transcript_ind]
      
    } else {
      species_dref[[s]]
    }
    
  }
    
  #Vector of taxids for each reference sequence.
  taxids <- foreach(r = reference_accessions, .combine = 'c') %do% {
    
    des_table[which(des_table[,9]==r),2]
    
  }
  
  #Vector of species for each reference sequence.
  species <- foreach(r = reference_accessions, .combine = 'c') %do% {
    
    des_table[which(des_table[,9]==r),1]
    
  }
  
  #Vector of gene descriptions for each reference sequence.
  gene_description <- foreach(r = reference_accessions, .combine = 'c') %do% {
    
    blast_results[[which(des_table[,9]==r)]]$DEFINITION
    
  }
  
  #Vector of gene symbols for each reference sequence.
  gene_symbol <- foreach(q = gene_description, .combine = 'c') %do% {
    
    gsub(".*\\(|\\).*", "", q)
    
  }
  
  #Output metadata for each species reference transcript as a .csv.
  gene_info <- data.frame(GeneAC = reference_accessions, 
                          GeneSymbol = gene_symbol, 
                          GeneDescription = gene_description, 
                          Species = species, 
                          Taxid = taxids
                          )

  write.csv(gene_info, file = paste0(g, " - Supplemental File 1.csv"))
  
}
