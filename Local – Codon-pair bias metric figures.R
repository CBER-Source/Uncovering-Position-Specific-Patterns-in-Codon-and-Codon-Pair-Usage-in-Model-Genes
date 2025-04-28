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

#Character vector of gene names to analyze.
genes <- c()

#Set the results working directory.
genewd <- paste0()

#Dataframe column names for reference during loops.
df_colnames <- c("g_name","Codon-pair",             
                    "Amino acid pair",                 "Position"  ,           
                    "Alignment-specific usage",        "Position-specific usage" ,           
                    "Alignment-specific frequency",   "Position-specific frequency" , 
                    "Alignment-specific CPS",                  "Position-specific CPS"   ,              
                    "Alignment-specific RSCPF" ,               "Position-specific RSCPF"  ,             
                    "Alignment-specific Observed to Expected Ratio" ,              "Position-specific Observed to Expected Ratio" ,             
                    "Domain",                         "Disease associated synonymous SNV",
                    "PTM", "Synonymous?", "Amino acid pair conservation", 
                    "Codon-pair Symmetric Uncertainty", "Amino acid pair Symmetric Uncertainty")

#Customized domain colors for the four human reference sequences.
domain_color_df <- list("vwf_75" = c(red = 4,	blue = 4,	green = 3,	purple = 3,	orange = 1), 
                        "f8_75" = c(red = 6,	blue = 2,	green = 0,	purple = 0,	orange = 0), 
                        "f9_75" = c(red = 1,	blue = 2,	green = 1,	purple = 0,	orange = 0), 
                        "ADAMTS13_75" = c(red = 1,	blue = 1,	green = 8,	purple = 2,	orange = 0))

domain_hues_df <- foreach(gd = 1:length(domain_color_df)) %do% {
  
  c("white", foreach(u = 1:length(domain_color_df[[gd]]), .combine = 'c') %do% {
    
    domain_group_n <- domain_color_df[[gd]][u]
    
    if (domain_group_n == 1) {
      names(domain_color_df[[gd]][u])
    } else if (domain_group_n > 1) {
      
      domain_group_color <- names(domain_color_df[[gd]][u])
      
      percentages <- 1:domain_group_n/(domain_group_n+0.5*domain_group_n)
      
      colorspace::darken(domain_group_color, percentages)
      
    } 
    
  })}

names(domain_hues_df) <- c("vwf_75", "f8_75", "f9_75", "ADAMTS13_75")

#Reorder the vwf domain colors to match the sequence.
new_order <- c(1, 2, 6, 3, 7, 8, 4, 9,10,11, 12, 5, 13, 14, 15, 16)
domain_hues_df$vwf_75 <- domain_hues_df$vwf_75[new_order]

#Generate dataframe for figures.
gene_df <- foreach(g = genes, .combine = 'rbind') %do% {
  
  #Set the working directory.
  pwd <- paste0()
  setwd(pwd)
  
  if (g=="vwf_75") {
    g_name <- "VWF"
  } else if (g=="f8_75") {
    g_name <- "F8"
  } else if (g=="f9_75") {
    g_name <- "F9"
  } else if (g=="ADAMTS13_75") {
    g_name <- "ADAMTS13"
  }
  
  #Load in alignment-specific codon-pair bias metric results.
  load(file = "aln_cps.df")
  aln_cps_df$CodonPair <- paste0(aln_cps_df$Codon1, aln_cps_df$Codon2)
  aln_cps_df$SequentialCPS <- aln_cps_df$Counts
  aln_cps_df$Codon1 <- aln_cps_df$Codon2<- aln_cps_df$Counts <- NULL
  
  load(file = "aln_obsexp.df")
  aln_obsexp_df$CodonPair <- paste0(aln_obsexp_df$Codon_1, aln_obsexp_df$Codon_2)
  aln_obsexp_df$SequentialObsExp <- aln_obsexp_df$ObservedExpected
  aln_obsexp_df$Codon_1 <- aln_obsexp_df$Codon_2 <- aln_obsexp_df$ObservedExpected <- NULL
  
  load(file = "aln_RSCPF.df")
  aln_RSCPF_df$CodonPair <- paste0(aln_RSCPF_df$Codon1, aln_RSCPF_df$Codon2)
  aln_RSCPF_df$SequentialRSCPF <- aln_RSCPF_df$Counts
  aln_RSCPF_df$Codon1 <- aln_RSCPF_df$Codon2 <- aln_RSCPF_df$Counts <- NULL
  
  #Load and format position-specific codon-pair bias metric results.
  load(paste0(g, "_cpbm.calc"))
  cpbm_calc <- as.data.frame(cpbm_calc)
  cpbm_calc$CodonPairsObserved <- as.numeric(cpbm_calc$CodonPairsObserved)
  cpbm_calc$CPS <- as.numeric(cpbm_calc$CPS)
  cpbm_calc$RSCPF <- as.numeric(cpbm_calc$RSCPF)
  cpbm_calc$ObsExp <- as.numeric(cpbm_calc$ObsExp)
  names(cpbm_calc)[4:6] <- paste0("Positional", names(cpbm_calc)[4:6])
  
  #Add alignment-specific results to the dataframe.
  calc_plot_df <- merge(x = cpbm_calc, y = aln_cps_df)
  calc_plot_df <- merge(x = calc_plot_df, y = aln_obsexp_df)
  calc_plot_df <- merge(x = calc_plot_df, y = aln_RSCPF_df)
  
  ############################################
  #Load the human reference alignment and index by codons.
  load(paste0(g, "_c.mat"))
  
  cod_mat <- foreach(n = 1:(ncol(cod_mat)/3), .combine = 'cbind') %do% {
    
    paste0(cod_mat[,1+3*(n-1)], cod_mat[,2+3*(n-1)], cod_mat[,3+3*(n-1)])
    
  }
  
  cod_mat <- cod_mat[,1:(ncol(cod_mat)-1)]

  #Read feature info
  gff <- read.delim(paste0(g, " - domains - uniprot - 3-27-24.gff"), header=F, comment.char="#")
  
  #Add domains to the dataframe.
  domain_rows <- which(gff[,3]=="Domain")
  gff <- gff[domain_rows,]
  
  domain_vec <- foreach(r = 1:nrow(gff), .combine = 'c') %do% {
    
    c(gff[r,4]:gff[r,5])
    
  }
  
  #Identify non-domain positions.
  no_domain_inds <- setdiff(1:(ncol(cod_mat)-1), domain_vec)
  no_domain_inds <- unname(split(no_domain_inds, cumsum(c(TRUE, diff(no_domain_inds)!=1))))
  
  no_domain_df <- foreach(l = 1:length(no_domain_inds), .combine = 'rbind') %do% {
    
    no_name <- paste0("No domain - ", l)
    cbind(no_domain_inds[[l]], no_name)
    
  }
  
  #Label each position with its domains.
  domain_df <- foreach(r = 1:nrow(gff), .combine = 'rbind') %do% {
    
    note <- gff[r,9]
    
    if (grepl(";", note)) {
      
      note <- sub("\\;.*", "", note)
      
    } 
    
    note <- unlist(strsplit(note, split = character(0)))
    note <- note[6:length(note)]
    
    domain <- paste0(note, collapse = "")
    cbind(gff[r,4]:gff[r,5], domain)
    
  }
  
  domain_df <- data.frame(Position = as.numeric(domain_df[,1]), Domain = domain_df[,2])
  no_domain_df <- data.frame(Position = as.numeric(no_domain_df[,1]), Domain = no_domain_df[,2])
  
  #Remove stop from domain data.
  if (max(domain_df$Position)>(ncol(cod_mat)-1)) {
    stop_ind <- which(domain_df$Position==max(domain_df$Position))
    domain_df <- domain_df[-stop_ind,]
  }
  
  domain_df <- rbind(domain_df, no_domain_df)
  position_plot_df <- data.frame(Position = 1:(ncol(cod_mat)-1))
  
  position_plot_df <- merge(position_plot_df, domain_df)
  
  #Reformat the domains to factors to maintain their order during plotting.
  dnames <- unique(position_plot_df$Domain)
  domain_order <- unique(position_plot_df$Domain)
  position_plot_df$DomainFactor <- factor(position_plot_df$Domain, domain_order)
  dnumber <- unique(position_plot_df$Domain)
  dpair <- rbind(dnames, dnumber)
  
  DomainPanelName <- foreach(d = position_plot_df$Domain, .combine = 'c') %do% {
    
    CheckForNoDomain <- grepl("No domain*", d)  
    
    if (CheckForNoDomain==T) {
      paste("No domain")
    } else {
      paste(d)
    }
  }
  
  domain_order <- domain_order[-grep("No", domain_order)]
  domain_order <- c("No domain", domain_order)
  position_plot_df$DomainPanelName <- factor(DomainPanelName, domain_order, ordered = T)
  
  names(position_plot_df)[1] <- "CodonPairPosition"

  #Add domains and non-domains to the dataframe.
  calc_plot_df <- merge(calc_plot_df, position_plot_df, all = T)
  
  #######################################################
  
  #Calculate codon-pair usage quartiles from alignment-specific usages.
  cp_counts <- calc_plot_df$CodonPairsObserved
  names(cp_counts) <- calc_plot_df$CodonPair
  cp_counts <- tapply(cp_counts, names(cp_counts), sum)
  cp_counts <- sort(cp_counts, decreasing = T)
  
  #4 is 76-100% quartile
  #1 is 0-25% quartile
  bins <- infotheo::discretize(cp_counts, disc = "equalfreq", nbins = 4)
  names(bins) <- "CodonPairUsageBin"
  
  ####
  #Extract the most observed codon-pair from each quartile as a case study.
  plot_cp <- c(which((bins$CodonPairUsageBin==4)==T)[1],
               which((bins$CodonPairUsageBin==3)==T)[1],
               which((bins$CodonPairUsageBin==2)==T)[1])
  
  plot_cp <- data.frame(CodonPair = names(cp_counts)[plot_cp], PlotCP = c(4,3,2))
  
  #Add quartile results to the dataframe.
  calc_plot_df <- merge(x = calc_plot_df, y = plot_cp, by = "CodonPair", all = T)
  bin_info <- data.frame(CodonPair = names(cp_counts), TotalCodonPairsObserved = cp_counts, CodonPairUsageBin = bins)
  calc_plot_df <- merge(x = calc_plot_df, y = bin_info, by = "CodonPair")
  calc_plot_df$CodonPairUsageBin <- as.character(calc_plot_df$CodonPairUsageBin)
  
  #######################################################
  
  #Add amino acid pairs to the dataframe.
  aa_pairs <- seqinr::translate(unlist(strsplit(paste0(calc_plot_df$CodonPair, collapse = ""), split = character(0))))
  aa_pairs <- paste0(aa_pairs, collapse = "")
  aa_pairs <- gsub("(.{2})", "\\1 ", aa_pairs)
  aa_pairs <- unlist(strsplit(aa_pairs, split = " "))
  calc_plot_df$AminoAcidPairs <- aa_pairs
  
  #######################################################

  #Add amino acid pair conservation to the dataframe.
  setwd(pwd)
  consurf_aa <- as.data.frame(readxl::read_excel('consurf_grades-aa.xlsx', skip = 27, col_names = T))
  
  if (ncol(cod_mat)!=nrow(consurf_aa)) {
    stop("Codon numbers don't match. ")
  }
 
  #Estimate amino acid pair conservation as the average of adjacent amino acid conservations.
  aap_conservation <- foreach(a = 1:(nrow(consurf_aa)-1), .combine = 'c') %do% {
    
    mean(c(consurf_aa[a,4], consurf_aa[a+1,4]))
    
  }
  
  #Add amino acid pair conservation to the dataframe.
  consurf_aa <- data.frame(Position = 1:length(aap_conservation), Conservation = aap_conservation)
  names(consurf_aa) <- c("CodonPairPosition", "Conservation")
  calc_plot_df <- merge(x = calc_plot_df, y = consurf_aa, by = "CodonPairPosition")

  #######################################################

  #read feature info
  gff <- read.delim(paste0(g, " - domains - uniprot - 3-27-24.gff"), header=F, comment.char="#")
  
  #Extract relevant metadata for PTMs.
  ptm_types <- c("Disulfide bond", "Glycosylation", "Modified residue")
  ptm_df <- foreach(p = ptm_types, .combine = 'rbind') %do% {
    
    gff[which(gff[,3]==p),]
    
  }
  
  #Identify positions with ptms.
  ptm_vec <- foreach(r = 1:nrow(ptm_df), .combine = 'c') %do% {
    
    if (ptm_df[r,3]=="Disulfide bond") {
      c(ptm_df[r,4],ptm_df[r,5])
    } else {
      c(ptm_df[r,4]:ptm_df[r,5])
    }
    
  }
  
  #Identify non-ptms positions.
  no_ptm_inds <- setdiff(1:(ncol(cod_mat)-1), unique(ptm_vec))
  no_ptm_inds <- unname(split(no_ptm_inds, cumsum(c(TRUE, diff(no_ptm_inds)!=1))))
  
  no_ptm_df <- foreach(l = 1:length(no_ptm_inds), .combine = 'rbind') %do% {
    
    no_name <- paste0("No ptm - ", l)
    cbind(no_ptm_inds[[l]], no_name)
    
  }
  
  #Format non-ptms and ptms for merging with the dataframe.
  no_ptm_df <- data.frame("CodonPairPosition"=as.numeric(no_ptm_df[,1]), 
                          "PTM"=NA)
  
  ptm_df <- foreach(r = 1:nrow(ptm_df), .combine = 'rbind') %do% {
    
    ptm <- ptm_df[r,3]

    if (ptm=="Disulfide bond") {
      inds <- c(ptm_df[r,4],ptm_df[r,5])
    } else if (ptm == "Glycosylation") {
      inds <- c(ptm_df[r,4]:ptm_df[r,5])
    } else if (ptm == "Modified residue") {
      inds <- c(ptm_df[r,4]:ptm_df[r,5])
    } else {
      stop("Unrecognized PTM")
    }
    
    cbind(inds, ptm)
    
  }
  
  #Add ptms and non-ptms to the dataframe.
  ptm_df <- data.frame("CodonPairPosition"=as.numeric(ptm_df[,1]), "PTM"=ptm_df[,2])
  ptm_df <- merge(ptm_df, no_ptm_df, all = T)
  calc_plot_df <- merge(x = calc_plot_df, y = ptm_df, by = "CodonPairPosition")
  
  #######################################################

  #Calculate sequence-specific codon-pair usage frequencies.
  calc_plot_df$SequentialCodonPairFrequency <- calc_plot_df$TotalCodonPairsObserved/sum(calc_plot_df$TotalCodonPairsObserved)
  
  #Calculate position-specific codon-pair usage frequencies.
  unique_cp_positions <- unique(calc_plot_df$CodonPairPosition)
  calc_plot_df$PositionalCodonPairFrequency <- NA
  
  foreach(u = unique_cp_positions) %do% {
    
    dat_inds <- which(calc_plot_df$CodonPairPosition==u)
    dat <- calc_plot_df[dat_inds, ]
    dat <- dat$CodonPairsObserved/sum(dat$CodonPairsObserved)
    calc_plot_df$PositionalCodonPairFrequency[dat_inds] <- dat
    
  }
  
  #Add alignment- and position-specific codon-pair usage frequencies to the dataframe.
  calc_plot_df <- cbind(g_name, calc_plot_df)
  
  #######################################################
  
  #Extract the human sequence from the human reference alignment.
  if (g=="vwf_75") {
    human_row <- 150
  } else if (g=="f8_75") {
    human_row <- 111
  } else if (g=="f9_75") {
    human_row <- 75
  } else if (g=="ADAMTS13_75") {
    human_row <- 106
  } else {
    stop("Error: This gene is not one of this study's human reference genes.")
  }
  
  human_seq <- cod_mat[human_row,]
  
  #Convert the human codon sequence into amino acid, amino acid pair, and codon-pair sequences.
  human_aa_seq <- seqinr::translate(unlist(strsplit(human_seq, split = character(0))))

  human_cp_seq <- foreach(p = 1:(length(human_seq)-1), .combine = 'c') %do% {
    paste0(human_seq[p], human_seq[p+1], collapse = "")
  }
  
  human_aap_seq <- foreach(p = 1:(length(human_cp_seq)-1), .combine = 'c') %do% {
    paste0(seqinr::translate(unlist(strsplit(human_cp_seq[p], split = character(0)))), collapse = "")
  }
  
  #Calculate percent identity and percent coverage for each homolog.
  pid_mat <- cod_mat[-human_row,]
  
  pid_per_seq <- foreach(r = 1:nrow(pid_mat), .combine = 'cbind') %do% {
  
    #codon identities
    row_seq <- pid_mat[r,]
    pid_vec <- row_seq==human_seq
    id <- plyr::count(pid_vec)
    id_names <- id[,1]
    id <- id[,2]
    names(id) <- id_names
    id <- id/sum(id)*100
    codon_id <- id[which(id_names=="TRUE")]
    
    #codon coverages
    codon_coverage <- (length(row_seq)-length(which(row_seq=="---")))/length(row_seq)*100
    
    #amino acid percent identity
    row_aa_seq <- seqinr::translate(unlist(strsplit(row_seq, split = character(0))))
    aa_pid_vec <- row_aa_seq==human_aa_seq
    aa_id <- plyr::count(aa_pid_vec)
    aa_id_names <- aa_id[,1]
    aa_id <- aa_id[,2]
    names(aa_id) <- aa_id_names
    aa_id <- aa_id/sum(aa_id)*100
    amino_acid_id <- aa_id[which(aa_id_names=="TRUE")]
    
    #codon_pair identities
    row_cp_seq <- foreach(s = 1:(length(human_seq)-1), .combine = 'c') %do% {
      paste0(row_seq[s], row_seq[s+1], collapse = "")
    }
    
    cp_pid_vec <- row_cp_seq==human_cp_seq
    cp_id <- plyr::count(cp_pid_vec)
    cp_id_names <- cp_id[,1]
    cp_id <- cp_id[,2]
    names(cp_id) <- cp_id_names
    cp_id <- cp_id/sum(cp_id)*100
    codon_pair_id <- cp_id[which(cp_id_names=="TRUE")]
    
    #codon_pair coverages
     cp_coverage <- (length(row_cp_seq)-length(grep("---", row_cp_seq)))/length(row_cp_seq)*100
    
    #amino acid pair percent identity
     row_aap_seq <- foreach(p = 1:(length(row_cp_seq)-1), .combine = 'c') %do% {
       paste0(seqinr::translate(unlist(strsplit(row_cp_seq[p], split = character(0)))), collapse = "")
     }
     
     aap_pid_vec <- row_aap_seq==human_aap_seq
     aap_id <- plyr::count(aap_pid_vec)
     aap_id_names <- aap_id[,1]
     aap_id <- aap_id[,2]
     names(aap_id) <- aap_id_names
     aap_id <- aap_id/sum(aap_id)*100
     amino_acid_pair_id <- aap_id[which(aap_id_names=="TRUE")]
     
    c(codon_id, codon_pair_id, codon_coverage, cp_coverage, amino_acid_id, amino_acid_pair_id)
    
  }
  
  rownames(pid_per_seq) <- c("Codon Percent Identity", "Codon-pair Percent Identity", "Codon Coverage", "Codon-pair Coverage", "Amino Acid Percent Identity", "Amino Acid Pair Percent Identity")
  
  #######################################################

  #How many codons in the codon-pair are synonymous to human for each entry at each position?
  syn_count_by_position <- foreach(r = 1:nrow(calc_plot_df), .combine = 'c') %do% {
    
    cp_pos <- as.numeric(calc_plot_df$CodonPairPosition[r])
    table_aap <- unlist(strsplit(calc_plot_df$AminoAcidPairs[r], split = character(0)))
    human_aap <- unlist(strsplit(human_aap_seq[cp_pos], split = character(0)))
    syn_count <- table_aap == human_aap
    
    syn_count <- plyr::count(syn_count)
    syn_count_names <- syn_count[,1]
    syn_count <- syn_count[,2]
    names(syn_count) <- syn_count_names
    
    out <- syn_count[which(syn_count_names=="TRUE")]
    if (length(out)==0) {
      0
    } else {
      out
    }
    
  }
  
  #Add synonymous codon counts to the dataframe.
  calc_plot_df$Synonymous <- syn_count_by_position
  calc_plot_df <- calc_plot_df[order(calc_plot_df$CodonPair),]
  
  #######################################################
  #Add original SU to data.frame
  load(file = "ref.su")
  load(file = "aap_ref.su")
  
  su_df <- data.frame(CodonPairPosition = 1:length(ref_SU), OriginalSU = ref_SU, AapOriginalSU = aap_ref_SU)
  calc_plot_df <- merge(x = calc_plot_df, y = su_df, by = "CodonPairPosition")
  
  #Output the gene's dataframe to be combined with the other genes.
  return(calc_plot_df)
} #end gene_df loop

#Add alignment- and position-specific amino acid pair observed to expected ratios to the dataframe.
new_data <- foreach(g = genes, .combine = 'rbind') %do% {
  
  pwd <- paste0()
  setwd(pwd)
  
  #position-specific amino acid pair observed to expected ratio
  original_aaoe <- read.csv(file = paste0(g, "_position-specific_original_aaoe.csv"))
  original_aaoe <- original_aaoe[,2:ncol(original_aaoe)]
  original_aaoe$aap <- paste0(original_aaoe[,3], original_aaoe[,4])
  original_aaoe <- original_aaoe[,c(1,2,6,5)]
  colnames(original_aaoe) <- c(df_colnames[c(4, 1, 3)],
                               "Position-specific amino acid pair Observed to Expected Ratio")
  
  #alignment-specific amino acid pair observed to expected ratio
  aln_original_aaoe <- read.csv(file = paste0(g, "_alignment-specific_original_aaoe.csv"))
  aln_original_aaoe <- aln_original_aaoe[,2:ncol(aln_original_aaoe)]
  aln_original_aaoe$aap <- paste0(aln_original_aaoe[,2], aln_original_aaoe[,3])
  aln_original_aaoe <- aln_original_aaoe[,c(1,5,4)]
  
  colnames(aln_original_aaoe) <- c(df_colnames[c(1, 3)],
                                   "Alignment-specific amino acid pair Observed to Expected Ratio")
  
  #Make a temporary dataframe to continue adding results to.
  dat <- merge(x = aln_original_aaoe, y = original_aaoe, by = df_colnames[c(1, 3)])
  
  
  #Load optimization results for each codon-pair.
  load(file = paste0(g, " supplemental file 6 - aln-pos.update"))
  
  #Add metadata for merging with other dataframes.
  cp_characters <- strsplit(supp6_update[,3], split = character(0))
  
  aap <- foreach(cp = cp_characters, .combine = 'c') %do% {
    
    paste0(seqinr::translate(unlist(cp)), collapse = "")
    
  }
 
  supp6_update <- cbind(supp6_update, 
                        aap)
  
  colnames(supp6_update) <- c(colnames(supp6_update)[1:(ncol(supp6_update)-1)],
                              df_colnames[3])
  supp6_update <- as.data.frame(supp6_update)
  
  #Add the codon-pair optimization results to the temporary dataframe.
  dat <- merge(x = dat, y = supp6_update, by = df_colnames[c(3, 4)])
  dat <- dat[,c(1,2,4:8)]
  colnames(dat) <- c("Amino.acid.pair",                                              
                     "Position",                                                     
                     "Alignment-specific.amino.acid.pair.Observed.to.Expected.Ratio",
                     "Position-specific.amino.acid.pair.Observed.to.Expected.Ratio", 
                     "g_name",                                                     
                     "Codon.pair",                                                   
                     "Compare")
    
    
  return(dat)
   
}

names(new_data) <- c("Aminoacidpair",                                              
  "CodonPairPosition",                                                     
  "Alignment-specific.amino.acid.pair.Observed.to.Expected.Ratio",
  "Position-specific.amino.acid.pair.Observed.to.Expected.Ratio", 
  "g_name",                                                     
  "Codonpair",                                                   
  "Compare")

gene_df <- merge(gene_df, new_data, by = c("g_name",
                                              "CodonPair",
                                              "CodonPairPosition"))

#######################################################


#Plot distributions of each human reference alignment’s relative synonymous codon-pair frequencies
#Gene results are separated so they can be plotted individually.
fig6_dat <- foreach(u = unique(gene_df$g_name), .combine = 'rbind') %do% {
  
  g_dat <- gene_df[which(gene_df$g_name==u),]
  
  #Alignment-specific results are repeated.
  #This loop extracts one result per codon-pair.
  alignment_specific_rscpf <- foreach(a = unique(g_dat$CodonPair), .combine = 'rbind') %do% {
    
    first_entry_ind <- which(g_dat$CodonPair==a)[1]
    
    c("Gene" = g_dat$g_name[first_entry_ind], 
      "CodonPairPosition" = NA, 
      "CodonPair" = a, 
      "RSCPF" = g_dat$SequentialRSCPF[first_entry_ind], 
      "CPOE" = g_dat$SequentialObsExp[first_entry_ind], 
      "CPS" = g_dat$SequentialCPS[first_entry_ind], 
      "Method" = "Alignment-specific",
      "Compare" = NA)
    
  }
  
  position_specific_rscpf <- g_dat[,c(1,2,3,6,7,5, ncol(g_dat))]
  position_specific_rscpf$Method = "Position-specific"

  names(position_specific_rscpf) <- c("Gene", "CodonPairPosition", "CodonPair", "RSCPF", "CPOE", "CPS", "Compare", "Method")
  
  position_specific_rscpf <- position_specific_rscpf[,c(1:6,8,7)]
  
  #Output in long format.
  rbind(alignment_specific_rscpf, position_specific_rscpf)
  
}

#Reformat data to numeric for plotting.
fig6_dat$RSCPF <- as.numeric(fig6_dat$RSCPF)
fig6_dat$CPOE <- as.numeric(fig6_dat$CPOE)
fig6_dat$CPS <- as.numeric(fig6_dat$CPS)

#Plot alignment- and position-specific RSCPF distributions as smoothed histograms.
ggplot(fig6_dat, aes(x=RSCPF, fill=Method, colour=Method)) +
  geom_histogram(aes(y=after_stat(density)), breaks=seq(0,1,.05), alpha=0.6, 
                 position="dodge", lwd=0.2) +
  ylab("Codon-pair density") +
  
  coord_cartesian(expand = F) +
  facet_wrap(~Gene) +
  theme(panel.spacing = unit(25, "pt"),
        legend.position = "bottom") +
  ggtitle(paste("Histogram of Normalized RSCPF"))

setwd(genewd)
ggsave(filename = "Figure 6.tiff", width = 7.5, height = 5.25, units = "in")

#Amino acid degeneracies for plotting.
aa_degen <- c(A=4,
              R=6,
              N=2,
              D=2,
              C=2,
              Q=2,
              E=2,
              G=4,
              H=2,
              I=3,
              L=6,
              K=2,
              M=1,
              "F"=2,
              P=4,
              S=6,
              "T"=4,
              W=1,
              Y=2,
              V=4
)

#######################################################

#One from each quartile for each gene (Figure 7 and Supplemental File 12)
for (u in unique(gene_df$g_name)) {
 
  #Convert plot names to file names.
  if (u=="VWF") {
    u_name <- "vwf_75"
  } else if (u=="F8") {
    u_name <- "f8_75"
  } else if (u=="F9") {
    u_name <- "f9_75"
  } else if (u=="ADAMTS13") {
    u_name <- "ADAMTS13_75"
  }
  
  #Set the results and working directory.
  genewd <- paste0()
  pwd <- paste0()
  
  #Extract data for plot.
  g_dat <- gene_df[which(gene_df$g_name==u),]
  plot_aap <- g_dat[!is.na(g_dat$PlotCP),]
  plot_aap <- unique(plot_aap$AminoAcidPairs)
  
  #Plot the case study for each quartile.
  for (a in plot_aap) {
    
    sep_aa <- unlist(strsplit(a, split = character(0)))
    degen1 <- aa_degen[which(names(aa_degen)==sep_aa[1])]
    degen2 <- aa_degen[which(names(aa_degen)==sep_aa[2])]
    aap_degen <- degen1 * degen2
    
    fig7_dat <- g_dat
    aap_inds <- which(fig7_dat$AminoAcidPairs==a)
    
    fig7_dat <- fig7_dat[aap_inds,]
    fig7_dat$MaxCodonPair <- !is.na(fig7_dat$PlotCP)
    
    #Labels for plots.
    most_abundant_cp <- unique(fig7_dat$CodonPair[which(fig7_dat$MaxCodonPair==T)])
    most_abundant_cp_seq_rscpf <- unique(fig7_dat$SequentialRSCPF[which(fig7_dat$MaxCodonPair==T)])
    most_abundant_cp_quartile <- unique(fig7_dat$CodonPairUsageBin[which(fig7_dat$MaxCodonPair==T)])
    
    #Convert bin values to quartiles.
    if (most_abundant_cp_quartile=="2") {
      most_abundant_cp_quartile <- "25-50%"
    } else if (most_abundant_cp_quartile=="3") {
      most_abundant_cp_quartile <- "50-75%"
    } else if (most_abundant_cp_quartile=="4") {
      most_abundant_cp_quartile <- "75-100%"
    }
    
    #Extract case study data.
    fig7_dat <- g_dat
    aap_inds <- which(fig7_dat$AminoAcidPairs==a)
    inds <- 1:length(fig7_dat$AminoAcidPairs)
    inds <- inds[-aap_inds]
    
    #Color points by if they are the most biased human codon-pair.
    fig7_dat$MaxCodonPair <- !is.na(fig7_dat$PlotCP)
    fig7_dat$PositionalRSCPF[inds] <- fig7_dat$MaxCodonPair[inds] <- NA
    
    nlabels <- unique(fig7_dat$DomainPanelName)
    fig7_dat$CodonPairPosition <- as.numeric(fig7_dat$CodonPairPosition)
    
    fig7_dat$CompareColor <- "black"
    fig7_dat$CompareColor[c(which(fig7_dat$Compare=="Optimal position-specific codon-pair."),
                            which(fig7_dat$Compare=="Optimal alignment- and position-specific codon-pair."),
                            which(fig7_dat$Compare=="Human codon-pair was alignment- and position-specific optimal."),
                            which(fig7_dat$Compare=="Human codon-pair was only position-specific optimal."))] <- "red"
    
    
    #######################################################
    #Plots
    
    #Code for RSCPF panel.
    fig7_rscpf <- 
      ggplot(data = fig7_dat) +
      
      geom_point(aes(x = CodonPairPosition, y = PositionalRSCPF, shape = MaxCodonPair, color = CompareColor), na.rm = T) +
      geom_hline(yintercept = 0) +
      geom_hline(yintercept = 1/aap_degen-0.05*1/aap_degen, linetype="dotted") +
      geom_hline(yintercept = 1/aap_degen+0.05*1/aap_degen, linetype="dotted") +
      geom_hline(yintercept = most_abundant_cp_seq_rscpf, linetype="dashed", linewidth = 1) +
      
      coord_cartesian(expand = F) +
      
      geom_hline(yintercept = 1) +

      scale_shape_manual(breaks = c(T, F),
                         labels = c(most_abundant_cp, "Other Synonymous \n Codon-pairs"), 
                         values = c(2, 15),
                         name = paste0("'", a, "' codon-pair"),
                         na.translate = F) +
      
      scale_color_manual(values = c("black", "red"), 
                         breaks = c("black", "red"), 
                         labels = c("No", "Yes"),
                         name = "Most biased \n position-specific \n codon-pair in \n human sequence") +
      scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.2)) +
      scale_x_continuous(limits = c(0, max(fig7_dat$CodonPairPosition)), breaks = seq(0, max(fig7_dat$CodonPairPosition), 200)) +
      theme(plot.margin = unit(c(0, 0, 0, 0), "null"),
            panel.background = element_blank(),
        axis.line.y = element_line(color = "grey", linewidth = 0.5),
        legend.margin = margin(c(0,0,0,0)),
            strip.background = element_blank(),
            strip.text = element_blank(), 
            text=element_text(size=11),
            panel.grid.major = element_line(color = "gray"),
            panel.grid.minor = element_line(color = "gray")
      ) +
    
      ylab("Relative synonymous codon pair frequency") +
      ggtitle(paste0(u, " Positional RSCPF for Amino Acid Pair ", most_abundant_cp, " (", a, ") - Quartile ", most_abundant_cp_quartile))+ 
      xlab("Codon-pair position")

    #Code for domain panel.
    fig7_domains <- ggplot(data = fig7_dat) +

      geom_tile(aes(x = CodonPairPosition, y = 1, fill = DomainPanelName), show.legend = T, na.rm = T) +
      scale_fill_manual(values = domain_hues_df[[which(names(domain_hues_df)==u_name)]][1:length(nlabels)], na.value = NA, name = "Domains") +

      scale_x_continuous(limits = c(0, max(fig7_dat$CodonPairPosition)), breaks = seq(0, max(fig7_dat$CodonPairPosition), 200)) +
      
      coord_cartesian(expand = F) +
      theme(
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        legend.margin = margin(c(0,0,0,0)),
            panel.background = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.y = element_line(color = "grey", linewidth = 0.5),
            strip.background = element_blank(),
            strip.text = element_blank(),
            text=element_text(size=11),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
      ) +
      ylab("")
    
    #Code for amino acid pair conservation panel.
    cons_plot <- 
      ggplot(data = fig7_dat) +
      geom_tile( aes(x = CodonPairPosition, y = 1, fill = Conservation), show.legend = T) +
      coord_cartesian(expand = F) +
      scale_x_continuous(limits = c(0, max(fig7_dat$CodonPairPosition)), breaks = seq(0, max(fig7_dat$CodonPairPosition), 200)) +
      scale_fill_continuous(name = "Amino acid pair \nconservation") +
      
      theme(
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y = element_line(color = "grey", linewidth = 0.5),
        legend.margin = margin(c(0,0,0,0)),
        plot.background = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),

        panel.spacing = unit(0, "lines"),
        panel.border = element_blank(),
        text=element_text(size=11))  + 
      ylab("") 
    
    if (u_name == "f8_75") {
      ptm_shape <- c(1, 17, 6)
    } else if (u_name == "f9_75") {
      ptm_shape <- c(1, 17, 6)
    } else {
      ptm_shape  <- c(1, 17)
    }
    
    #Code for PTM panel.
    ptm_plot <- 
      ggplot(data = fig7_dat) +
      geom_point(aes(x = CodonPairPosition, y = "0.00", shape = PTM), size = 4, show.legend = T) +
      geom_hline(yintercept = 0) +
      coord_cartesian(expand = F) +
      scale_x_continuous(limits = c(0, max(fig7_dat$CodonPairPosition)), breaks = seq(0, max(fig7_dat$CodonPairPosition), 200)) +
      scale_shape_manual(values = ptm_shape, 
                         na.translate = F) +
      theme(
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "null"),
            legend.title = element_blank(),
        legend.margin = margin(c(0,0,0,0)),
            axis.line.y = element_line(color = "grey", linewidth = 0.5),
            plot.background = element_blank(),
            strip.background = element_blank(),
            strip.text = element_blank(),
        axis.title.x = element_blank(),
        
        panel.background = element_blank(),
            panel.grid.minor = element_blank(),
            panel.spacing = unit(0, "lines"),
            panel.border = element_blank(),
            text=element_text(size=9))+ 
      ylab("") 
    
    #Extract grobs for each panel and prepare them for combination.
    g1grob <- ggplotGrob(fig7_rscpf)
    gdomaingrob <- ggplotGrob(fig7_domains)
    g2grob <- ggplotGrob(cons_plot)
    g3grob <- ggplotGrob(ptm_plot)
    
    g1grob_t <- ggplotGrob(fig7_rscpf+theme(legend.position = "none"))
    gdomaingrob_t <- ggplotGrob(fig7_domains+theme(legend.position = "none"))
    g2grob_t <- ggplotGrob(cons_plot+theme(legend.position = "none"))
    g3grob_t <- ggplotGrob(ptm_plot+theme(legend.position = "none"))

    legend1 <- g1grob$grobs[[which(sapply(g1grob$grobs, function(x) x$name) == "guide-box")]]
    legenddomain <- gdomaingrob$grobs[[which(sapply(gdomaingrob$grobs, function(x) x$name) == "guide-box")]]
    legend2 <- g2grob$grobs[[which(sapply(g2grob$grobs, function(x) x$name) == "guide-box")]]
    legend3 <- g3grob$grobs[[which(sapply(g3grob$grobs, function(x) x$name) == "guide-box")]]

    model_label <- legend1$grobs[[2]]
    ml <- cowplot::plot_grid(model_label)
    model_legend <- legend1$grobs[[1]]
    mp <- cowplot::plot_grid(model_legend)
    ldomain <- cowplot::plot_grid(legenddomain)
    l2 <- cowplot::plot_grid(legend2)
    l3 <- cowplot::plot_grid(legend3)
    
    #Make a legend title for domains for the multi-panel plots.
    domain_label_plot <- ggplot() + 
      annotate("text", x = 0, y = 0, size=4, label = "Domain") + 
      theme_void()
    
    domain_grob <- ggplotGrob(domain_label_plot)

    #Combine legend and plot grobs, respectively.
    legend <- cowplot::plot_grid(ml, mp, l2, domain_grob,   l3, align = "hv", nrow = 5, rel_heights = c(5,5,5,1,1))
    pgrid <- cowplot::plot_grid(g1grob_t, g2grob_t, gdomaingrob_t,  g3grob_t, align = "v", nrow = 4, rel_heights = c(20, 1, 1, 2))
    
    #Combine legend and plot grobs into a single figure.
    cowplot::plot_grid(pgrid, legend, ncol = 2, rel_widths = c(4.5, 1))
    
    setwd(pwd)
    ggsave(filename = paste0(u, " Figure 7", a, ".tiff"), width = 9, height = 5.25, units = "in")
    
  }
  
}