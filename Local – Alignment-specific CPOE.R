library(foreach)
library(doParallel)
library(seqinr)
library(plyr)
library(ggplot2)

#Character vector of gene names to analyze.
genes <- c()

for (g in genes) {
  
  #Set the working directory.
  pwd <- paste0()  
  setwd(pwd)
  
  #Load the nucleotide alignment.
load(file = paste0(g, "_c.mat"))

#Count each position's codon usages.
c_count <- foreach(c = 1:(ncol(cod_mat)/3)) %do% {
  
  plyr::count(cod_mat[,(1:3)+3*(c-1)])
  
}

#Convert the nucleotide alignment into a codon matrix
cod_mat <- foreach(c = 1:(ncol(cod_mat)/3), .combine = 'cbind') %:% foreach(r = 1:nrow(cod_mat), .combine = 'c') %do% {
  
  paste0(cod_mat[r,(1:3)+3*(c-1)], collapse = "")
  
}

#Count each position's codon-pair usages.
cp_counts <- foreach(h = 1:(ncol(cod_mat)-1), .combine = 'rbind') %do% {
  
  plyr::count(cbind(cod_mat[,h], cod_mat[,h+1]))
  
}

save(cp_counts, file = "cp.counts")

#Vectors of codon and codon-pair names.
codon_names <- toupper(seqinr::words())

cp_names <- foreach(y = 1:nrow(cp_counts), .combine = 'c') %do% {
  
  paste0(cp_counts[y, 1:2], collapse = "")
  
}

#Sum the codon-pair counts for the alignment.
cp_counts <- cp_counts[,3]
names(cp_counts) <- cp_names
cp_counts <- tapply(cp_counts, names(cp_counts), sum)  

cp_count_names <- foreach(c = 1:length(cp_counts), .combine = 'rbind') %do% {
  
  cp_string <- unlist(strsplit(names(cp_counts)[c], split = character()))
  c(paste0(cp_string[1:3], collapse = ""), paste0(cp_string[4:6], collapse = ""))
  
}

#Remove irrelevant codon-pairs.
c1_good_names <- foreach(c = 1:nrow(cp_count_names), .combine = "c") %do% {
  
  cp_count_names[c,1] %in% codon_names
  
}

c2_good_names <- foreach(c = 1:nrow(cp_count_names), .combine = "c") %do% {
  
  cp_count_names[c,2] %in% codon_names
  
}


cp_good_names <- foreach(c = 1:nrow(cp_count_names), .combine = 'c') %do% {
  
  (c1_good_names[c] == T) && (c2_good_names[c] == T)
  
}

cp_counts <- cp_counts[which(cp_good_names==T)]

#Make an empty dataframe then fill it with codon-pair counts.
cp_counts_df <- expand.grid(Codon1 = codon_names, Codon2 = codon_names)
cp_vec <- rep(NA, times = nrow(cp_counts_df))


foreach(e = 1:length(names(cp_counts))) %do% {
  
  cp <- unlist(strsplit(names(cp_counts)[e], split = character(0)))
  
  c1 <- which(cp_counts_df$Codon1 == paste0(cp[1:3], collapse = ""))
  c2 <- which(cp_counts_df$Codon2 == paste0(cp[4:6], collapse = ""))
  cp_vec[intersect(c1, c2)] <- cp_counts[e]
  
}

cp_counts_df$Counts <- cp_vec
cp_counts_df$Counts <- cp_counts_df$Counts/sum(cp_counts_df$Counts, na.rm = T)*1000
cp_counts_table <- sort(test, decreasing = T)

#Calculate codon-pair usage frequencies.
cp_freqs <- cp_counts/sum(cp_counts)

#Count each position's codon usages.
c_counts <- foreach(c = 1:length(c_count), .combine = 'rbind') %do% {
  
  c_count[[c]]
  
}

c_names <- foreach(c = 1:nrow(c_counts), .combine = 'c') %do% {
  
  paste0(c_counts[c,1:3], collapse = "")
  
}

#Sum the codon usages for the alignment and calculate their frequencies.
c_counts <- c_counts[,4]
names(c_counts) <- c_names
c_counts <- tapply(c_counts, names(c_counts), sum) 
c_freqs <- c_counts/sum(c_counts)

#Calculate alignment-specific codon-pair observed to expected ratio.
cp_oe <- foreach(c = 1:length(cp_counts), .combine = 'c') %do% {
  
  cp_split <- unlist(strsplit(names(cp_counts)[c], split = character(0)))
  c1 <- paste0(cp_split[1:3], collapse = "")
  c2 <- paste0(cp_split[4:6], collapse = "")
  
  c1_freq <- c_freqs[which(names(c_freqs)==c1)]
  c2_freq <- c_freqs[which(names(c_freqs)==c2)]
  
  oe <- cp_freqs[c]/(c1_freq*c2_freq)
  
  oe
  
}

#Make an empty dataframe then fill it with CPOE values (in amino acid pair order)
aa_names <- foreach(r = codon_names, .combine = 'c') %do% {
  
  seqinr::translate(unlist(strsplit(r, split = character(0))))
  
}

aa_order <- order(aa_names)
aa_names <- aa_names[aa_order]
aa_names <- aa_names[4:length(aa_names)]

codon_names <- codon_names[aa_order]

cp_grid <- expand.grid(Codon_1 = codon_names, Codon_2 = codon_names)
cp_grid$ObservedExpected <- NA

cp_oe <- log10(cp_oe)

foreach(p = 1:length(cp_oe)) %do% {
  
  cp <- unlist(strsplit(names(cp_oe)[p], split = character(0)))
  
  c1 <- which(cp_grid$Codon_1 == paste0(cp[1:3], collapse = ""))
  c2 <- which(cp_grid$Codon_2 == paste0(cp[4:6], collapse = ""))
  cp_grid$ObservedExpected[intersect(c1, c2)] <- cp_oe[p]
  
}

cp_grid$Codon_1 <- factor(cp_grid$Codon_1, levels = unique(cp_grid$Codon_1))
cp_grid$Codon_2 <- factor(cp_grid$Codon_2, levels = unique(cp_grid$Codon_2))

#Remove indices with stop codons.
stop_codons <- c("TAA", "TGA", "TAG")

stop_positions <- foreach(col = 1:2, .combine = 'rbind') %:% foreach(s = stop_codons, .combine = 'rbind') %do% {
  
  cp_grid[which(cp_grid[,col]==s),]
  
}

ic_obsexp_df <- cp_grid
save(ic_obsexp_df, file = "ic_obsexp.df")

}