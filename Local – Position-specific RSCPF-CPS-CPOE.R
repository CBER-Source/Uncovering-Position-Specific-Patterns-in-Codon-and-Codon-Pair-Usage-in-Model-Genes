library(foreach)
library(doParallel)

#Set up a parallel backend for the codon pair count function.
#Choose an integer for the number of cores to use during parallelization.
clustnum <- 
co <- clustnum
cl <- makeForkCluster(nnodes = co)
registerDoParallel(cl, clustnum)

#Character vector of gene names to analyze.
genes <- c()

for (g in genes) {
  
  #Set the working directory.
  pwd <- paste0()  
  setwd(pwd)
  
  #Load the codon matrix
  load(file = paste0(g, "_codon_alignment.mat"))
  
  cod_mat <- cod_mat[,1:(ncol(cod_mat)-1)]
  
  #Count codon-pair usages by position.
  aln_counts <- foreach(u = 1:(ncol(cod_mat)-1), .combine = 'rbind') %do% {
    
    count <-  cbind(cod_mat[,u], cod_mat[,u+1])
    
    #Remove irrelevant codon-pairs.
    bad_inds <- c(grep("-", count[,1]),
                  grep("-", count[,2]),
                  grep("N", count[,1]),
                  grep("N", count[,2]),
                  grep("R", count[,1]),
                  grep("R", count[,2]),
                  grep("Y", count[,1]),
                  grep("Y", count[,2]))
    
    bad_inds <- unique(bad_inds)
    
    if (length(bad_inds) == nrow(count)) {
      NA
    } else if (length(bad_inds) == (nrow(count)-1)) {
      NA
    } else {
      
      if (length(bad_inds) > 0) {
        count <- count[-bad_inds,]
      }
      

      cp <- plyr::count(count)
      
      cp <- cbind(cp, cp[,3]/sum(cp[,3]))
      colnames(cp) <- c("Codon1",
                        "Codon2",
                        "NumberObserved",
                        "ObservedFrequency")
      cp
      
    }}
  
  #Sum the position-specific counts for the alignment.
  m1_cp_names <- paste0(aln_counts[,1], aln_counts[,2])
  m1_cp_counts <- aln_counts[,3]
  names(m1_cp_counts) <- m1_cp_names
  m1_cp_counts <- tapply(m1_cp_counts, names(m1_cp_counts), sum)
  
  #Position-specific codon, codon-pair, amino acid, and amino acid pair counts for codon-pair bias metric calculations.
  cpbm_counts <- foreach(u = 1:(ncol(cod_mat)-1)) %do% {
    
    count <-  cbind(cod_mat[,u], cod_mat[,u+1])
    
    #Remove irrelevant codon-pairs.
    bad_inds <- c(grep("-", count[,1]),
                  grep("-", count[,2]),
                  grep("N", count[,1]),
                  grep("N", count[,2]),
                  grep("R", count[,1]),
                  grep("R", count[,2]),
                  grep("Y", count[,1]),
                  grep("Y", count[,2]))
    
    bad_inds <- unique(bad_inds)
    
    if (length(bad_inds) == nrow(count)) {
      NA
    } else if (length(bad_inds) == (nrow(count)-1)) {
      NA
    } else {
      
      if (length(bad_inds) > 0) {
        count <- count[-bad_inds,]
      }
      
      ##########
      
      #Amino acid and amino acid pair usage counts and frequencies.
      a1 <- seqinr::translate(unlist(strsplit(count[,1], split = character(0))))
      a2 <- seqinr::translate(unlist(strsplit(count[,2], split = character(0))))
      a12 <- cbind(a1, a2)
      
      a1 <- plyr::count(a1)
      a2 <- plyr::count(a2)
      a12 <- plyr::count(a12)
      
      a1 <- cbind(a1, a1[,2]/sum(a1[,2]))
      colnames(a1) <- c("AminoAcid 1",
                        "NumberObserved",
                        "ObservedFrequency")
      
      a2 <- cbind(a2, a2[,2]/sum(a2[,2]))
      colnames(a2) <- c("AminoAcid 2",
                        "NumberObserved",
                        "ObservedFrequency")
      
      a12 <- cbind(a12, a12[,3]/sum(a12[,3]))
      colnames(a12) <- c("AminoAcid 1",
                         "AminoAcid 2",
                         "NumberObserved",
                         "ObservedFrequency")
      
      #Codon and codon-pair usage counts and frequencies.
      cp <- plyr::count(count)
      c1 <- plyr::count(count[,1])
      c2 <- plyr::count(count[,2])
      
      
      c1 <- cbind(c1, c1[,2]/sum(c1[,2]))
      colnames(c1) <- c("Codon1",
                        "NumberObserved",
                        "ObservedFrequency")
      
      c2 <- cbind(c2, c2[,2]/sum(c2[,2]))
      colnames(c2) <- c("Codon2",
                        "NumberObserved",
                        "ObservedFrequency")
      
      cp <- cbind(cp, cp[,3]/sum(cp[,3]))
      colnames(cp) <- c("Codon1",
                        "Codon2",
                        "NumberObserved",
                        "ObservedFrequency")
      
      list("CodonPairs" = cp, 
           "Codon1" = c1, 
           "Codon2" = c2,
           "AAPairs" = a12,
           "AA1" = a1,
           "AA2" = a2)
      
      ##########
    } #end bad inds else
  } #end cpbm loop
  
  #Calculate position-specific RSCPF, CPS, and CPOE.
  cpbm_calc <- foreach(l = 1:length(cpbm_counts), .combine = 'rbind') %do% {
    
    na_check <- is.na(cpbm_counts[[l]])
    
    #If all codon-pairs are removed during filtration in the last step return NA.
    if (T %in% na_check) {
      NA
    } else {
      
      cpbm <- foreach(m = 1:nrow(cpbm_counts[[l]]$CodonPairs), .combine = 'rbind', .multicombine = T) %do% {
        
        loop_ref <- paste(cpbm_counts[[l]]$CodonPairs[m,])
        loop_aa <- seqinr::translate(unlist(strsplit(loop_ref[1:2], split = character(0))))
        
        c1 <- cpbm_counts[[l]]$Codon1[which(cpbm_counts[[l]]$Codon1[,1]==loop_ref[1]),3]
        c2 <- cpbm_counts[[l]]$Codon2[which(cpbm_counts[[l]]$Codon2[,1]==loop_ref[2]),3]
        c1c2 <- c1*c2
        
        a1 <- cpbm_counts[[l]]$AA1[which(cpbm_counts[[l]]$AA1[,1]==loop_aa[1]),3]
        a2 <- cpbm_counts[[l]]$AA2[which(cpbm_counts[[l]]$AA2[,1]==loop_aa[2]),3]
        a1a2 <- a1*a2
        
        
        a1_ref <- which(cpbm_counts[[l]]$AAPairs[,1]==loop_aa[1])
        a2_ref <- which(cpbm_counts[[l]]$AAPairs[,2]==loop_aa[2])
        a12_ref <- intersect(a1_ref, a2_ref)
        
        a12 <- cpbm_counts[[l]]$AAPairs[a12_ref,4]
        cp <- as.numeric(loop_ref[4])
        cp_count <- as.numeric(loop_ref[3])
        
        #Calculate and log10 normalize CPS.
        loop_cps <- cp/((c1c2/a1a2)*a12)
        log_loop_cps <- log10(loop_cps)
        
        #Calculate RSCPF.
        loop_rscpu <- cp/a12
        
        #Calculate and log10 normalize CPOE.
        loop_oe <- cp/c1c2
        log_loop_oe <- log10(loop_oe)
        
        c(cp_count, log_loop_cps, loop_rscpu, log_loop_oe)
        
      } #end good cps_rscpu loop
      
      #Different methods to label the results depending on the size of the results.
      if (nrow(cpbm_counts[[l]]$CodonPairs)>1) {
        cpbm <- cbind(l,
                      paste0(cpbm_counts[[l]]$CodonPairs[,1], cpbm_counts[[l]]$CodonPairs[,2]),
                      cpbm)
        colnames(cpbm) <- c("CodonPairPosition", "CodonPair", "CodonPairsObserved", "CPS", "RSCPU", "ObsExp")
        
      } else if (nrow(cpbm_counts[[l]]$CodonPairs)==1) {
        cpbm <- c(l,
                      paste0(cpbm_counts[[l]]$CodonPairs[,1], cpbm_counts[[l]]$CodonPairs[,2]),
                      cpbm)
        names(cpbm) <- c("CodonPairPosition", "CodonPair", "CodonPairsObserved", "CPS", "RSCPU", "ObsExp")
        
      } else {
        names(cpbm) <- c("CodonPairPosition", "CodonPair", "CodonPairsObserved", "CPS", "RSCPU", "ObsExp")
      }
      
      return(cpbm)
      
    } #end is.na else
    
  } #end cpbm_calc foreach
  
  save(cpbm_counts, file = paste0(g, "_cpbm.counts"))
  save(cpbm_calc, file = paste0(g, "_cpbm.calc"))
  
} #end gene loop

#Close parallel backend.
stopCluster(cl)
