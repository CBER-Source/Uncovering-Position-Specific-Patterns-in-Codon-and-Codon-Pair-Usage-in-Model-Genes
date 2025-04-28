library(seqinr)
library(foreach)
library(doParallel)
library(ggplot2)
library(stringr)
library(gridExtra)
library(cowplot)

#Set the parent working directory.
home_wd <- paste0()
setwd(home_wd)
source("SU functions.R")

#Character vector of gene names to analyze.
genes <- c()

for (g in genes) {
  
  #Set the working directory.
  pwd <- paste0()
  setwd(pwd)
  
  #codon pair sequence-specific permutation SU
  load(paste0(g, "_seq.SU"))
  seq_SU <- permutation_SU  
  
  #codon pair position-specific permutation SU
  load(paste0(g, "_spos.SU"))
  pos_SU <- permutation_SU
  rm("permutation_SU")
  
  #amino acid pair sequence-specific permutation SU
  load(paste0(g, "_seq_aap.SU"))
  seq_aap_SU <- permutation_aap_SU

  #amino acid pair position-specific permutation SU
  load(paste0(g, "_spos_aap.SU"))
  pos_aap_SU <- permutation_aap_SU
  rm("permutation_aap_SU")
  
  #load codon-oriented alignments
  load(paste0(g, "_c.mat"))
  aa_mat <- nt_2_aa_mat(cod_mat)
  cod_mat <- nt_2_codon_mat(cod_mat)
  
  #Remove stop.
  cod_mat <- cod_mat[,1:(ncol(cod_mat)-1)]
  aa_mat <- aa_mat[,1:(ncol(aa_mat)-1)]
  
  #Count codon usages by position.
  c_counts <- foreach(codon = 1:ncol(cod_mat), .combine = 'c') %do% {
    
    counts <- plyr::count(cod_mat[,codon])
    codon_names <- counts[,1]
    counts <- counts[,2]
    names(counts) <- codon_names
    return(counts)
    
  }
  
  #Sum positional usages to get values for the entire alignment.
  c_counts <- tapply(c_counts, names(c_counts), sum)
  
  #Count codon-pair usages by position.
  cp_counts <- foreach(h = 1:(ncol(cod_mat)-1), .combine = 'c') %do% {
    
    counts <- plyr::count(cbind(cod_mat[,h], cod_mat[,h+1]))
    codon_pair_names <- paste0(counts[,1], counts[,2])
    counts <- counts[,3]
    names(counts) <- codon_pair_names
    return(counts)
    
  }
  
  #Sum positional usages to get values for the entire alignment.
  cp_counts <- tapply(cp_counts, names(cp_counts), sum)
  
  #amino acid pair SU
  aap_ref_SU <- aap_SU(aa_mat)
  save(aap_ref_SU, file = "aap_ref.su") 

  #codon-pair SU
  ref_SU <- cp_SU(cod_mat)
  save(ref_SU, file = "ref.su")
  
  #Calculate p-values for original SU vs seq/pos distributions (alignment/position-)
  seq_pvals <- foreach(cp = 1:length(ref_SU), .combine = 'rbind') %do% {
    
    top5_range <- quantile(seq_SU[,cp], c(0.96, 1), na.rm = T) 
    
    if (ref_SU[cp]>=top5_range[1]) {
      
      if (ref_SU[cp]>top5_range[2]) {
        return(c(ref_SU[cp], "Original su > 100% (sequence-specific CP SU)", top5_range))
      } else if (ref_SU[cp]<=top5_range[2]) {
        return(c(ref_SU[cp], "96% < Original su < 100% (sequence-specific CP SU)", top5_range))
      } 
      
    } else {
      return(c(ref_SU[cp], "Original SU < 95% (sequence-specific CP SU)", top5_range))
    }
    
  }
  
  seq_pvals <- cbind(1:length(ref_SU), seq_pvals)
  colnames(seq_pvals) <- c("Codon-pair position", 
                           "Original CP SU",
                           "Original CP SU vs sequence-specific CP SU",
                           "95% sequence-specific CP SU",
                           "100% sequence-specific CP SU")
  
  pos_pvals <- foreach(cp = 1:length(ref_SU), .combine = 'rbind') %do% {
    
    top5_range <- quantile(pos_SU[,cp], c(0.96, 1), na.rm = T) 
    
    
    if (ref_SU[cp]>=top5_range[1]) {
      
      if (ref_SU[cp]>top5_range[2]) {
        return(c(ref_SU[cp], "Original su > 100% (position-specific CP SU)", top5_range))
      } else if (ref_SU[cp]<=top5_range[2]) {
        return(c(ref_SU[cp], "96% < Original su < 100% (position-specific CP SU)", top5_range))
      } 
      
    } else {
      return(c(ref_SU[cp], "Original SU < 95% (position-specific CP SU)", top5_range))
    }
    
  }

  pos_pvals <- cbind(1:length(ref_SU), pos_pvals)
  colnames(pos_pvals) <- c("Codon-pair position", 
                           "Original CP SU",
                           "Original CP SU vs position-specific CP SU",
                           "95% position-specific CP SU",
                           "100% position-specific CP SU")
  
  pvals <- cbind(seq_pvals, pos_pvals[,3:5])

  #Summarize the p-value results and output as a .csv.
  pos_pval_summary <- plyr::count(pvals[,6])
  pos_pval_percents <- pos_pval_summary[,2]/sum(pos_pval_summary[,2])*100
  pos_pval_summary$percents <- pos_pval_percents
  pos_pval_summary$method <- "synonymous position-specific"
  
  seq_pval_summary <- plyr::count(pvals[,3])
  seq_pval_percents <- seq_pval_summary[,2]/sum(seq_pval_summary[,2])*100
  seq_pval_summary$percents <- seq_pval_percents
  seq_pval_summary$method <- "sequence-specific"
  
  pval_summary <- rbind(pos_pval_summary, old_pval_summary, seq_pval_summary)
  write.csv(pval_summary, paste0(g, " 96-100 summary.csv"))
  
  #Calculate average SU values for permuted alignments.
  seq_avg_su <- colMeans(seq_SU)
  pos_avg_su <- colMeans(pos_SU)
  seq_aap_avg_su <- colMeans(seq_aap_SU)
  pos_aap_avg_su <- colMeans(pos_aap_SU)

  seq_length <- length(ref_SU)
  
  #Set the order of factors for plotting.
  model_order <- c("Original Codon-pair",
                   "Sequential Method Codon-pair", 
                   "Positional Method Codon-pair", 
                   "Original Amino Acid Pair",
                   "Sequential Method Amino Acid Pair",
                   "Positional Method Amino Acid Pair"
                   )
  
  model_colors <- c("black", 
                    "blue", 
                    "red", 
                    "black",
                    "blue",
                    "red"
                    )
  
  #Make a dataframe for plotting.
  su_dt <- data.frame(Position = 1:seq_length, Model = rep(model_order, 
                                                           each = seq_length), 
                      SU = c(ref_SU, 
                             seq_avg_su, 
                             pos_avg_su, 
                             aap_ref_SU,  
                             seq_aap_avg_su,
                             pos_aap_avg_su
                             
                             ))
  
  #Read consurf amino acid conservation 
  consurf_aa <- as.data.frame(readxl::read_excel('consurf_grades-aa.xlsx', skip = 27, col_names = T))
  
  #Estimate amino acid pair conservation from amino acid conservation.
  aap_conservation <- foreach(a = 1:(nrow(consurf_aa)-1), .combine = 'c') %do% {
    
    mean(c(as.numeric(consurf_aa[a,4]), as.numeric(consurf_aa[a+1,4])), na.rm = T)
    
  }
  
  cons_dt <- data.frame(Position = 1:length(aap_conservation), Conservation = aap_conservation)

  #Read gene feature information.
  gff <- read.delim(paste0(g, " - domains - uniprot - 3-27-24.gff"), header=F, comment.char="#")

  #extract post translational modification data
  ptm_types <- c("Disulfide bond", "Glycosylation", "Modified residue")
  ptm_df <- foreach(p = ptm_types, .combine = 'rbind') %do% {
    
    gff[which(gff[,3]==p),]
    
  }
  
  ptm_vec <- foreach(r = 1:nrow(ptm_df), .combine = 'c') %do% {
    
    if (ptm_df[r,3]=="Disulfide bond") {
      c(ptm_df[r,4],ptm_df[r,5])
    } else {
      c(ptm_df[r,4]:ptm_df[r,5])
    }
    
  }
  
  #Add labels to positions with no ptms.
  no_ptm_inds <- setdiff(1:(ncol(cod_mat)-1), unique(ptm_vec))
  no_ptm_inds <- unname(split(no_ptm_inds, cumsum(c(TRUE, diff(no_ptm_inds)!=1))))
  
  no_ptm_df <- foreach(l = 1:length(no_ptm_inds), .combine = 'rbind') %do% {
    
    no_name <- paste0("No ptm - ", l)
    cbind(no_ptm_inds[[l]], no_name)
    
  }
  
  no_ptm_df <- data.frame("Position"=as.numeric(no_ptm_df[,1]), 
                          "PTM"=NA)
  
  #Add labels to positions with ptms.
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
  
  #Combine non-ptm and ptm dataframes.
  ptm_df <- data.frame("Position"=as.numeric(ptm_df[,1]), "PTM"=ptm_df[,2])
  ptm_df <- merge(ptm_df, no_ptm_df, all = T)
  
  
  #Extract domain information.
  domain_rows <- which(gff[,3]=="Domain")
  gff <- gff[domain_rows,]
  
  domain_vec <- foreach(r = 1:nrow(gff), .combine = 'c') %do% {
    
    c(gff[r,4]:gff[r,5])
    
  }
  
  #Label positions with no domain.
  no_domain_inds <- setdiff(1:(ncol(cod_mat)-1), domain_vec)
  no_domain_inds <- unname(split(no_domain_inds, cumsum(c(TRUE, diff(no_domain_inds)!=1))))
  
  no_domain_df <- foreach(l = 1:length(no_domain_inds), .combine = 'rbind') %do% {
    
    no_name <- paste0("No domain - ", l)
    cbind(no_domain_inds[[l]], no_name)
    
  }
  
  #Label positions in a domain.
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
  
  #Add the results to a dataframe.
  domain_df <- data.frame(Position = as.numeric(domain_df[,1]), Domain = domain_df[,2])
  no_domain_df <- data.frame(Position = as.numeric(no_domain_df[,1]), Domain = no_domain_df[,2])
  
   if (max(domain_df$Position)>(ncol(cod_mat)-1)) {
    stop_ind <- which(domain_df$Position==max(domain_df$Position))
    domain_df <- domain_df[-stop_ind,]
  }
  
  #Combine positions with and without domains.
  domain_df <- rbind(domain_df, no_domain_df)

  #Add domain informaiton to SU data.
  su_dt <- merge(su_dt, domain_df, all = T)
  dnames <- unique(su_dt$Domain)
  
  #Prepare domain data for plotting.
  su_dt$DomainFactor <- factor(su_dt$Domain, unique(su_dt$Domain))
  dnumber <- unique(su_dt$Domain)
  dpair <- rbind(dnames, dnumber)
  DomainPanelName <- foreach(d = su_dt$Domain, .combine = 'c') %do% {
    
  CheckForNoDomain <- grepl("No domain*", d)  
  
  if (CheckForNoDomain==T) {
    paste("No domain")
  } else {
    paste(d)
  }
  }
  
  su_dt$DomainPanelName <- factor(DomainPanelName, unique(DomainPanelName), ordered = T)
  nlabels <- levels(su_dt$DomainPanelName)

  #Add ptms and conservation to SU data.
  su_dt <- merge(su_dt, ptm_df, all = T)
  su_dt <- merge(su_dt, cons_dt, all = T)

  #Labels for different types of SU to be used in plotting.
  su_dt$Model <- factor(su_dt$Model, levels = model_order, ordered = T)
  su_dt$SUvariable <- su_dt$SUdata <- NA
  aap_inds <- grep("Amino", su_dt$Model)
  su_dt$SUvariable[aap_inds] <- "Amino Acid Pair"
  cp_inds <- grep("Codon", su_dt$Model)
  su_dt$SUvariable[cp_inds] <- "Codon-pair"
  
  su_dt$SUvariable <- factor(x = su_dt$SUvariable,
                             levels = c("Codon-pair",
                               "Amino Acid Pair"
                                        ),
                             ordered = T)
  
  
  original_inds <- grep("Original", su_dt$Model)
  su_dt$SUdata[original_inds] <- "Original"
  sequential_inds <- grep("Sequential", su_dt$Model)
  su_dt$SUdata[sequential_inds] <- "Sequential"
  positional_inds <- grep("Positional", su_dt$Model)
  su_dt$SUdata[positional_inds] <- "Positional"
  su_dt$SUdata <- factor(x = su_dt$SUdata,
                         levels = c("Original",  
                                    "Positional",  
                                    "Sequential" ),
                         ordered = T)
  
  
  model_names <- as_labeller(c("Original Codon-pair" = "Original",
                               "Sequential Method Codon-pair" = "Sequence-specific", 
                               "Positional Method Codon-pair" = "Position-specific", 
                               "Original Amino Acid Pair" = "Original", 
                               "Positional Method Amino Acid Pair" = "Position-specific",
                               "Sequential Method Amino Acid Pair" = "Sequence-specific"))
  
  
  #Make a legend title for domains for the multi-panel plots.
  domain_label_plot <- ggplot() + 
    annotate("text", x = 0, y = 0, size=8, label = "Domain") + 
    theme_void()
  
  domain_grob <- ggplotGrob(domain_label_plot)
  
  ####
  
  #Convert file titles to plot titles.
  if (g == "f8_75") {
    plot_title <- "F8"
  } else if (g == "f9_75") {
    plot_title <- "F9"
  } else if (g == "ADAMTS13_75") {
    plot_title <- "ADAMTS13"
  } else if (g == "vwf_75") {
    plot_title <- "VWF"
  }
  
  ####
  #Normal SU with domain, conservation, and ptm
  normal_su_plot <- su_dt
  
  #Code for SU panel.
  su_plot <- 
    ggplot(normal_su_plot) + 

    geom_point(aes(x=Position, y=SU, color = SUdata), show.legend = T) +
    ylab("Symmetric Uncertainty") +
    
    xlab("Codon Pair Position") +
    coord_cartesian(expand = F) +
    scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.05)) +
    scale_x_continuous(limits = c(0, max(su_dt$Position)), breaks = seq(0, max(su_dt$Position), 200)) +
    scale_color_manual(labels = c("Original",   "Position-specific", "Sequence-specific"), values = model_colors[c(1,5,4)], na.value = NA, name = "Alignment type") +
    
    facet_wrap(~SUvariable, ncol = 1, nrow = 2) +

    ggtitle(paste0(plot_title, " Permutation Test Results")) +
    

    theme(                    axis.line.y = element_line(color = "grey", linewidth = 0.5),
plot.margin = unit(c(0, 0, 0, 0), "null"),plot.title = element_text(hjust = 0.5),
          axis.line = element_line(color = "black", linewidth = 0.5),
          plot.background = element_blank(),
          strip.background = element_blank(),
          panel.spacing = unit(0, "lines"),
          panel.border = element_blank(),
          panel.background = element_blank(),
          text=element_text(size=22),
          strip.text = element_text(size=22),
          panel.grid.major = element_line(color = "gray"),
          panel.grid.minor = element_blank()) 
  
  #Code for domain panel.
  domain_plot <-     
    ggplot(normal_su_plot) + 
    geom_tile( aes(x = Position, y = 1, fill = DomainPanelName), show.legend = F) +
    scale_fill_manual(labels = nlabels, values = domain_hues_df[[which(names(domain_hues_df)==g)]][1:length(nlabels)]) +
      coord_cartesian(expand = F) +
    scale_x_continuous(limits = c(0, max(su_dt$Position)), breaks = seq(0, max(su_dt$Position), 200)) +
      theme(                    axis.line.y = element_line(color = "grey", linewidth = 0.5),
plot.margin = unit(c(0, 0, 0, 0), "null"),axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
axis.title.x = element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.line = element_line(color = "black", linewidth = 0.5),
            plot.background = element_blank(),
            strip.background = element_blank(),
            strip.text.x = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.spacing = unit(0, "lines"),
            panel.border = element_blank(),
            text=element_text(size=22))  + 
      ylab("") 
  
  #Code for amino acid pair conservation panel.
  cons_plot <- 
    ggplot(data = normal_su_plot) +
    geom_tile( aes(x = Position, y = 1, fill = Conservation), show.legend = T) +
    coord_cartesian(expand = F) +
    scale_x_continuous(limits = c(0, max(su_dt$Position)), breaks = seq(0, max(su_dt$Position), 200)) +
    scale_fill_continuous(name = "Amino acid pair \nconservation") +
    theme(                    axis.line.y = element_line(color = "grey", linewidth = 0.5),
plot.margin = unit(c(0, 0, 0, 0), "null"),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
axis.title.x = element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.line = element_line(color = "black", linewidth = 0.5),
          plot.background = element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.spacing = unit(0, "lines"),
          panel.border = element_blank(),
          text=element_text(size=22))  + 
    ylab("") 
  
  #Custom point shapes for PTM panels, because the human reference genes have different numbers of PTMs.
  if (g == "f8_75") {
    ptm_shape <- c(1, 17, 6)
  } else if (g == "f9_75") {
    ptm_shape <- c(1, 17, 6)
  } else if (grepl("samp", g)) {
    ptm_shape <- c(1, 17, 6)
  } else {
    ptm_shape  <- c(1, 17)
  }
  
  #Code for PTM panel.
  ptm_plot <- 
    ggplot(data = normal_su_plot) +
    geom_point(aes(x = Position, y = "0.00", shape = PTM), size = 4, show.legend = T) +
    geom_hline(yintercept = 0) +
    coord_cartesian(expand = F) +
    scale_x_continuous(limits = c(0, max(su_dt$Position)), breaks = seq(0, max(su_dt$Position), 200)) +
    scale_shape_manual(values = ptm_shape, 
                       na.translate = F) +
    theme(                    axis.line.y = element_line(color = "grey", linewidth = 0.5),
plot.margin = unit(c(0, 0, 0, 0), "null"),axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.title = element_blank(),
          axis.title.x = element_blank(),
          plot.background = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.spacing = unit(0, "lines"),
          panel.border = element_blank(),
          text=element_text(size=22))+ 
    ylab("") 
  
  #Extract grobs for each panel and prepare them for combination.
  g1grob <- ggplotGrob(su_plot)
  gdomaingrob <- ggplotGrob(domain_plot)
  g2grob <- ggplotGrob(cons_plot)
  g3grob <- ggplotGrob(ptm_plot)
  
  g1grob_t <- ggplotGrob(su_plot+theme(legend.position = "none"))
  gdomaingrob_t <- ggplotGrob(domain_plot+theme(legend.position = "none"))
  g2grob_t <- ggplotGrob(cons_plot+theme(legend.position = "none"))
  g3grob_t <- ggplotGrob(ptm_plot+theme(legend.position = "none"))

  legend1 <- g1grob$grobs[[which(sapply(g1grob$grobs, function(x) x$name) == "guide-box")]]
  legend2 <- g2grob$grobs[[which(sapply(g2grob$grobs, function(x) x$name) == "guide-box")]]
  legend3 <- g3grob$grobs[[which(sapply(g3grob$grobs, function(x) x$name) == "guide-box")]]
  
  l1 <- cowplot::plot_grid(legend1)
  l2 <- cowplot::plot_grid(legend2)
  l3 <- cowplot::plot_grid(legend3)

  #Combine legend and plot grobs, respectively.
  legend <- cowplot::plot_grid(l1, l2, domain_grob,   l3, align = "v", nrow = 4, rel_heights = c(4,5,1,1))
  pgrid <- cowplot::plot_grid(g1grob_t, g2grob_t, gdomaingrob_t,  g3grob_t, align = "v", nrow = 4, rel_heights = c(20, 1, 1, 2))

  #Combine legend and plot grobs into a single figure.
  cowplot::plot_grid(pgrid, legend, ncol = 2, rel_widths = c(5, 1))
  ggsave(filename = paste(g, "Figure 6 - SU.tiff"), width = 17.5, height = 10.5, units = "in")

  #Function to calculate a vector's mode.
  Mode <- function(x, na.rm = FALSE) {
    if(na.rm){
      x = x[!is.na(x)]
    }
    
    ux <- unique(x)
    return(ux[which.max(tabulate(match(x, ux)))])
  }
  
  #Summary statistics of the gene's SU. 
  SU_stat_table <- foreach(u = unique(normal_su_plot$Model), .combine = 'rbind') %do% {
    
    inds <- which(normal_su_plot$Model==u)
    dat <- normal_su_plot$SU[inds]
    
    c(g,
      paste(u),
      mean(dat, na.rm = T),
      median(dat, na.rm = T),
      Mode(dat, na.rm = T)
    )
    
  }
  
  colnames(SU_stat_table) <- c("Gene",
                               "Model",
                               "SU Mean",
                               "SU Median",
                               "SU Mode")
  
  write.csv(SU_stat_table, file = paste0("SU statistic table for Supplemental File 8 - ", g, ".csv"))
  
  #Summary statistics of the gene's SU by domain.
  domain_stat_dat <- normal_su_plot[which(normal_su_plot$Model=="Original Codon-pair"),]
  
  domain_stat_table <- foreach(u = unique(domain_stat_dat$Domain), .combine = 'rbind') %do% {
    
    inds <- which(domain_stat_dat$Domain==u)
    dat <- domain_stat_dat$SU[inds]
    cons_dat <- domain_stat_dat$Conservation[inds]
    
    c(g,
      paste(u),
      paste(unique(domain_stat_dat$DomainPanelName[inds])),
      mean(dat, na.rm = T),
      median(dat, na.rm = T),
      Mode(dat, na.rm = T),
      mean(cons_dat, na.rm = T),
      median(cons_dat, na.rm = T),
      Mode(cons_dat, na.rm = T)
    )
    
  }
  
  colnames(domain_stat_table) <- c("Gene",
                               "Domain",
                               "Domain Group",
                               "SU Mean",
                               "SU Median",
                               "SU Mode",
                               "Amino acid pair conservation Mean",
                               "Amino acid pair conservation Median",
                               "Amino acid pair conservation Mode"
                               )
  
  write.csv(domain_stat_table, file = paste0("SU statistic table for Supplemental File 8 - ", g, "-domains.csv"))
  
  #Modify column names for easier interaction during regression analysis.
  colnames(domain_stat_table) <- c("Gene",
                                   "Domain",
                                   "DomainGroup",
                                   "SUMean",
                                   "SUMedian",
                                   "SUMode",
                                   "AAPconsMean",
                                   "AAPconsMedian",
                                   "AAPconsMode")
  
  #Format data for regression.
  domain_stat_table <- as.data.frame(domain_stat_table)
  domain_stat_table$Domain <- factor(domain_stat_table$Domain, levels = unique(domain_stat_table$Domain), ordered = T)
  
  domain_stat_table$SUMedian <- as.numeric(domain_stat_table$SUMedian)
  domain_stat_table$SUMean <- as.numeric(domain_stat_table$SUMean)
  domain_stat_table$AAPconsMedian <- as.numeric(domain_stat_table$AAPconsMedian)
  domain_stat_table$AAPconsMean <- as.numeric(domain_stat_table$AAPconsMean)
  
  domain_stat_table$DomainGroup <- factor(domain_stat_table$DomainGroup, levels = unique(domain_stat_table$DomainGroup), ordered = T)
  nlabels_domain_su <- levels(domain_stat_table$DomainGroup)
  
  #SU and amino acid pair conservation regression by domain mean/median.
  model <- lm(formula = SUMean ~ AAPconsMean, 
              data = domain_stat_table)
  
  sink(paste0(g, " SU vs conservation by domain - mean.txt"))
  print(summary(model))
  sink() 
  
  model <- lm(formula = SUMedian ~ AAPconsMedian, 
              data = domain_stat_table)
  
  sink(paste0(g, " SU vs conservation by domain - median.txt"))
  print(summary(model))
  sink() 
  
}