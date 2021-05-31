## Miscellaneous functions used for KARMA1 master thesis pipeline


# Function for counting unique proteins
get_unique_prot <- function(binfo) {
  # Count unique proteins
  # binfo is the binder info data frame, subset as pleased
  unpro <- unique(binfo$ENSG.ID)
  unpro <- unpro[unpro != ""]  # Remove empty ENSG.IDs
  return(unpro)
}


# Plotting box charts to inspect variation, used in normalisation.Rmd
plot_96well_box <- function(m, s, norm=FALSE) {
  # Function for plotting box plots to look at variation
  # m is log-transformed MFI data, s is sample info, set norm is TRUE if data is normalised
  if (norm) {
    plt_name1 <- "384-well plate 1 (normalised)"
    plt_name2 <- "384-well plate 2 (normalised)"
  } else {
    plt_name1 <- "384-well plate 1"
    plt_name2 <- "384-well plate 2"
  }
  
  melty_df1 <- data.frame(m[s$plate384 == 1, ], 
                          "plate"=s$plate96[s$plate384 == 1])
  melty_df2 <- data.frame(m[s$plate384 == 2, ],
                          "plate"=s$plate96[s$plate384 == 2])
  
  melt96_1 <- melt(melty_df1,
                   id.vars="plate",
                   measure.vals=subset(melty_df1, select=-"plate"))
  
  melt96_2 <- melt(melty_df2,
                   id.vars="plate",
                   measure.vals=subset(melty_df2, select=-"plate"))
  
  par(mfrow=c(1,2))
  bp1 <- boxplot(value~plate, data=melt96_1,
                 main=plt_name1,
                 xlab="96-well plate", ylab="log(MFI)",
                 ylim=c(4, 9.5))
  abline(h=bp1$stats[3, 1], lty="dashed", col="red")
  
  bp2 <- boxplot(value~plate, data=melt96_2,
                 main=plt_name2, 
                 xlab="96-well plate", ylab="log(MFI)",
                 ylim=c(4, 9.5))
  abline(h=bp2$stats[3, 1], lty="dashed", col="red")
}


# UMAP on cases/controls/doubles, used in normalisation.Rmd and lin_normalisation.Rmd
umap_casecontrol <- function(m, s) {
  # Function for performing UMAP on cases, controls and doubles ("Other" samples)
  # m is the mfi data frame to do the UMAP on with sample names as rownames, s is 
  # the corresponding sample info, containing information about the classes of the 
  # samples (case, control etc)
  s_umap <- s[s$Class == "Case" | s$Class == "Control" | s$Class == "Other", ]
  m_umap <- m[rownames(s_umap), ]
  
  the_umap <- umap(m_umap[, -1], random_state=123)
  return(the_umap)
}


# Plotting several scatter plots of UMAP-reduced data, used in conjunction with umap_casecontrol()
plot_umap <- function(umap_mfi, dimred_s, new_clin, printall=T) {
  # umap_mfi is a umap object, dimred_s is the corresponding sample info.
  # Plot UMAP values and colour by different variables
  umap_data <- data.frame("UMAP1"=umap_mfi$layout[, 1],
                          "UMAP2"=umap_mfi$layout[, 2],
                          "sampname"=dimred_s$Sample.id, 
                          "Hospital"=dimred_s$hospital,
                          "Age"=dimred_s$Sampling.age,
                          "BMI"=dimred_s$bmi,
                          "Entry_date"=as.Date(dimred_s$entry.date),
                          "Sampling_date"=as.Date(dimred_s$Sampling.date),
                          "Class"=dimred_s$Class,
                          "tubelabs"=dimred_s$Tube.label,
                          "repl"=dimred_s$replicate,
                          "postmenop"=dimred_s$postmenop)
  
  # Add info to connect replicates with each other and doubles with each other
  umap_data$repl.groups = umap_data$other.groups <- NA
  umap_data[umap_data$repl == 1 | umap_data$repl == 2, "repl.groups"] <-
    umap_data[umap_data$repl == 1 | umap_data$repl == 2, "sampname"]
  umap_data[umap_data$Class == "Other", "other.groups"] <-
    umap_data[umap_data$Class == "Other", "tubelabs"]
  
  q <- ggplot(data = umap_data, mapping = aes(x = UMAP1, y = UMAP2))
  
  # By hospital
  hosp_q <- q + geom_point(mapping = aes(colour = Hospital)) +
          labs(title = "UMAP coloured by hospital") +
          scale_colour_brewer(palette = "Accent", na.value = "red")
  
  # By age
  age_q <- q + geom_point(mapping = aes(colour = Age)) +
          scale_colour_viridis_c(na.value = "red") +
          labs(title = "UMAP coloured by age",
               colour = "Age (yrs)")
  
  # By BMI
  bmi_q <- q + geom_point(mapping = aes(colour = log(BMI))) +
          scale_colour_viridis_c(na.value = "transparent") +
          labs(title = "UMAP coloured by BMI",
               colour = "log(BMI)")
  

  # By study entry date
  date_labs <- pretty(umap_data$Entry_date)
  entrydate_q <- q + geom_point(mapping = aes(colour = as.numeric(Entry_date))) +
          scale_colour_viridis_c(breaks = as.numeric(date_labs),
                                 labels = date_labs,
                                 na.value = "transparent") +
          labs(title = "UMAP coloured by study entry date",
               colour = "Date")
  
  # If using new clinical data, plot the density distributions too
  if (new_clin) {
    umap_data$pdens <- dimred_s$percent.dens
    umap_data$adens <- dimred_s$area.dens
    q <- ggplot(data = umap_data, mapping = aes(x = UMAP1, y = UMAP2))
    pdens_q <- q + geom_point(mapping = aes(colour = pdens)) +
            scale_colour_viridis_c(na.value = "transparent") +
            labs(title = "UMAP coloured by percent density",
                 colour = "Percent density")

    adens_q <- q + geom_point(mapping = aes(colour = adens)) + 
            scale_colour_viridis_c(na.value = "transparent") + 
            labs(title = "UMAP coloured by absolute dense area",
                 colour = "Absolute dense area")
  }
  
  # Plot distribution of premenopausal/postmenopausal
  postmenop_q <- q + geom_point(mapping = aes(colour = postmenop)) + 
          scale_colour_brewer(palette = "Accent", na.value = "transparent") + 
          labs(title = "UMAP coloured by postmenopausal status",
               colour = "Postmenopausal (1) or not (0)")
  
  # Also by class since distracting to use as shape
  class_q <- q + geom_point(mapping = aes(colour = Class)) +
          scale_colour_viridis_d(na.value = "transparent") +
          labs(title = "UMAP coloured by class",
               colour = "Sample class",
               caption = "Replicates are connected by plum-coloured lines \nDoubles connected by black dotted lines") +
          # Connect samples that are replicates of eachother
          geom_line(data = na.omit(umap_data[, c("UMAP1", "UMAP2", "repl.groups")]),
                    mapping = aes(group = repl.groups,
                                  x = UMAP1, y = UMAP2),
                    colour = "plum", show.legend = F) +
          # Connect samples that are doubles of eachother
          geom_line(data = na.omit(umap_data[, c("UMAP1", "UMAP2", "other.groups")]),
                    mapping = aes(group = other.groups,
                                  x = UMAP1, y = UMAP2),
                    colour = "black", lty = "dotted", show.legend = F)
  
  umap_list <- list("hospital"=hosp_q, 
                    "age"=age_q,
                    "bmi"=bmi_q,
                    "entrydate"=entrydate_q,
                    "pdens"=pdens_q,
                    "adens"=adens_q,
                    "postmenop"=postmenop_q,
                    "class"=class_q)
  
  if (printall) {
    for (u in 1:length(umap_list)) {
      print(umap_list[[u]])
    }
  } else {
    return(umap_list)
  }
  
}


# Plotting correlation density to compare replicates, doubles and other random sample pairs, 
# used in normalisation.Rmd
cor_dens <- function(m, s, xrange=NULL, yrange=NULL, plotname="") {
  # Function for plotting density of correlations between replicates, 
  # between doubles and between random other samples
  # m is a data frame containing MFI values and s is a data frame with 
  # the corresponding sample informaion. Both data frames should have the sample names as row names
  # Specify x limits and y limits using 
  
  # Correlation between replicate samples and their originals
  # originals labelled with a 1
  ori_samp <- m[rownames(s[s$replicate == 1 & s$Class != "Other", ]), -1]       
  # replicates have a 2
  repl_samp <- m[rownames(s[s$replicate == 2 & s$Class != "Other", ]), -1]    
  
  repl_cor <- diag(            # Pick out diagonal to get correlation between corresponding samples
    cor(t(ori_samp), t(repl_samp), method="spearman")    
  )
  
  # Correlation between doubles, omit the replicates. Add tube labels to know which are pairs
  doubles <- m[rownames(s[s$Class == "Other" & s$replicate != 2, ]),
               -1]
  doubles$labels <- s[rownames(doubles), "Tube.label"]
  
  # Order to get pairs together and split into different matrices assume an even number of doubles)
  doubles <- doubles[order(doubles$labels), ]
  double_cor <- diag(
    cor(t(doubles[seq(1, nrow(doubles), 2), -which(colnames(doubles) == "labels")]),
        t(doubles[seq(2, nrow(doubles), 2), -which(colnames(doubles) == "labels")]),
        method="spearman")
  )
  
  # Getting correlation between random pairs of samples (non-replicate, non-double)
  rest_samp <- m[rownames(s[s$replicate == 0 & s$Class != "Other", ]), -1]
  rand_order <- sample(nrow(rest_samp), nrow(rest_samp))
  rest_cor <- diag(                # This looks terrible, think of a better way?
    cor(t(rest_samp[1:floor(length(rand_order)/2), ]),
        t(rest_samp[(floor(length(rand_order)/2) + 1):length(rand_order), ]),
        method="spearman")
  )
  
  # Plot densities in different layers
  to_plot <- data.frame("Correlation_Coefficient"=c(rest_cor, double_cor, repl_cor),
                        "Sample"=c(rep("Random", length(rest_cor)),
                                   rep("Double", length(double_cor)),
                                   rep("Replicate", length(repl_cor))
                        ))
  
  pl <- ggplot(data = to_plot) +
    geom_density(aes(Correlation_Coefficient, colour = Sample)) +
    scale_colour_discrete(breaks=c("Random", "Double", "Replicate")) + 
    # Vertical lines for medians, but kind of just in the way when the peaks are thin
    # geom_vline(aes(xintercept = median(rest_cor)), lty="dotted") +
    # geom_vline(aes(xintercept = median(double_cor)), lty = "dotted") +
    # geom_vline(aes(xintercept = median(repl_cor)), lty = "dotted") +
    coord_cartesian(xlim = c(0.85, 1)) + # Just hard-coding limits for now, should think of better way
    labs(title = plotname,
         x = "Correlation coefficient",
         y = "Density",
         caption = paste0("Medians,",
                          "\n Random: ", round(median(rest_cor), 3),
                          "\n Double: ", round(median(double_cor), 3),
                          "\n Replicate: ", round(median(repl_cor), 3))) +
    theme(plot.caption = element_text(size = 12),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14),
          axis.text = element_text(size = 11))
  
  # If x limits or y limits are specified (will reset the previous hard-coded limits for x axis)
  if (is.null(xrange) == F & is.null(yrange) == F) {
    pl <- pl + coord_cartesian(xlim = xrange,
                               ylim = yrange)
  } else if (is.null(yrange) == F) {
    pl <- pl + coord_cartesian(ylim = yrange)
  } else if (is.null(xrange) == F) {
    pl <- pl + coord_cartesian(xlim = xrange)
  }
  
  return(pl)
  
}


plot_pc123 <- function(in_df, colourby, only12=F, 
                       plot_palette="Paired", weighted=F,
                       dimred_arch=NULL, printplots=T) {
  ### Function for plotting PCs 1, 2 and 3 from PCA of archetypal analysis results ###
  # in_df is a data frame containing PCs 1-3 (called "PC1", "PC2" and "PC3") 
  # and the variable to colour by for each row.
  # only12, logical, if TRUE plots only first two PCs. 
  # colourby, character, is the name of the variable to colour by (column in in_df), e.g. archetype.
  # plot_palette, character, is the brewer discrete colour palette to use for colouring.
  # weighted, logical, indicates whether the points are weighted. If TRUE, in_df should
  # contain a column named "weights" that specifies the weight for each row.
  # dimred_arch, data frame, is a data frame containing archetypes in the coordinates of the PCA. 
  # If NULL, will skip
  
  q <- ggplot(data = in_df)
  
  # Scatter with archetypes coloured
  pc_xy <- q + theme_bw() + scale_colour_brewer(palette = plot_palette)
  
  if (weighted) {
    pc12 <- pc_xy + 
      geom_point(mapping = aes_string(x = "PC1", y = "PC2", colour = colourby, alpha = "weights")) + 
      labs(title = "PC1 vs PC2")
  
    pc13 <- pc_xy + 
      geom_point(mapping = aes_string(x = "PC1", y = "PC3", colour = colourby, alpha = "weights")) + 
      labs(title = "PC1 vs PC3")
  
    pc23 <- pc_xy + 
      geom_point(mapping = aes_string(x = "PC2", y = "PC3", colour = colourby, alpha = "weights")) + 
      labs(title = "PC2 vs PC3")
  } else {
    pc12 <- pc_xy + 
      geom_point(mapping = aes_string(x = "PC1", y = "PC2", colour = colourby)) + 
      labs(title = "PC1 vs PC2")
    
    pc13 <- pc_xy +
      geom_point(mapping = aes_string(x = "PC1", y = "PC3", colour = colourby)) + 
      labs(title = "PC1 vs PC3")
    
    pc23 <- pc_xy + 
      geom_point(mapping = aes_string(x = "PC2", y = "PC3", colour = colourby)) + 
      labs(title = "PC2 vs PC3")
  }
  
  if (class(dimred_arch) == "data.frame") {
    pc12 <- pc12 + 
      geom_text(data = dimred_arch,
                 mapping = aes(x = PC1, y = PC2), 
                 label = as.character(1:nrow(dimred_arch)), 
                 size = 4)
    
    pc13 <- pc13 + 
      geom_text(data = dimred_arch,
                 mapping = aes(x = PC1, y = PC3), 
                 label = as.character(1:nrow(dimred_arch)), 
                 size = 4)
    
    pc23 <- pc23 + 
      geom_text(data = dimred_arch,
                 mapping = aes(x = PC2, y = PC3), 
                 label = as.character(1:nrow(dimred_arch)), 
                 size = 4)
  }
  
  
  dsgn <- "A#
           BC"
  
  if (only12) {
    if (printplots) {
      print(pc12)
    } else {
      return(pc12)
    }
  } else {
    pc123 <- wrap_plots(A = pc12, B = pc13, C = pc23, design = dsgn, guides = "collect")
    if (printplots) {
      print(pc123)
    } else {
      return(list("pc123"=pc123,
                  "pc12"=pc12,
                  "pc23"=pc23,
                  "pc13"=pc13,
                  "design"=dsgn))
    }
    
  }
  
}


print_packversions <- function() {
  # Function for printing the versions of the used packages
  
  pack <- .packages()
  pack_ver <- vector("character", length(pack))
  for (p in 1:length(pack)) {
    pack_ver[p] <- paste(pack[p], packageVersion(pack[p]), sep=": ")
  }

  print(c("Package versions:", pack_ver))
}


nice_barplot <- function(plotlist, legendtitle, legendlabels) {
  # Make a nicer looking barplots (for e.g. reports or presentations)
  # from the output list of the plot_results function. 
  # Only for categorical data (i.e. that get barplots)
  # plotlist is the output list of plot_results
  # legendtitle is the title of the legend
  # legendlabes is a vector containing the names of the legend labels
  (plotlist$leftplot +
     labs(fill = legendtitle) + 
     scale_fill_hue(labels = legendlabels)) + 
    (plotlist$rightplot + 
       theme(legend.position = "none")) + 
    plot_layout(guides = "collect")
}
