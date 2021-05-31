### Function for visualisation of results in Karma archetypal analysis pipeline ###
### Mainly for visualising archetypes, but can take other groupings ### 
# Leo Dahl
# 2020-03-29

plot_results <- function(input_df, groupby, clin_var, 
                         plot_colour="Paired", stars=T,
                         dens_plot=T, printplots=T, printp=T) {
  # Function for plotting clinical parameter to look at differences
  # input_df is a data frame containing at least the column to group by and the column with the clinical variable to analyse
  # groupby is the name of the column containing the grouping
  # clin_var is the name of the column with the clinical variable of interest
  # plot_colour is the brewer colour set to use for the groups
  # stars, logical, whether to change pairwise p-values to asterisks or not
  # dens_plot, logical, indicates whether to make a density plot for continuous variables or not
  # printplots, logical, whether to print the plots immediately. If FALSE, output is a list of different output plots and p-values
  # printp, logical, print a matrix of p-values if TRUE  
  
  # If clin_var is a continuous variable (numeric class), boxplots are made
  # If clin_var is a discrete variable (character class), bar charts are made
  select <- dplyr::select
  
  if (class(input_df[, groupby]) != "character") {
    print("Make sure the grouping variable is of the character class!")
  }
  if (class(input_df[, clin_var]) != "numeric" &
      class(input_df[, clin_var]) != "integer" &
      class(input_df[, clin_var]) != "character") {
    print("Make sure the clinical variable is of the integer or numeric classes (continuous variable) or of the character class (categorical variable)!")
  }
  
  # First letter to upper case for axes
  clin_ax <- clin_var
  substr(clin_ax, 1, 1) <- toupper(substr(clin_ax, 1, 1))
  groupax <- groupby
  substr(groupax, 1, 1) <- toupper(substr(groupax, 1, 1))
  
  # Get groups
  grp <- sort(unique(input_df[, groupby]))
  n_grp <- length(grp)
  
  # Combinations of groups
  
  pairs = pair_out <- data.frame(combn(grp, 2))
  
  q <- input_df %>%
    select(groupby, clin_var) %>%
    ggplot()
  
  if (class(input_df[, clin_var]) %in% c("numeric", "integer")) {
    ### Continuous variables (should be numeric) ###
    
    # Count number of samples in each group (not NA)
    group_n <- input_df %>%
      filter(is.na(get(clin_var)) == F) %>%
      count(get(groupby)) %>%
      as.data.frame()
    
    colnames(group_n) <- c(groupby, "n")
    group_n_axistext <- paste0(group_n[, groupby], "\nn=", group_n$n)
    
    ### Make density plot & count plot ###
    if (dens_plot) {
      d <- q +
        geom_density(mapping = aes_string(x = clin_var, colour = groupby), 
                     size = 0.7, show.legend = F) + 
        labs(title = paste(clin_ax, "distribution for", groupax)) + 
        scale_colour_brewer(palette = plot_colour)
      
      f <- q +
        geom_freqpoly(mapping = aes_string(x = clin_var, colour = groupby), size = 0.7) + 
        labs(title = paste(clin_ax, "counts for", groupax)) +
        scale_colour_brewer(palette = plot_colour)
    }
    ### Make boxplots ###

    # If more than one group, do Wilcoxon rank-sum test (Mann-Whitney U test)
    if (n_grp > 1) {
      
      pvals <- vector("numeric", ncol(pairs))
      for (i in 1:ncol(pairs)) {
        pvals[i] <- wilcox.test(
          x=input_df[which(input_df[, groupby] == pairs[1, i]), 
                     clin_var],
          y=input_df[which(input_df[, groupby] == pairs[2, i]),
                     clin_var]
        )$p.value
      }

      pval_out <- pvals
      # Pick out comparisons with p-value below 0.05
      pairs <- pairs[which(pvals < 0.05)]
      pvals <- pvals[which(pvals < 0.05)]
      
      if (length(pvals) > 0) {
        # Put p-values in plot if any are below 0.05
        # Swap out p-values for stars, should move somewhere else?
        # "*" if lower than 0.05, "**" if lower than 0.01, "***" if lower than 0.001
        if (stars) {
          pval_num <- pvals
          pvals <- symnum(pvals,
                          cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                          symbols = c("***", "**", "*", "n.s"),
                          corr = F, legend = F) %>%
            as.character()
        }
        
        q <- input_df %>%
          select(groupby, clin_var) %>%
          ggplot(aes_string(x = groupby, y = clin_var, fill = groupby)) +
          geom_boxplot() + 
          scale_fill_brewer(palette = plot_colour) + 
          geom_signif(comparisons = pairs, 
                      step_increase = 0.07,
                      annotations = as.character(pvals),
          						size = 1, 
          						textsize = 8,
          						vjust = 0.5) +
          labs(title = paste(clin_ax, "vs", groupax, "and significant Wilcoxon tests")
               , x = groupax, y = clin_ax, fill = groupax) + 
          scale_x_discrete(limits = grp, 
                           labels = group_n_axistext) 
      } else {
        # Just make the boxplot without p-values if no p-values below 0.05
        q <- input_df %>%
          select(groupby, clin_var) %>%
          ggplot(aes_string(x = groupby, y = clin_var, fill = groupby)) +
          geom_boxplot() +
          scale_fill_brewer(palette = plot_colour) + 
          labs(title = paste(clin_ax, "vs", groupax),
               x = groupax, y = clin_ax, fill = groupax) + 
          scale_x_discrete(limits = grp,
                           labels = group_n_axistext)
      }
    } else {
      # If only two groups, skip Wilcoxon test
      pval_out <- vector("character", nrow(pairs))
      q <- q +
        geom_boxplot(mapping = aes_string(x = groupby, y = clin_var, fill = groupby)) + 
        scale_x_discrete(limits = grp, 
                         labels = group_n_axistext) + 
        labs(x = groupax, y = clin_ax, fill = groupax) +
        scale_colour_brewer(palette = plot_colour)
    }
    
  } else if (class(input_df[, clin_var]) == "character") {
    ### Discrete variables (should be characters) ###
    ###  Make bar plots ###
    
    ## Count number of non-NA samples in each cluster (for plotting) ##
    group_n <- input_df %>%
      filter(is.na(get(clin_var)) == F) %>%
      count(get(groupby)) %>%
      as.data.frame()
    
    colnames(group_n) <- c(groupby, "n")
    group_n_axistext <- paste0(group_n[, groupby], "\nn=", group_n$n)
    
    ## Fisher exact test ##
    elem <- sort(unique(input_df[, clin_var]))
    
    if (length(elem) > 1) {
      # Perform Fisher exact tests if there are multiple classes
      
      # Make contingency table
      cont_tbl <- table(input_df[, groupby], input_df[, clin_var])
      
      # Pair-wise Fisher exact test if more than one group
      if (n_grp > 1) {
        
        pvals <- vector("numeric", ncol(pairs))

        # Pick out two rows at a time from contingency table and perform test
        
        for (pair in 1:ncol(pairs)) {
          tbl <- cont_tbl[pairs[, pair], ]
          pvals[pair] <- fisher.test(tbl)$p.value %>%
            signif(3)
        }

        pval_out <- pvals
        # Pick out comparisons with p-value below 0.05
        pairs <- pairs[which(pvals < 0.05)]
        pvals <- pvals[which(pvals < 0.05)]
        
        if (length(pvals) > 0) {
          # Put p-values in plot if any are below 0.05
          # Swap out p-values for stars, should move somewhere else?
          # "*" if lower than 0.05, "**" if lower than 0.01, "***" if lower than 0.001
          
          # Convert p-values to stars
          if (stars) {
            pval_num <- pvals
            pvals <- symnum(pvals,
                            cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                            symbols = c("***", "**", "*", "n.s"),
                            corr = F, legend = F) %>%
              as.character()
          }

          # Set y positions for p-values
          ymin <- 1.03
          ymax <- ymin + (length(pvals)*0.07)
          ypos <- seq(ymin, ymax, by=0.07)
          
          # Bar plot with significance values, 
          # make proportions before ggplot since does not work very well with geom_signif? 
          # use get(clin_var) and get(groupby) since the functions don't 
          # take strings as arguments
          qright <-
            input_df %>%
            filter(is.na(get(clin_var)) == F) %>%
            group_by(get(groupby)) %>%
            count(get(clin_var)) %>%
            mutate(Freq = n / sum(n)) %>%
            ggplot() +
            aes_string(groupby, "Freq", fill = clin_var) +
            geom_col(position = "fill") +
            geom_signif(
              comparisons = pairs,
              annotations = pvals,
              y_position = ypos,
              size = 1, 
              textsize = 8,
              vjust = 0.5) +
            labs(x = groupax, y = "Proportion") +
            scale_x_discrete(limits = grp,
                             labels = group_n_axistext)
          
          # Names are messed up!
          # instead of groupby and clin_var you get "get(groupby)" 
          # and "get(clin_var)" as list element names, so change that 
          names(qright$data)[1:2] <- c(groupby, clin_var)
          
        } else {
          qright <- q +
            geom_bar(data = input_df %>%
                       filter(is.na(get(clin_var)) == F),
                     mapping = aes_string(x = groupby, fill = clin_var), 
                     position = "fill") +
            labs(x = groupax, y = "Proportion") + 
            scale_x_discrete(limits = grp,
                             labels = group_n_axistext)
        }
      } else {
        
        # If no more than two groups, plot normally
        pval_out <- vector("character", nrow(pairs))
        
        qright <- q +
          geom_bar(data = input_df %>%
                     filter(is.na(get(clin_var)) == F),
                   mapping = aes_string(x = groupby, fill = clin_var), 
                   position = "fill") +
          labs(x = groupax, y = "Proportion") + 
          scale_x_discrete(limits = grp,
                           labels = group_n_axistext)
      }
      
    } else {
      # If only one class, plot without performing any tests
      
      qright <- q +
        geom_bar(data = input_df %>%
                   filter(is.na(get(clin_var)) == F),
                 mapping = aes_string(x = groupby, fill = clin_var), 
                 position = "fill") +
        labs(x = groupax, y = "Proportion") + 
        scale_x_discrete(limits = grp,
                         labels = group_n_axistext)
    }
    
    qleft <- q + 
      geom_bar(mapping = aes_string(x = groupby, fill = clin_var), 
               position = "dodge") +
      labs(x = groupax, y = "Count", caption = " ")
    
    q <- qleft + qright + 
      plot_annotation(title = paste(clin_ax, "counts and proportions in", groupax)) +
      plot_layout(guides = "collect")
  }
  
  
  # In case pval_out was never made due to conditions not being met
  if (!exists("pval_out")) {
    pval_out <- NULL
  } else {
    if (length(pval_out) == 1) {
      pmat <- pval_out
    } else {
      # P-value printing
      names(pair_out) <- NULL
      pmat <- data.frame(t(pair_out), pval_out)
      colnames(pmat) <- c(groupax, 
                          paste0(groupax, "2"),
                          "P_value")
      pmat$P_value <- signif(as.numeric(pmat$P_value), 3)
      pmatformula <- as.formula(paste0(groupax, " ~ ", groupax, "2"))
      pmat <- dcast(pmat, pmatformula, value.var = "P_value")
      pmat[is.na(pmat)] <- ""
      rownames(pmat) <- pmat[, 1]
    }
    
    if(printp) {
      print(paste0("P-values for ", clin_ax, " comparisons between ", groupax, " groups:"))
      ifelse(class(pmat) == "data.frame", print(pmat[, -1]), print(pmat))
    }
  }
  
  # Giving the output, either print the plots or give them as a list
  if (class(input_df[, clin_var]) == "character") { 
    # Categorical variable
    if (length(elem) > 1) {
      if (printplots) {
        print(q)
      } else {
        plot_out <- list("plot"=q,
                         "leftplot"=qleft,
                         "rightplot"=qright,
                         "pairs"=pair_out,
                         "p.value"=pval_out,
                         "pmatrix"=pmat)
        return(plot_out)
        }
    } else {
        if (printplots) {
          plot(q)
        }
      }
  } else { 
    # Continuous variable
    if (printplots) {
      if (dens_plot) {
        print(f / d + plot_layout(guides = "collect"))
      }
      print(q)
    } else {
      if (dens_plot) {
        plot_out <- list("boxplot"=q,
                         "density"=d,
                         "frequency"=f,
                         "pairs"=pair_out,
                         "p.value"=pval_out,
                         "pmatrix"=pmat)
      } else {
        plot_out <- q
      }
      return(plot_out)
    }
  }
}
