### Antibody QC, criterion 1, high MFI in empty wells
# Leo Dahl
# 2020-02-17

empty_idx <- which(sinfo$Class == "Empty")

# Get mean over all antibodies for each row with an empty well
empty_means <- lapply(mfi_list, 
                      function(assay) rowMeans(assay[empty_idx, ]))

# Get 3*(standard deviation) too
empty_sd <- lapply(mfi_list,
                   function(assay) apply(assay[empty_idx, ], 1, 
                                         function(w) 3*sd(w)))

# Calculating the cutoff MFI for each empty well for each assay
crit1_cutoffs <- vector("list", length(mfi_list))
for (i in 1:length(crit1_cutoffs)) {
  crit1_cutoffs[[i]] <- Map("+", empty_means[[i]], empty_sd[[i]])
}


# Check which antibodies have too high MFI
# Loop through each array with i
# Loop through each row with empty well with j
# Check each column (antibody) with the apply function
crit1_high_mfi <- vector("list", length(mfi_list))

for (i in 1:length(crit1_high_mfi)) {
  v <- vector("list", length(empty_idx))
  for (j in 1:length(empty_idx)) {
    v[[j]] <- which(apply(mfi_list[[i]], 2, 
                          function(column) column[empty_idx[j]] > crit1_cutoffs[[i]][j]))
  }
  crit1_high_mfi[[i]] <- v
}

