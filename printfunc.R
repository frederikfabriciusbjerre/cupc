flatten_two_sepsets_with_indices <- function(sepset1, sepset2, exclude_both_na = TRUE) {
  # Initialize empty vectors to store indices and values
  indices_i <- c()
  indices_j <- c()
  values1 <- c()
  values2 <- c()
  
  # Determine the maximum dimensions
  max_i <- max(length(sepset1), length(sepset2))
  
  # Iterate over the outer list indices
  for (i in seq_len(max_i)) {
    inner_list1_i <- if (i <= length(sepset1)) sepset1[[i]] else list()
    inner_list2_i <- if (i <= length(sepset2)) sepset2[[i]] else list()
    
    max_j_i <- max(length(inner_list1_i), length(inner_list2_i))
    
    # Iterate over the inner list indices
    for (j in seq_len(max_j_i)) {
      # Process (i, j)
      val1_ij <- if (j <= length(inner_list1_i)) inner_list1_i[[j]] else NULL
      val2_ij <- if (j <= length(inner_list2_i)) inner_list2_i[[j]] else NULL
      
      # Process val1_ij
      vals1_ij <- if (!is.null(val1_ij)) {
        if (is.atomic(val1_ij)) as.character(val1_ij) else "<complex>"
      } else {
        NA
      }
      
      # Process val2_ij
      vals2_ij <- if (!is.null(val2_ij)) {
        if (is.atomic(val2_ij)) as.character(val2_ij) else "<complex>"
      } else {
        NA
      }
      
      # Determine the maximum number of values at this position
      max_k_ij <- max(length(vals1_ij), length(vals2_ij))
      
      for (k in seq_len(max_k_ij)) {
        value1 <- if (k <= length(vals1_ij)) vals1_ij[k] else NA
        value2 <- if (k <= length(vals2_ij)) vals2_ij[k] else NA
        
        # If exclude_both_na is TRUE and both values are NA, skip this entry
        if (!(exclude_both_na && is.na(value1) && is.na(value2))) {
          indices_i <- c(indices_i, i)
          indices_j <- c(indices_j, j)
          values1 <- c(values1, value1)
          values2 <- c(values2, value2)
        }
      }
      
      # Process (j, i)
      inner_list1_j <- if (j <= length(sepset1)) sepset1[[j]] else list()
      inner_list2_j <- if (j <= length(sepset2)) sepset2[[j]] else list()
      
      val1_ji <- if (i <= length(inner_list1_j)) inner_list1_j[[i]] else NULL
      val2_ji <- if (i <= length(inner_list2_j)) inner_list2_j[[i]] else NULL
      
      # Process val1_ji
      vals1_ji <- if (!is.null(val1_ji)) {
        if (is.atomic(val1_ji)) as.character(val1_ji) else "<complex>"
      } else {
        NA
      }
      
      # Process val2_ji
      vals2_ji <- if (!is.null(val2_ji)) {
        if (is.atomic(val2_ji)) as.character(val2_ji) else "<complex>"
      } else {
        NA
      }
      
      # Determine the maximum number of values at this position
      max_k_ji <- max(length(vals1_ji), length(vals2_ji))
      
      for (k in seq_len(max_k_ji)) {
        value1 <- if (k <= length(vals1_ji)) vals1_ji[k] else NA
        value2 <- if (k <= length(vals2_ji)) vals2_ji[k] else NA
        
        # If exclude_both_na is TRUE and both values are NA, skip this entry
        if (!(exclude_both_na && is.na(value1) && is.na(value2))) {
          indices_i <- c(indices_i, j)
          indices_j <- c(indices_j, i)
          values1 <- c(values1, value1)
          values2 <- c(values2, value2)
        }
      }
    }
  }
  
  # Create a dataframe with indices and values from both sepsets
  df <- data.frame(
    i = indices_i,
    j = indices_j,
    values_sepset1 = values1,
    values_sepset2 = values2,
    stringsAsFactors = FALSE
  )
  
  return(df)
}
