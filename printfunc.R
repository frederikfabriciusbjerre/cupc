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


findDiffIndexes <- function(lst1, lst2){
  # make lists symmetric (like symmetric matrix)
  lst1 <- makeSymmetric(lst1)
  lst2 <- makeSymmetric(lst2)
  
  # init index list 
  differences <- list()
  n <- length(lst1)
  
  # check if they are comparable
  if (n != length(lst2)){
    cat("lists are not of same length")
    return(FALSE)
  }
  # compare only the upper triangular elements, including the diagonal
  for (i in seq_len(n)) {
    for (j in i:n) {
      if (!identical(lst1[[i]][[j]], lst2[[i]][[j]])) {
        differences <- append(differences, list(c(i, j)))
      }
    }
  }
  return(differences)
}

makeSymmetric <- function(lst){
  # assumes nxn list of lists
  for (i in seq_along(lst)){
    for (j in seq_along(lst)){
      elm.ij <- lst[[i]][[j]]
      elm.ji <- lst[[j]][[i]]
      if (identical(elm.ij, elm.ji) && !is.null(elm.ij) && !is.null(elm.ji)) {
        elm <- elm.ij
      }
      else{
        elm <- handle_union(elm.ij, elm.ji)
      }
      # set both to elm
      lst[[i]][[j]] <- elm
    }
  }
  return (lst)
}

handle_union <- function(x, y) {
  # return NULL if both is NULL
  if (is.null(x) && is.null(y)) {
    return(NA)
  }
  # return y if x is NULL
  if (is.null(x)) {
    return(y)
  }
  # return x if y is NULL
  if (is.null(y)) {
    return(x)
  }
  # return integer(0) if either or both are integer(0)
  if (identical(x, integer(0)) || identical(y, integer(0))) {
    return(integer(0))
  }
    # return integer(0) if either or both are integer(0)
  if (identical(x, 0) || identical(y, 0)) {
    return(0)
  }
  # else combine them 
  combined <- union(x, y)
  return(combined)
}
