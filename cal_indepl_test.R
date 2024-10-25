library(pcalg)
library(graph)
library(MASS)
library(tictoc)
library(igraph)
library(mice)
#library(micd)

source("cuPC.R")
source("cuPCMI.R")
gaussMItest <- function (x, y, S, suffStat) {
  # number of imputations
  M <- length(suffStat) - 1
  # sample size
  n <- suffStat[[M+1]]
  suffStat[[M+1]] <- NULL
  
  z <- sapply(suffStat, function(j) {
    zStatMI(x, y, S, C=j, n=n)
  })
  
  # 1. Average of M imputed data sets
  avgz <- mean(z)

  # 2. Average of completed-data variance
  W <- 1 / (n - 3)
  
  # 3. Between variance
  B <- sum( ( z - avgz )^2 ) / (M-1)
  
  # 4. Total variance
  TV <- W + (1 + 1 / M) * B
  
  # 5. Test statistic
  ts <- avgz / sqrt(TV)
  
  # 6. Degrees of freedom
  df <- (M - 1) * (1 + (W / B) * (M/(M + 1)))^2
  
  # 7. pvalue
  pvalue <- 2 * stats::pt(abs(ts), df = df, lower.tail = FALSE)

  #cat(sprintf("avgz: %f, W: %f, B: %f, TV: %f, ts: %f, df: %f, p_val: %f\n", avgz, W, B, TV, ts, df, pvalue))
  return(pvalue)
}
zStatMI <- function (x, y, S, C, n)
{
  r <- pcalg::pcorOrder(x, y, S, C)
  res <- 0.5 * log_q1pm(r)
  if (is.na(res))
    0
  else res
}

log_q1pm <- function(r) log1p(2 * r / (1 - r))

set.seed(11337)
p <- 100
prob_dag <- 0.55
prob_miss <- 0.25
n <- 3000
alpha <- 0.1
max_order <- 2

# Simulate random DAG and data
dag_true <- randomDAG(p, prob = prob_dag)
cpdag_true <- dag2cpdag(dag_true)
print(cpdag_true)
data <- rmvDAG(n, dag_true, errDist = "normal", mix = 0.3)

# missing at random data 
data_missing <- ampute(data, prop = prob_miss, 
                        mech = "MAR", 
                        bycases = TRUE)$amp

# naive mice imputation
tic()
imputed_data <- mice(data_missing, m = 10, method = 'pmm', printFlag = FALSE)
toc()

suffStatMI <- getSuffCU(imputed_data) 
suffStatMICD <- micd::getSuff(imputed_data, test="gaussMItest")

tic()
micd_PC <- pc(suffStatMICD, indepTest = gaussMItest, p = p, alpha = alpha, skel.method = "stable",
 m.max = Inf)
print("The total time consumed by micd_PC is:")
toc()
cat("\n")
cat("micd_PC\n")
print(micd_PC)
cat("\n")
micd_PC@max.ord

tic()
cuPCMI_fit <- cu_pc_MI(suffStatMI, p = p, alpha = alpha, m.max = max_order)
print("The total time consumed by cuPCMI is:")
toc()
cat("\n")
cat("cuPCMI\n")
print(cuPCMI_fit)
cat("\n")
# shdSkeleton <- function(fit1, fit2){
#   graph1 <- fit1 %>% getGraph() %>% ugraph()
#   graph2 <- fit2 %>% getGraph() %>% ugraph()
#   return (shd(graph1, graph2))
# }
# shdSkeleton(micd_PC, cuPCMI_fit)
# cuPCMI_fit@sepset
# micd_PC@sepset
# # if (require(Rgraphviz)) {
# #   ## show estimated CPDAG
# #   par(mfrow = c(1, 2))
# #   plot(micd_PC, main = "Estimated CPDAG (micd_PC)")
# #   plot(cuPCMI_fit, main = "Estimated CPDAG (cuPC)")
# # }



# #system("R")
# # Function to flatten two nested lists into a dataframe with indices
flatten_two_sepsets_with_indices <- function(sepset1, sepset2) {
  # Initialize empty vectors to store indices and values
  indices_i <- c()
  indices_j <- c()
  values1 <- c()
  values2 <- c()
  
  # Determine the maximum dimensions
  max_i <- max(length(sepset1), length(sepset2))
  
  # Iterate over the outer list indices
  for (i in seq_len(max_i)) {
    inner_list1 <- if (i <= length(sepset1)) sepset1[[i]] else list()
    inner_list2 <- if (i <= length(sepset2)) sepset2[[i]] else list()
    
    max_j <- max(length(inner_list1), length(inner_list2))
    
    # Iterate over the inner list indices
    for (j in seq_len(max_j)) {
      val1 <- if (j <= length(inner_list1)) inner_list1[[j]] else NULL
      val2 <- if (j <= length(inner_list2)) inner_list2[[j]] else NULL
      
      # Process val1
      if (!is.null(val1)) {
        if (is.atomic(val1)) {
          # Handle multiple values in val1
          vals1 <- as.character(val1)
        } else {
          vals1 <- "<complex>"
        }
      } else {
        vals1 <- NA
      }
      
      # Process val2
      if (!is.null(val2)) {
        if (is.atomic(val2)) {
          # Handle multiple values in val2
          vals2 <- as.character(val2)
        } else {
          vals2 <- "<complex>"
        }
      } else {
        vals2 <- NA
      }
      
      # Determine the maximum number of values at this position
      max_k <- max(length(vals1), length(vals2))
      
      for (k in seq_len(max_k)) {
        indices_i <- c(indices_i, i)
        indices_j <- c(indices_j, j)
        values1 <- c(values1, if (k <= length(vals1)) vals1[k] else NA)
        values2 <- c(values2, if (k <= length(vals2)) vals2[k] else NA)
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

# Example usage:
# Assuming 'sepset1' and 'sepset2' are your two nested lists
sepset1 <- cuPCMI_fit@sepset      # Replace with your first sepset
sepset2 <- micd_PC@sepset     # Replace with your second sepset

# Flatten the two sepsets into a dataframe
df_values_indices <- flatten_two_sepsets_with_indices(sepset1, sepset2)

# Print the dataframe
print(df_values_indices)

#sepset1
