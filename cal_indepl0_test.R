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
  S = NULL
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

set.seed(1)
p <- 20
prob_dag <- 0.55
prob_miss <- 0.25
n <- 10000
alpha <- 0.01

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
imputed_data <- mice(data_missing, m = 10, method='pmm', printFlag = FALSE)
imputed_data
suffStatMI <- getSuffCU(imputed_data) #hardcoded n=1000
suffStatMI
suffStatMICD <- micd::getSuff(imputed_data, test="gaussMItest")
suffStatMICD
suffStat <- list(C = cor(data), n = nrow(data))
# cat("\n")
# tic()
# cuPC_fit <- cu_pc(suffStat, p = p, alpha = alpha, m.max = 0)
# print("the total time consumed by cuPC is:")
# toc()
# cat("\n")
# cat("cuPC\n")
# print(cuPC_fit)
# cat("\n")

tic()
micd_PC <- pc(suffStatMICD, indepTest = gaussMItest, p = p, alpha = alpha, m.max = 0)
print("the total time consumed by micd_PC is:")
toc()
cat("\n")
cat("micd_PC\n")
print(micd_PC)
cat("\n")

tic()
cuPCMI_fit <- cu_pc_MI(suffStatMI, p = p, alpha = alpha, m.max = 0)
print("The total time consumed by cuPCMI is:")
toc()
cat("\n")
cat("cuPCMI\n")
print(cuPCMI_fit)
cat("\n")
system("R")
