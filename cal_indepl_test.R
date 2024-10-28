library(pcalg)
library(graph)
library(MASS)
library(tictoc)
library(igraph)
library(mice)
#library(tidyverse)
library(micd)

source("cuPC.R")
source("cuPCMI.R")

set.seed(9)
p <- 22
prob_dag <- 0.25
prob_miss <- 0.1
n <- 10000
alpha <- 0.1
max_order <- 10

# Simulate random DAG and data
dag_true <- randomDAG(p, prob = prob_dag)
cpdag_true <- dag2cpdag(dag_true)
print(cpdag_true)

# now scaled
data <- rmvDAG(n, dag_true, errDist = "normal", mix = 0.3) #%>% scale()

# missing at random data 
data_missing <- ampute(data, prop = prob_miss, 
                        mech = "MAR", 
                        bycases = TRUE)$amp

# naive mice imputation
tic()
imputed_data <- mice(data_missing, m = 10, method = 'norm', printFlag = TRUE)
toc()

suffStatMI <- getSuffCU(imputed_data) 
suffStatMICD <- micd::getSuff(imputed_data, test="gaussMItest")

tic()
micd_PC <- pc(suffStatMICD, indepTest = gaussMItest, p = p, alpha = alpha, skel.method = "stable",
 m.max = max_order)
print("The total time consumed by micd_PC is:")
toc()
cat("\n")
cat("micd_PC\n")
print(micd_PC)
cat("\n")

tic()
cuPCMI_fit <- cu_pc_MI(suffStatMI, p = p, alpha = alpha, m.max = max_order)
print("The total time consumed by cuPCMI is:")
toc()
cat("\n")
cat("cuPCMI\n")
print(cuPCMI_fit)
cat("\n")

cat("micdPC ord:", micd_PC@max.ord, "\n")
cat("cuPC ord:", cuPCMI_fit@max.ord, "\n")
shdSkeleton <- function(fit1, fit2){
  graph1 <- fit1 %>% getGraph() %>% ugraph()
  graph2 <- fit2 %>% getGraph() %>% ugraph()
  return (shd(graph1, graph2))
}

# if (require(Rgraphviz)) {
#   ## show estimated CPDAG
#   par(mfrow = c(1, 2))
#   plot(micd_PC, main = "Estimated CPDAG (micd_PC)")
#   plot(cuPCMI_fit, main = "Estimated CPDAG (cuPC)")
# }



source("printfunc.R")
sepset1 <- cuPCMI_fit@sepset
sepset2 <- micd_PC@sepset
findDiffIndexes(sepset1, sepset2)


# Flatten the two sepsets into a dataframe
df_values_indices <- flatten_two_sepsets_with_indices(sepset1, sepset2)

# # Print the dataframe
#print(df_values_indices)
# source("gaussMItestPrint.R")
# gaussMItest(10, 9, c(4,7,8,12), suffStatMICD)

cat("Hamming Distance            =", shdSkeleton(micd_PC, cuPCMI_fit), "\n")
cat("Structural Hamming Distance =", shd(micd_PC, cuPCMI_fit), "\n")

#system("R")
