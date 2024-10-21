library(pcalg)
library(graph)
library(MASS)
library(tictoc)
library(igraph)

source("cuPC.R")

p <- 
prob_dag <- 0.55
n <- 1000
alpha <- 0.01

# Simulate random DAG and data
dag_true <- randomDAG(p, prob = prob_dag)
cpdag_true <- dag2cpdag(dag_true)
data <- rmvDAG(n, dag_true, errDist = "normal", mix = 0.3)

suffStat <- list(C = cor(data), n = nrow(data))
cat("\n")
tic()
stable_fast_fit <- pc(suffStat, indepTest=gaussCItest, p=p, skel.method="stable.fast", alpha=alpha)
print("the total time consumed by stable.fast is:")
toc()
cat("\n")
print(stable_fast_fit)
cat("\n")


tic()
cuPC_fit <- cu_pc(suffStat, p=p, alpha=alpha)
print("The total time consumed by cuPC is:")
toc()
cat("\n")
print(cuPC_fit)
cat("\n")