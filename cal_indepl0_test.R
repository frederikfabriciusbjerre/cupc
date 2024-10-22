library(pcalg)
library(graph)
library(MASS)
library(tictoc)
library(igraph)
library(mice)

source("cuPC.R")
source("cuPCMI.R")


set.seed(1)
p <- 10
prob_dag <- 0.55
n <- 1000
alpha <- 0.01

# Simulate random DAG and data
dag_true <- randomDAG(p, prob = prob_dag)
cpdag_true <- dag2cpdag(dag_true)
print(cpdag_true)
data <- rmvDAG(n, dag_true, errDist = "normal", mix = 0.3)

suffStatMI <- getSuff(list(data, data, data, data, data, data, data, data))
suffStat <- list(C = cor(data), n = nrow(data))
cat("\n")
tic()
cuPC_fit <- cu_pc(suffStat, p = p, alpha = alpha, m.max = 0)
print("the total time consumed by stable.fast is:")
toc()
cat("\n")
print(cuPC_fit)
cat("\n")


tic()
cuPCMI_fit <- cu_pc_MI(suffStatMI, p = p, alpha = alpha, m.max = 0)
print("The total time consumed by cuPC is:")
toc()
cat("\n")
print(cuPCMI_fit)
cat("\n")
