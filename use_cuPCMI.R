# ================================= A test case =================================
library(pcalg)
library(graph)
library(MASS)
library(tictoc)
library(igraph)
library(mice)
library(micd)

source("cuPCMI.R")

# read data
dataset_path <- file.path("data/dataset.csv", fsep = .Platform$file.sep)
data <- read.table(dataset_path, sep = ",")

prob_miss <- 0.1
alpha <- 0.95
max_order <- 15

# missing at random data 
data_missing <- ampute(data, prop = prob_miss, 
                        mech = "MAR", 
                        bycases = TRUE)$amp

# naive mice imputation
tic()
imputed_data <- mice(data_missing, m = 10, method = "norm", printFlag = FALSE, remove.collinear = TRUE)
toc()

suffStatMI <- getSuffCU(imputed_data) 

# suffStatMICD <- micd::getSuff(imputed_data, test="gaussMItest")
# tic()
# micd_PC <- pc(suffStatMICD, indepTest = gaussMItest, p = p, alpha = alpha, skel.method = "stable.fast", m.max = max_order)
# print("The total time consumed by micd_PC is:")
# toc()
# cat("\n")
# cat("micd_PC\n")
# print(micd_PC)
# cat("\n")

tic()
cuPCMI_fit <- cu_pc_MI(suffStatMI, p = p, alpha = alpha, m.max = max_order)
print("The total time consumed by cuPCMI is:")
toc()
cat("\n")
cat("cuPCMI\n")
print(cuPCMI_fit)
cat("\n")

# cat("micdPC ord:", micd_PC@max.ord, "\n")
cat("cuPC ord:", cuPCMI_fit@max.ord, "\n")