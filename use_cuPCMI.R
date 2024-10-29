# ================================= A test case =================================
library(pcalg)
library(graph)
library(MASS)
library(tictoc)
library(igraph)
library(mice)
library(micd)
library(miceadds)

source("cuPCMI.R")

# read data as imputed_data
dataset_path <- file.path("dataset_imputed/dataset_imputed.Rdata", fsep = .Platform$file.sep)
load.Rdata(dataset_path, "imputed_data")

# make suffStat
suffStatMI <- getSuffCU(imputed_data) 

# input params to pc
p <- imputed_data[[1]] %>% length()
alpha <- 0.875
max_order <- 20

cat("Fitting with alpha =", alpha, "\n")
tic()
cuPCMI_fit <- cu_pc_MI(suffStatMI, p = p, alpha = alpha, m.max = max_order)
cat("\n")
cat("cuPCMI\n")
print(cuPCMI_fit)
cat("\n")
cat("cuPC ord =", cuPCMI_fit@max.ord, "\n")
cat("alpha    =", alpha, "\n\n")
cat("Time consumed:\n")
toc()
cat("\n")

# suffStatMICD <- micd::getSuff(imputed_data, test="gaussMItest")
# tic()
# micd_PC <- pc(suffStatMICD, indepTest = gaussMItest, p = p, alpha = alpha, skel.method = "stable.fast", m.max = max_order)
# print("The total time consumed by micd_PC is:")
# toc()
# cat("\n")
# cat("micd_PC\n")
# print(micd_PC)
# cat("\n")
# cat("cuPC ord:", micd_PC@max.ord, "\n")