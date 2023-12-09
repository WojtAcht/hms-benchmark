library(smoof)
library(GA)
library(hmsr)
library(readr)
source("benchmark.R")

N <- 10

fitness_functions <- lapply(1:20, function(id) smoof::makeBBOBFunction(N, id, 1L))

budgets <- list(10000)

run_experiment("bbob", fitness_functions, budgets)

bbob_results <- read_csv("bbob/results.csv")
View(bbob_results)
