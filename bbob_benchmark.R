library(smoof)
library(GA)
library(hmsr)
library(readr)
source("benchmark.R")

N <- 10

fitness_functions <- lapply(1:20, function(id) smoof::makeBBOBFunction(N, id, 1L))

budgets <- list(10000)

run_experiment("bbob_cma_es_two_levels", fitness_functions, budgets, "cmaes_two_levels")

bbob_cmaes_results <- read_csv("bbob_cma_es_two_levels/results.csv")
high_mutation_results <- read_csv("bbob_high_mutation_60/results.csv")
View(bbob_results)
