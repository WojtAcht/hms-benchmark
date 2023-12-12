library(smoof)
library(hmsr)
library(GA)
library(readr)
library(foreach)
library(doParallel)
source("hms_config.R")

registerDoParallel(cores = 6)

hms_vs_ga <-
  function(evaluations,
           fitness_function,
           config = "sea_two_levels") {
    seeds <- 1141:1171

    params <-
      attr(fitness_function, "par.set")$pars$x
    lower <- params[["lower"]]
    upper <- params[["upper"]]

    fitness <- function(x) {
      -1 * fitness_function(x)
    }

    run_hms <- function(seed) {
      set.seed(seed)
      get_params_fn <- hms_configs[[config]]
      params <- get_params_fn(lower, upper)
      params[["fitness"]] <- fitness
      params[["gsc"]] <- gsc_max_fitness_evaluations(evaluations)
      result <- do.call(hmsr::hms, params)
      fitness_value <- -1 * result@best_fitness
      list(
        "fitness" = fitness_value,
        "time" = result@total_time_in_seconds
      )
    }

    evals <- 0
    ga_f <- function(x) {
      evals <<- evals + 1

      if (evals > evaluations) {
        -Inf
      } else {
        fitness(x)
      }
    }

    run_ga <- function(seed) {
      set.seed(seed)
      evals <<- 0
      start <- Sys.time()
      result <-
        GA::ga(
          type = "real-valued",
          fitness = ga_f,
          lower = lower,
          upper = upper,
          monitor = FALSE,
          maxiter = (evaluations / 150) * 5,
          popSize = 150
        )
      end <- Sys.time()
      fitness_value <- -1 * result@fitnessValue
      list("fitness" = fitness_value, "time" = (end - start))
    }

    NULL_RESULT <- list("fitness" = NA, "time" = NA)

    hms <- mapply(function(seed) {
      tryCatch(run_hms(seed), error = function(e) {
        message(e)
        NULL_RESULT
      })
    }, seeds)
    ga <- mapply(function(seed) {
      tryCatch(run_ga(seed), error = function(e) {
        message(e)
        NULL_RESULT
      })
    }, seeds)

    hms <- hms[, !is.na(hms["fitness", ])]
    ga <- ga[, !is.na(ga["fitness", ])]

    list(hms = hms, ga = ga)
  }


get_hms_vs_ga_res <- function(budgets, fitness_functions, config) {
  all_results <-
    foreach(fitness_function = fitness_functions) %dopar% {
      results <- list()
      for (budget in budgets) {
        result <- hms_vs_ga(budget, fitness_function, config)
        results[[paste("ga_", budget, sep = "")]] <-
          result$ga
        results[[paste("hms_", budget, sep = "")]] <-
          result$hms
        results[["global.opt.value"]] <-
          attr(fitness_function, "global.opt.value")
      }
      results
    }
  names(all_results) <-
    lapply(fitness_functions, function(fitness_function) {
      attr(fitness_function, "name")
    })
  all_results
}

get_at <- function(budgets) {
  at <- c()
  index <- 1
  for (budget_index in seq_along(budgets)) {
    at <- c(at, c(index, index + 1))
    index <- index + 4
  }
  at
}

get_names <- function(budgets) {
  names <- c()
  for (budget in budgets) {
    if (budget >= 1000) {
      names <- c(names, rep(paste(budget %/% 1000, "k", sep = ""), 2))
    } else {
      names <- c(names, rep(as.character(budget), 2))
    }
  }
  names
}

plot_hms_vs_ga <-
  function(results, title, dir_path, budgets) {
    path <- paste(dir_path, title, ".jpg", sep = "")
    jpeg(path, width = 1600, height = 1200)
    cols <- tail(rainbow(3, s = 0.1), n = 2)
    at <- get_at(budgets)
    results$global.opt.value <- NULL
    fitness_values <-
      lapply(results, function(result) {
        -1 * unlist(result["fitness", ])
      })
    boxplot(
      fitness_values,
      at = at,
      col = rep(cols, times = 8),
      ylab = "f(x)",
      xlab = "Fitness evaluations",
      main = title,
      names = get_names(budgets)
    )
    points(
      x = at,
      y = -1 * mapply(mean, results),
      col = rep(tail(rainbow(3), n = 2), times = 6),
      pch = 16
    )
    legend(
      "topleft",
      fill = tail(rainbow(3), n = 2),
      legend = c("GA", "HMS"),
      horiz = T
    )
    dev.off()
  }

generate_data_frame <-
  function(all_results, budgets) {
    mean_times <- c()
    mean_values <- c()
    q_25_values <- c()
    q_50_values <- c()
    q_75_values <- c()
    min_values <- c()
    max_values <- c()
    function_name_values <- c()
    global_opt_values <- c()
    budget_values <- c()
    for (fitness_function_name in names(all_results)) {
      fitness_function_values <- all_results[[fitness_function_name]]
      for (budget in budgets) {
        budget_name <- paste("hms_", budget, sep = "")
        values <-
          fitness_function_values[[budget_name]]
        fitness_values <- unlist(values["fitness", ])
        mean_values <-
          c(mean_values, mean(fitness_values))
        q_25_values <-
          c(q_25_values, quantile(fitness_values, 0.25))
        q_50_values <-
          c(q_50_values, quantile(fitness_values, 0.5))
        q_75_values <-
          c(q_75_values, quantile(fitness_values, 0.75))
        min_values <-
          c(min_values, min(fitness_values))
        max_values <-
          c(max_values, max(fitness_values))
        time_values <- unlist(values["time", ])
        mean_times <- c(mean_times, mean(time_values))
        function_name_values <-
          c(function_name_values, fitness_function_name)
        budget_values <-
          c(budget_values, budget)
        global_opt_values <-
          c(global_opt_values, fitness_function_values[["global.opt.value"]])
      }
    }
    data.frame(
      function_name = function_name_values,
      budget = budget_values,
      mean_time = mean_times,
      min = min_values,
      q_25 = q_25_values,
      q_50 = q_50_values,
      mean = mean_values,
      q_75 = q_75_values,
      max = max_values,
      global_opt_values = global_opt_values
    )
  }

run_experiment <-
  function(experiment_name,
           fitness_functions,
           budgets,
           config = "sea_two_levels") {
    experiment_dir_path <- paste("./", experiment_name, "/", sep = "")
    dir.create(experiment_dir_path)
    all_results <- get_hms_vs_ga_res(budgets, fitness_functions, config)
    data_frame <-
      generate_data_frame(all_results, budgets)
    write.csv(
      data_frame,
      paste(experiment_dir_path, "results.csv", sep = "")
    )
    for (fitness_function_name in names(all_results)) {
      fitness_function_results <- all_results[[fitness_function_name]]
      plot_hms_vs_ga(
        fitness_function_results,
        fitness_function_name,
        experiment_dir_path,
        budgets
      )
    }
  }

# Run experiment:
BUDGETS <- c(1000, 5000, 10000, 15000)
FITNESS_FUNCTIONS <- list(
  smoof::makeAckleyFunction(10),
  smoof::makeEggholderFunction(),
  smoof::makeGriewankFunction(10),
  smoof::makeRastriginFunction(5)
)

# run_experiment("custom-cache", FITNESS_FUNCTIONS, BUDGETS)
#
# results <- read_csv("main/results.csv")
# View(results)
