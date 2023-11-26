library(smoof)
library(hmsr)
library(GA)
library(readr)

sprouting_default_euclidean_distances <- function(sigma) {
  sprouting_condition_distance_ratio <- 0.6
  lapply(sigma, function(x) {
    sum(x * sprouting_condition_distance_ratio)
  })
}

default_sigma <- function(lower, upper, tree_height) {
  sigma_ratio <- 0.04
  sigma_exponent <- 0.5
  domain_length <- upper - lower
  sigma <- list()
  for (height in 1:tree_height) {
    sigma <- c(sigma, list(domain_length * sigma_ratio))
    sigma_ratio <- sigma_ratio * sigma_exponent
  }
  sigma
}

hms_vs_ga <-
  function(evaluations,
           fitness_function) {
    seeds <- 1141:1171

    params <-
      attr(fitness_function, "par.set")$pars$x
    lower <- params[["lower"]]
    upper <- params[["upper"]]
    dimensions <- length(lower)
    sigma <- default_sigma(lower, upper, 3)
    sprouting_distances <-
      sprouting_default_euclidean_distances(sigma)

    fitness <- function(x) {
      -1 * fitness_function(x)
    }

    run_hms <- function(seed) {
      ga_config <- list(
        list(pmutation = 0.4),
        list(
          pmutation = 0.2,
          mutation = rtnorm_mutation(lower, upper, sigma[[2]])
        ),
        list(
          pmutation = 0.2,
          mutation = rtnorm_mutation(lower, upper, sigma[[3]])
        )
      )

      set.seed(seed)

      result <- hmsr::hms(
        fitness = fitness,
        tree_height = 3,
        lower = lower,
        upper = upper,
        run_metaepoch = ga_metaepoch(ga_config),
        population_sizes = c(50, 25, 10),
        sigma = sigma,
        gsc = gsc_max_fitness_evaluations(evaluations),
        sc = sc_max_metric(euclidean_distance, sprouting_distances),
        lsc = lsc_metaepochs_without_improvement(7),
        monitor_level = "none"
      )
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

    hms <- mapply(function(seed) {
      run_hms(seed)
    }, seeds)
    ga <- mapply(function(seed) {
      run_ga(seed)
    }, seeds)

    list(hms = hms, ga = ga)
  }


get_hms_vs_ga_res <- function(budgets, fitness_functions) {
  all_results <- list()
  for (fitness_function in fitness_functions) {
    results <- list()
    for (budget in budgets) {
      result <- hms_vs_ga(budget, fitness_function)
      results[[paste("ga_", budget, sep = "")]] <-
        result$ga
      results[[paste("hms_", budget, sep = "")]] <-
        result$hms
    }
    all_results[[attr(fitness_function, "name")]] <-
      results
  }
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
  function(results, title, dir_path) {
    path <- paste(dir_path, title, ".jpg", sep = "")
    jpeg(path, width = 1600, height = 1200)
    cols <- tail(rainbow(3, s = 0.1), n = 2)
    at <- get_at(budgets)
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
  function(all_results) {
    mean_times <- c()
    mean_values <- c()
    q_25_values <- c()
    q_50_values <- c()
    q_75_values <- c()
    min_values <- c()
    max_values <- c()
    function_name_values <- c()
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
      max = max_values
    )
  }

run_experiment <-
  function(experiment_name,
           fitness_functions,
           budgets) {
    experiment_dir_path <- paste("./", experiment_name, "/", sep = "")
    dir.create(experiment_dir_path)
    all_results <- get_hms_vs_ga_res(budgets, fitness_functions)
    for (fitness_function_name in names(all_results)) {
      fitness_function_results <- all_results[[fitness_function_name]]
      plot_hms_vs_ga(
        fitness_function_results,
        fitness_function_name,
        experiment_dir_path
      )
    }
    data_frame <-
      generate_data_frame(all_results)
    write.csv(
      data_frame,
      paste(experiment_dir_path, "results.csv", sep = "")
    )
  }

# Run experiment:
BUDGETS <- c(1000, 5000, 10000, 15000)
FITNESS_FUNCTIONS <- list(
  smoof::makeAckleyFunction(10),
  smoof::makeEggholderFunction(),
  smoof::makeGriewankFunction(10),
  smoof::makeRastriginFunction(5)
)

run_experiment("main", FITNESS_FUNCTIONS, BUDGETS)

results <- read_csv("main/results.csv")
View(results)