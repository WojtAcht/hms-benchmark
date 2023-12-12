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

get_hms_config_2_levels_sea <- function(lower, upper) {
  dimensions <- length(lower)
  sigma <- default_sigma(lower, upper, 2)
  sprouting_distances <-
    sprouting_default_euclidean_distances(sigma)
  ga_config <- list(
    list(pmutation = 0.6, mutation = rtnorm_mutation(lower, upper, sigma[[1]])),
    list(
      pmutation = 0.2,
      mutation = rtnorm_mutation(lower, upper, sigma[[2]])
    )
  )
  list(
    tree_height = 2,
    lower = lower,
    upper = upper,
    run_metaepoch = ga_metaepoch(ga_config),
    population_sizes = c(50, 25),
    sigma = sigma,
    sc = sc_max_metric(euclidean_distance, sprouting_distances),
    lsc = lsc_metaepochs_without_improvement(10),
    monitor_level = "none"
  )
}

get_hms_config_3_levels_sea <- function(lower, upper) {
  dimensions <- length(lower)
  tree_height <- 3
  sigma <- default_sigma(lower, upper, tree_height)
  sprouting_distances <-
    sprouting_default_euclidean_distances(sigma)
  ga_config <- list(
    list(pmutation = 0.6, mutation = rtnorm_mutation(lower, upper, sigma[[1]])),
    list(
      pmutation = 0.2,
      mutation = rtnorm_mutation(lower, upper, sigma[[2]])
    ),
    list(
      pmutation = 0.2,
      mutation = rtnorm_mutation(lower, upper, sigma[[3]])
    )
  )
  list(
    tree_height = tree_height,
    lower = lower,
    upper = upper,
    run_metaepoch = ga_metaepoch(ga_config),
    population_sizes = c(50, 25, 10),
    sigma = sigma,
    sc = sc_max_metric(euclidean_distance, sprouting_distances),
    lsc = lsc_metaepochs_without_improvement(7),
    monitor_level = "none"
  )
}

get_hms_config_2_levels_cmaes <- function(lower, upper) {
  dimensions <- length(lower)
  sigma <- default_sigma(lower, upper, 2)
  sprouting_distances <-
    sprouting_default_euclidean_distances(sigma)
  ga_config <- list(
    list(pmutation = 0.6, mutation = rtnorm_mutation(lower, upper, sigma[[1]])),
    list()
  )
  cmaes_config <- list(list(), list())
  run_cmaes_metaepoch <- cma_es_metaepoch(cmaes_config)
  run_ga_metaepoch <- ga_metaepoch(ga_config)
  run_metaepoch <- function(fitness,
                            suggestions,
                            lower,
                            upper,
                            tree_level,
                            minimize) {
    if (tree_level == 1) {
      return(run_ga_metaepoch(
        fitness,
        suggestions,
        lower,
        upper,
        tree_level,
        minimize
      ))
    } else {
      return(run_cmaes_metaepoch(
        fitness,
        suggestions,
        lower,
        upper,
        tree_level,
        minimize
      ))
    }
  }
  list(
    tree_height = 2,
    lower = lower,
    upper = upper,
    run_metaepoch = run_metaepoch,
    population_sizes = c(50, 25),
    sigma = sigma,
    sc = sc_max_metric(euclidean_distance, sprouting_distances),
    lsc = lsc_metaepochs_without_improvement(7),
    monitor_level = "none"
  )
}

hms_configs <- list(
  "sea_two_levels" = get_hms_config_2_levels_sea,
  "sea_three_levels" = get_hms_config_3_levels_sea,
  "cmaes_two_levels" = get_hms_config_2_levels_cmaes
)
