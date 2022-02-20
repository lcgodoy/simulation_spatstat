library(dplyr)

my_power <- function(x, sig_level = .05) {
    mean(x <= sig_level)
}

my_power_ci <- function(x, sig_level = .05) {
    pp <- mean(x <= sig_level)
    sd_pp <- sqrt(pp * (1 - pp) / length(x))
    sprintf("[%.3f; %.3f]",
            pmax(pp - qnorm(1 - (.5 * sig_level),
                            sd = sd_pp), 0),
            pmin(pp + qnorm(1 - (.5 * sig_level),
                            sd = sd_pp), 1))
}

my_files <- list.files("data/results",
                       pattern = "^(indep|att|rep)",
                       full.names = TRUE)

my_results <- lapply(my_files, readRDS)

my_results <- bind_rows(my_results)

my_results |>
    filter(sim_id <= 1000) |>
    tidyr::pivot_longer(cols = 3:4,
                        names_to  = "method",
                        values_to = "pval") |>
    group_by(scenario, method) |>
    summarise(power = my_power(pval),
              ci_pow = my_power_ci(pval))
