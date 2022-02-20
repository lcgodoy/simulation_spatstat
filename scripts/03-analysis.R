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

if(!file.exists("data/analysis_dt.csv")) {
    my_files <- list.files("data/results",
                           pattern = "^(indep|att|rep)",
                           full.names = TRUE)
    
    my_results <- lapply(my_files, readRDS)
    
    my_results <- bind_rows(my_results)

    readr::write_csv(my_results, "data/analysis_dt.csv")
} else
    my_results <- readr::read_csv("data/analysis_dt.csv")

##--- sig_level = .05 ----

my_results |>
    tidyr::pivot_longer(cols = 3:4,
                        names_to  = "method",
                        values_to = "pval") |>
    mutate(scenario = factor(scenario,
                             levels = c("indep", "att",
                                        "rep"),
                             ordered = TRUE)) |>
    group_by(scenario, method) |>
    summarise(power = my_power(pval, sig_level = .05),
              ci_pow = my_power_ci(pval, sig_level = .05))

##--- evaluating at different significance levels ----

my_results <- mutate(my_results,
                     alpha_01 = .01,
                     alpha_05 = .05,
                     alpha_10 = .1)

my_results |>
    tidyr::pivot_longer(cols = 3:4,
                        names_to  = "method",
                        values_to = "pval") |>
    tidyr::pivot_longer(cols = 3:5,
                        names_to  = "alpha",
                        values_to = "sig_level") |>
    select(-alpha) |>
    mutate(scenario = factor(scenario,
                             levels = c("indep", "att",
                                        "rep"),
                             ordered = TRUE)) |>
    group_by(scenario, method, sig_level) |>
    summarise(power = my_power(pval, sig_level = unique(sig_level)),
              ci_pow = my_power_ci(pval, sig_level = unique(sig_level)))


##--- test statistic ----

if(!file.exists("data/test_stat")) {
    test_stat <- list.files("data/results",
                            pattern = "^ts",
                            full.names = TRUE)
    
    test_stat <- sapply(test_stat, readRDS,
                        USE.NAMES = FALSE)
    
    write(test_stat, file = "data/test_stat")
} else
    test_stat <- scan("data/test_stat")

hist(test_stat, freq = FALSE,
     xlim = c(-4, max(c(test_stat, 4))))
curve(dnorm, from = -5, to = 5, col = 2,
      lwd = 2, add = TRUE)

qqnorm(test_stat, pch = 19,
       xlim = c(-4, 4),
       ylim = c(-4, 4))
abline(a = 0, b = 1, lty = 2, col = 2, lwd = 2)

ks.test(test_stat, "pnorm")
