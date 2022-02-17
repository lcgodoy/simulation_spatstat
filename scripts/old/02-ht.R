sinput <- commandArgs(trailingOnly = TRUE)
sinput <- as.integer(sinput[1])

box_sim <- matrix(rep(c(0, 1), 2),
                  ncol = 2, byrow = T)

source("scripts/utils.R")

##--- indep

file_indep <- sprintf("data/raw/indep_%s.rds",
                      formatC(sinput, width = 4,
                              flag = "0"))

poly_indep <- readRDS(file_indep)

indep_res <- compute_pvals(poly_indep, bbox_sim,
                           sinput, "indep")

saveRDS(x = indep_res,
        file = sprintf("data/results/indep_%s.rds",
                       formatC(sinput, width = 4,
                               flag = "0")))

file.remove(file_indep)

rm(poly_indep); gc()

##--- repulsion ----

file_rep <- sprintf("data/raw/rep_%s.rds",
                    formatC(sinput, width = 4,
                            flag = "0"))

poly_rep <- readRDS(file_rep)

indep_rep <- compute_pvals(poly_rep, bbox_sim,
                           sinput, "rep")

saveRDS(x = indep_rep,
        file = sprintf("data/results/rep_%s.rds",
                       formatC(sinput, width = 4,
                               flag = "0")))

file.remove(file_rep)

rm(poly_rep); gc()

##--- attraction ----

file_att <- sprintf("data/raw/att_%s.rds",
                    formatC(sinput, width = 4,
                            flag = "0"))

poly_att <- readRDS(file_att)

indep_att <- compute_pvals(poly_att, bbox_sim,
                           sinput, "att")

saveRDS(x = indep_att,
        file = sprintf("data/results/att_%s.rds",
                       formatC(sinput, width = 4,
                               flag = "0")))

file.remove(file_att)

rm(poly_att); gc()
