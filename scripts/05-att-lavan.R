source("scripts/utils_sequential.R")

box_sim <- matrix(rep(c(0, 1), 2L),
                  ncol = 2L,
                  byrow = TRUE)

n <- 50L ## controls the number of polygons within each set
size <- .02 ## parameter that control the size of the polygons
shape <- 5L ## shape of the polygons
rel <- "attraction"
r_at <- seq(from = .8,
            to   = 1e-4,
            length.out = 8)[7]
n_repl <- 100L
## n_repl <- 1000

output <- vector(mode = "numeric", length = n_repl)

for(i in seq_along(output)) {
    set.seed(i)
    aux <- tpsautils::sim_data(n_1 = n, n_2 = n,
                               relation = rel,
                               r = size, points_chull = shape,
                               r_att = r_at,
                               bbox = box_sim)

    output[[i]] <- compute_lavan_pvals(sp_lst = aux, bb = box_sim,
                                       path_gcops = "lavancier/bin")
}

saveRDS(output,
        sprintf("data/raw/lavan_%s.rds", rel))
