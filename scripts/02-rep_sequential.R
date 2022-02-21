source("scripts/utils_sequential.R")

box_sim <- matrix(rep(c(0, 1), 2L),
                  ncol = 2L,
                  byrow = TRUE)

n <- 50L ## controls the number of polygons within each set
size <- .02 ## parameter that control the size of the polygons
shape <- 5L ## shape of the polygons
rel <- "repulsion"
m_hc <- seq(from = 0,
            to   = tpsautils::max_hc(n, box_sim),
            length.out = 8)[6]
## n_repl <- 100L
n_repl <- 1000

output <- vector(mode = "list", length = n_repl)

for(i in seq_along(output)) {
    set.seed(i)
    aux <- tpsautils::sim_data(n_1 = n, n_2 = n,
                               relation = rel,
                               r = size, points_chull = shape,
                               max_hc = m_hc[i],
                               bbox = box_sim)

    aux_lavan <- compute_lavan_pvals(sp_lst = aux, bb = box_sim,
                                     return_ts = FALSE,
                                     path_gcops = "lavancier/bin")
    output[[i]] <- data.frame(sim_id    = i,
                              scenario  = rel,
                              lavancier = aux_lavan,
                              mc_test   = tpsa::gof_mc(aux[[1]], aux[[2]],
                                                       n_sim = 499L,
                                                       unique_bbox = bb,
                                                       alpha = 0.05,
                                                       H = 'L',
                                                       ts = 'SMAD',
                                                       distances = NULL,
                                                       fixed = FALSE,
                                                       method = 'hausdorff')$p_value)
}

saveRDS(output, "data/raw/rep.rds")
