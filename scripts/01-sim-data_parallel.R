##--- simulating independent data ----

sinput <- commandArgs(trailingOnly = TRUE)
sinput <- as.integer(sinput[1])

box_sim <- matrix(rep(c(0, 1), 2),
                  ncol = 2, byrow = T)

n <- 50 ## controls the number of polygons within each set
size <- .02 ## parameter that control the size of the polygons
shape <- 5 ## shape of the polygons
rel <- "independence"
## n_repl <- 1000

set.seed(sinput)
saveRDS(
    tpsautils::sim_data(n_1 = n, n_2 = n,
                        relation = rel,
                        r = size, points_chull = shape,
                        bbox = box_sim),
    file = sprintf("data/raw/indep_%s.rds",
                   formatC(sinput, width = 4,
                           flag = "0"))
)

##--- simulating (moderate) repulsion data ----

rel <- "repulsion"
m_hc <- seq(from = 0,
            to   = tpsautils::max_hc(n, box_sim),
            length.out = 8)[6]

set.seed(sinput)
saveRDS(
    tpsautils::sim_data(n_1 = n, n_2 = n,
                        relation = rel,
                        r = size, points_chull = shape,
                        max_hc = m_hc,
                        bbox = box_sim),
    file = sprintf("data/raw/rep_%s.rds",
                   formatC(sinput, width = 4,
                           flag = "0"))
)

##--- simulating (moderate) attraction data ----

rel <- "attraction"
r_at <- seq(from = .8,
            to   = 1e-4,
            length.out = 8)[7]

set.seed(sinput)
saveRDS(
    tpsautils::sim_data(n_1 = n, n_2 = n,
                        relation = rel,
                        r = size, points_chull = shape,
                        r_att = r_at,
                        bbox = box_sim),
    file = sprintf("data/raw/att_%s.rds",
                   formatC(sinput, width = 4,
                           flag = "0"))
)
