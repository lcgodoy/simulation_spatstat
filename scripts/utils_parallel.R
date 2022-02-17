compute_pvals <- function(sp_lst, bb, sid, scen, verbose = TRUE) {
    rst_sim <- lapply(sp_lst,
                      function(x) {
                          ## ncol and nrow will control the raster resolution
                          r <- raster::raster(ncol = 450,
                                              nrow = 450)
                          raster::extent(r) <- raster::extent(x)
                          raster::rasterize(x = x, y = r,
                                            getCover = TRUE)
                      })
    
    tmp_dir <- tempdir()
    tmp_files <- vector(mode   = "character",
                        length = 2L)
    
    tmp_files <- sapply(tmp_files,
                        function(x)
                            tempfile(fileext = ".tif",
                                     tmpdir  = tmp_dir))
    
    raster::writeRaster(x = rst_sim[[1]], filename = tmp_files[1],
                        format = "GTiff")
    raster::writeRaster(x = rst_sim[[2]], filename = tmp_files[2],
                        format = "GTiff")
    
    sh_cmd <- sprintf("/opt/GcoPS/bin/GcoPS -i1 %s -i2 %s -nbt 1",
                      tmp_files[1], tmp_files[2])

    ## check OS
    if(.Platform$OS.type != "unix")
        lavan_ts <- shell(command = sh_cmd, intern = TRUE)
    else
        lavan_ts <- system(command = sh_cmd, intern = TRUE)
        
    invisible(lapply(tmp_files, unlink))

    if(verbose) {
        lavan_ts
        cat("\n")
        stringr::str_extract(string  = lavan_ts[5],
                             pattern = "-?[\\d,.]+")
    }

    ## extracting test statistic from output
    lavan_ts <- as.numeric(
        stringr::str_extract(string  = lavan_ts[5],
                             pattern = "-?[\\d,.]+")
    ) 

    ## p-value according to the asymptotic dist
    lavan_pv <- 2 * pnorm(abs(lavan_ts),
                          lower.tail = FALSE)
    
    mc_pv <- tpsa::gof_mc(sp_lst[[1]], sp_lst[[2]], n_sim = 499L,
                          unique_bbox = box_sim, alpha = 0.05,
                          H = 'L', ts = 'SMAD', distances = NULL,
                          fixed = FALSE, method = 'hausdorff')$p_value

    out <- data.frame(sim_id    = sid,
                      scenario  = scen,
                      lavancier = lavan_pv,
                      mc_test   = mc_pv)


    return(out)
}
