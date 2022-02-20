compute_lavan_pvals <- function(sp_lst, bb,
                                return_ts = FALSE,
                                verbose = FALSE,
                                path_gcops = "lavancier/bin") {
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
    
    sh_cmd <- sprintf("%s -i1 %s -i2 %s",
                      file.path(path_gcops, "GcoPS"),
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
    
    if(return_ts)
        return(c(lavan_ts, lavan_pv))
    else
        return(lavan_pv)
}
