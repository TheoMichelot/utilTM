
mahalanobis_pairwise <- function(data) {
    out <- lapply(1:nrow(data), function(i) {
        mahalanobis(x = data, 
                    center = do.call("c", data[i, ]),
                    cov = cov(data))
    })
    return(as.dist(do.call("rbind", out)))
}
