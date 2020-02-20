# Convert epsg to epsg KM ------------------------------------------------------
epsgKM <- function(x) paste(paste0("+init=epsg:", x), "+units=km")


# Envelope for variogram -------------------------------------------------------
variog_envelope <- function (geodata, coords = geodata$coords, 
                             data = geodata$data, 
                             obj.variog, nsim = 99, save.sim = FALSE, messages) 
{
  call.fc <- match.call()
  if (missing(geodata)) 
    geodata <- list(coords = coords, data = data)
  if (missing(messages)) 
    messages.screen <- as.logical(ifelse(is.null(getOption("geoR.messages")), 
                                         TRUE, getOption("geoR.messages")))
  else messages.screen <- messages
  obj.variog$v <- NULL
  if ((is.matrix(data) | is.data.frame(data))) 
    if (ncol(data) > 1) 
      stop("envelops can be computed for only one data set at once")
  if (!is.null(obj.variog$estimator.type)) 
    estimator.type <- obj.variog$estimator.type
  else estimator.type <- "classical"
  if (abs(obj.variog$lambda - 1) > 1e-04) {
    if (abs(obj.variog$lambda) < 1e-04) 
      data <- log(data)
    else data <- ((data^obj.variog$lambda) - 1)/obj.variog$lambda
  }
  xmat <- unclass(trend.spatial(trend = obj.variog$trend, geodata = geodata))
  if (obj.variog$trend != "cte") {
    if (is.vector(data)) {
      data <- lm(data ~ xmat + 0)$residuals
      names(data) <- NULL
    }
    else {
      only.res <- function(y, x) {
        lm(y ~ xmat + 0)$residuals
      }
      data <- apply(data, 2, only.res, x = xmat)
    }
  }
  if (messages.screen) 
    cat(paste("variog.env: generating", nsim, "simulations by permutating data values\n"))
  simula <- list(coords = coords)
  n.data <- length(data)
  perm.f <- function(i, data, n.data) {
    return(data[sample(1:n.data)])
  }
  simula$data <- apply(as.matrix(1:nsim), 1, perm.f, data = data, 
                       n.data = n.data)
  if (messages.screen) 
    cat(paste("variog.env: computing the empirical variogram for the", 
              nsim, "simulations\n"))
  nbins <- length(obj.variog$bins.lim) - 1
  if (obj.variog$direction == "omnidirectional") {
    bin.f <- function(sim) {
      cbin <- vbin <- sdbin <- rep(0, nbins)
      temp <- .C("binit", as.integer(obj.variog$n.data), 
                 as.double(as.vector(coords[, 1])), as.double(as.vector(coords[, 
                                                                               2])), as.double(as.vector(sim)), as.integer(nbins), 
                 as.double(as.vector(obj.variog$bins.lim)), as.integer(estimator.type == 
                                                                         "modulus"), as.double(max(obj.variog$u)), as.double(cbin), 
                 vbin = as.double(vbin), as.integer(FALSE), as.double(sdbin), 
                 PACKAGE = "geoR")$vbin
      return(temp)
    }
    simula.bins <- apply(simula$data, 2, bin.f)
  }
  else {
    variog.vbin <- function(x, ...) {
      variog(geodata = geodata, 
             data = x, uvec = obj.variog$uvec, estimator.type = obj.variog$estimator.type, 
             nugget.tolerance = obj.variog$nugget.tolerance, max.dist = obj.variog$max.dist, 
             pairs.min = obj.variog$pairs.min, direction = obj.variog$direction, 
             tolerance = obj.variog$tolerance, messages.screen = FALSE,...)$v
    }
    simula.bins <- apply(simula$data, 2, variog.vbin)
  }
  simula.bins <- simula.bins[obj.variog$ind.bin, ]
  if (save.sim == FALSE) 
    simula$data <- NULL
  if (messages.screen) 
    cat("variog.env: computing the envelops\n")
  limits <- apply(simula.bins, 1, quantile, prob = c(0.025, 0.975))
  res.env <- list(u = obj.variog$u, v.lower = limits[1, ], 
                  v.upper = limits[2, ])
  if (save.sim) 
    res.env$simulations <- simula$data
  res.env$call <- call.fc
  oldClass(res.env) <- "variogram.envelope"
  return(res.env)
}

# Calculate and plot the variogram ---------------------------------------------
ggvario <- function(coords, 
                    data, 
                    bins = 15, 
                    maxdist = max(dist(coords))/3, 
                    uvec = NULL, 
                    nsim = 999,
                    color = "royalblue1", 
                    xlab = "distance", 
                    show_nbins = T) {
  require(ggplot2)
  require(geoR)
  coords <- as.matrix(coords)
  min_dist <- min(dist(coords))
  if(is.null(uvec)) uvec <- seq(min_dist, maxdist, l = bins)
  empvario <- variog(coords = coords, data = data, uvec = uvec, messages = F)
  envmc <- variog_envelope(coords = coords, data = data, 
                           obj.variog = empvario, nsim = nsim, messages = F)
  dfvario <- data.frame(distance = empvario$u, empirical = empvario$v,
                        lowemp = envmc$v.lower, upemp = envmc$v.upper, 
                        nbins = empvario$n)
  p1 <- ggplot(dfvario, aes(y = empirical, x = distance, label = nbins)) +
    geom_ribbon(aes(ymin = lowemp, ymax = upemp), fill = color, alpha = .3) +
    geom_point(aes(y = empirical), col = "black", fill = color, shape = 21, size = 3) +
    scale_x_continuous(name = xlab, limits = c(0, uvec[length(uvec)]),
                       breaks = round(seq(0, uvec[length(uvec)], l = 6))) +
    scale_y_continuous(name = "semivariance", 
                       #breaks = round(seq(0, max(dfvario$upemp, dfvario$empirical), l = 6), 1), 
                       limits = c(0, max(dfvario$upemp, dfvario$empirical))) +
    ggtitle("Empirical semivariogram") +
    theme_classic()
  p2 <- p1 + geom_text(vjust = 1, nudge_y = - diff(range(dfvario$empirical)) / 22)
  if(show_nbins) p2 else p1
}

# Empirical logit --------------------------------------------------------------
elogit <- function(num, den) {
  log((num + 0.5) / (den - num + 0.5))
}

# Check duplicated rows --------------------------------------------------------
check_dups <- function(cols, data) {
  data$row_id <- 1:nrow(data)
  dups <- data %>%
    group_by_at(cols) %>% 
    filter(n() > 1)  %>% 
    ungroup()
  if(nrow(dups) == 0) {
    dups <- NULL
  } else {
    dups$dup_id <- unclass(factor(apply(dups[, cols], 1, paste0, collapse = "")))
    dups <- dups[order(dups$dup_id), ]
  }
  return(dups)
}


# Calculate and extract UTM km coordinates from data ---------------------------

extract_coords <- function(x) {
  x[, c("utm_x", "utm_y")] <- st_coordinates(st_transform(x, crs = epsgKM(unique(x$crs))))
  x$geometry <- NULL
  x$crs <- NULL
  as.data.frame(x)
}

# Fit binomial geostatistical models to list of datasets -----------------------

fitALL <- function(data, path) {
  
  all_countries <- unique(data$country)
  data_list <- split(data, f = data$country)
  
  # Fit a linear geostatistical model to the empirical logit to get starting
  # value for the variance covariance parameters
  
  map(data_list, safely(function(x) {
    x <- as.data.frame(x)
    country <- unique(x$country)
    cat("Analysing", nrow(x), "survey data from", country, "(",
        match(country, all_countries), "out of", length(all_countries), " countries)\n\n") 
    phi_start <- median(dist(x[, c("utm_x", "utm_y")])) / 3
    f <- elogit ~ 1
    cat("Fitting Linear Geostatistical Model to", country, "data\n")
    geo_linear <- linear.model.MLE(formula = f,
                                   coords =  ~ utm_x + utm_y,
                                   data = x, start.cov.pars = c(phi_start, 1),
                                   kappa = .5, messages = T, method = "BFGS")
    
    # Fitting binomial geostatistical model
    cmcmc <- list()
    cmcmc[[1]] <- control.mcmc.MCML(n.sim = 10000, burnin = 2000, thin = 8)
    cmcmc[[2]] <- control.mcmc.MCML(n.sim = 10000, burnin = 2000, thin = 8)
    cmcmc[[3]] <- control.mcmc.MCML(n.sim = 65000, burnin = 5000, thin = 6)
    
    par0 <- as.numeric(coef(geo_linear))
    p <- length(par0) - 3
    f <- positive ~ 1
    niter <- 3
    for(i in 1:niter) {
      cat("\nFitting Geostatistical Binomial Model to", country, "data:",
          "Iteration", i, "of", 3, "\n")
      init_cov_pars <- c(par0[p + 2], par0[p + 3] / par0[p + 1])
      geo_binomial <- binomial.logistic.MCML(formula = f, units.m = ~ examined,
                                             coords = ~ utm_x + utm_y,
                                             data = x,
                                             par0 = par0,
                                             start.cov.pars = init_cov_pars,
                                             control.mcmc = cmcmc[[i]], 
                                             kappa = 0.5,
                                             method = "nlminb", messages = T)
      par0 <- as.numeric(coef(geo_binomial))
      llik <- geo_binomial$log.lik
      cat("Estimated parameters", par0, "Log-likelihood:", llik, "\n\n")
    }
    cat("Saving final results for", country, "\n\n")
    fpath <- paste0(path, "geo_binomial_", gsub(" ", "_", country), ".rds")
    geo_binomial$country <- country
    saveRDS(geo_binomial, fpath)
  }))
}

ggdist <- function(coords) {
  require(ggplot2)
  dd <- sort(as.numeric(dd <- dist(x[, c("utm_x", "utm_y")], diag = F, upper = F)))
  df <- data.frame(prop = 1:length(dd) / length(dd), dists = dd)
  ggplot(df, aes(x = dists, y = prop)) +
    geom_line() +
    scale_y_continuous(labels = scales::percent) +
    labs(x = "Distance between pairs of points (km)", 
         y = "Proportions of pairs of points within distance x")
    
}

extract_params <- function(x) {
  fit <- readRDS(x)
  params <- as.numeric(PrevMap:::coef.PrevMap(fit))
  n <- length(fit$y)
  return(c(n, params))
}
