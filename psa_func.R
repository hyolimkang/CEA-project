function (params = NULL, dists = c("normal", "log-normal", "truncated-normal", 
                                   "beta", "gamma", "dirichlet", "bootstrap", "constant", "triangle"), 
          parameterization_types = c("mean, sd", "a, b", "shape, scale", 
                                     "value, mean_prop, sd", "value, n", "value, alpha", 
                                     "mean, sd, ll, ul", "val", "meanlog, sdlog", "ll, ul, mode"), 
          dists_params = NULL, nsamp = 100) 
{
  dists <- match.arg(dists, several.ok = TRUE)
  parameterization_types <- match.arg(parameterization_types, 
                                      several.ok = TRUE)
  n_params <- length(params)
  params_df <- vector(mode = "list", length = n_params)
  for (i in 1:n_params) {
    if (dists[i] == "normal") {
      params_df[[i]] <- data.frame(param_val = rnorm(nsamp, 
                                                     mean = dists_params[[i]][1], sd = dists_params[[i]][2]))
      names(params_df[[i]]) <- paste0(params[i])
    }
    if (dists[i] == "log-normal") {
      if (parameterization_types[i] == "mean, sd") {
        mu <- lnorm_params(dists_params[[i]][1], (dists_params[[i]][2])^2)[[1]]
        sd <- lnorm_params(dists_params[[i]][1], (dists_params[[i]][2])^2)[[2]]
      }
      else if (parameterization_types[i] == "meanlog, sdlog") {
        mu <- dists_params[[i]][1]
        sd <- dists_params[[i]][2]
      }
      params_df[[i]] <- data.frame(param_val = rlnorm(nsamp, 
                                                      meanlog = mu, sdlog = sd))
      names(params_df[[i]]) <- paste0(params[i])
    }
    if (dists[i] == "truncated-normal") {
      sample_mean <- dists_params[[i]][1]
      sample_sd <- dists_params[[i]][2]
      lowerbound <- ifelse(!is.na(dists_params[[i]][3]), 
                           dists_params[[i]][3], -Inf)
      upperbound <- ifelse(!is.na(dists_params[[i]][4]), 
                           dists_params[[i]][4], Inf)
      params_df[[i]] <- data.frame(param_val = rtruncnorm(nsamp, 
                                                          a = lowerbound, b = upperbound, mean = sample_mean, 
                                                          sd = sample_sd))
      names(params_df[[i]]) <- paste0(params[i])
    }
    if (dists[i] == "beta") {
      if (parameterization_types[i] == "mean, sd") {
        a <- beta_params(dists_params[[i]][1], dists_params[[i]][2])[[1]]
        b <- beta_params(dists_params[[i]][1], dists_params[[i]][2])[[2]]
        params_df[[i]] <- as.data.frame(rbeta(nsamp, 
                                              a, b))
      }
      else if (parameterization_types[i] == "a, b") {
        a <- dists_params[[i]][1]
        b <- dists_params[[i]][2]
        params_df[[i]] <- as.data.frame(rbeta(nsamp, 
                                              a, b))
      }
      names(params_df[[i]]) <- paste0(params[i])
    }
    if (dists[i] == "gamma") {
      if (parameterization_types[i] == "mean, sd") {
        shape <- gamma_params(dists_params[[i]][1], 
                              dists_params[[i]][2], scale = TRUE)[[1]]
        scale <- gamma_params(dists_params[[i]][1], 
                              dists_params[[i]][2], scale = TRUE)[[2]]
        params_df[[i]] <- as.data.frame(rgamma(nsamp, 
                                               shape = shape, scale = scale))
      }
      else if (parameterization_types[i] == "shape, scale") {
        shape <- dists_params[[i]][1]
        scale <- dists_params[[i]][2]
        params_df[[i]] <- as.data.frame(rgamma(nsamp, 
                                               shape = shape, scale = scale))
      }
      names(params_df[[i]]) <- paste0(params[i])
    }
    if (dists[i] == "dirichlet") {
      if (parameterization_types[i] == "value, mean_prop, sd") {
        alpha <- dirichlet_params(dists_params[[i]][, 
                                                    2], dists_params[[i]][, 3])
        params_df[[i]] <- as.data.frame(rdirichlet(nsamp, 
                                                   alpha))
      }
      else if (parameterization_types[i] == "value, n") {
        val_n <- as.data.frame(dists_params[[i]])
        total <- sum(val_n[, 2])
        p_mean <- val_n[, 2]/total
        sd <- sqrt((p_mean * (1 - p_mean))/total)
        alpha <- dirichlet_params(p_mean, sd)
        params_df[[i]] <- as.data.frame(rdirichlet(nsamp, 
                                                   alpha))
      }
      else if (parameterization_types[i] == "value, alpha") {
        alpha <- dists_params[[i]][, 2]
        params_df[[i]] <- as.data.frame(rdirichlet(nsamp, 
                                                   alpha))
      }
      names(params_df[[i]]) <- paste0(dists_params[[i]][, 
                                                        1])
    }
    if (dists[i] == "bootstrap") {
      samp_vec <- vector(mode = "numeric", length = nsamp)
      for (k in 1:nsamp) {
        samp_vec[k] <- mean(sample(x = dists_params[[i]][, 
                                                         1], size = sum(dists_params[[i]][, 2]), replace = TRUE, 
                                   prob = dists_params[[i]][, 2]))
      }
      params_df[[i]] <- as.data.frame(samp_vec)
      names(params_df[[i]]) <- paste0(params[i])
    }
    if (dists[i] == "triangle") {
      if (parameterization_types[i] == "ll, ul, mode") {
        a <- dists_params[[i]][1]
        b <- dists_params[[i]][2]
        c <- dists_params[[i]][3]
        params_df[[i]] <- as.data.frame(rtriangle(n = nsamp, 
                                                  a = a, b = b, c = c))
      }
      names(params_df[[i]]) <- paste0(params[i])
    }
    if (dists[i] == "constant") {
      val <- dists_params[[i]]
      params_df[[i]] <- data.frame(param_val = rep(val, 
                                                   nsamp))
      names(params_df[[i]]) <- paste0(params[i])
    }
  }
  params_df <- do.call(cbind, params_df)
  nsamp <- 1:nsamp
  params_df <- cbind(nsamp, params_df)
  return(params_df)
}