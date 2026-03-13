# Custom functions used for analysis and figure production

# Original code developed by:
#   Harvey, B.J., R.A. Andrus, and S.C. Anderson. 2019. Incorporating biophysical gradients and 
#   uncertainty into burn severity maps in a temperate fire‐prone forested region. Ecosphere 10.
# Code modified by M. Buonanduci


# Helper functions
rescale <- function(x) {
  x.obs <- x[!is.na(x)]
  (x - mean(x.obs))/(2 * stats::sd(x.obs))
}
log_sum_exp <- function(u) {
  max_u <- max(u)
  max_u + log(sum(exp(u - max_u)))
}
# and then include the −logM term to make it log_mean_exp:
log_mean_exp <- function(u) {
  M <- length(u)
  -log(M) + log_sum_exp(u)
}


# Format data: generate y, x; scale predictors
get_dat <- function(d, # data frame
                    response, # response variable
                    main_predictor, # primary predictor (i.e., satellite index)
                    other_predictors, # additional predictors (i.e., modulating covariates)
                    scale_other_predictors) { # logical T/F for whether to center & scale predictors; same length as other_predictors
  
  if(length(other_predictors) != length(scale_other_predictors)){
    warning("Length of scale_other_predictors should match length of other_predictors.")
  }
  
  d <- d %>%
    select( all_of(c(response, main_predictor, "Fire_ID", other_predictors)) ) %>%
    drop_na(1) %>% drop_na(2) %>% # drop NA values in response and main predictor columns
    rename(y = 1, x = 2, group = 3) %>%
    mutate(
      y0 = ifelse(y == 0, 1, 0),
      y1 = ifelse(y == 1 & y != 0, 1, ifelse(y < 1 & y != 0, 0, NA)),
      yp = ifelse(y == 0 | y == 1, NA, y),
      xraw_mean = mean(x),
      xraw_sd = sd(x),
      xscaled = rescale(x),
      group = as.factor(tolower(group)),
      group_id = as.integer(group)
    ) %>%
    arrange(group_id, x)
  
  names(d) <- tolower(names(d))
  other_predictors <- tolower(other_predictors)
  
  if (sum(scale_other_predictors) > 0) {
    d <- d %>%
      dplyr::mutate_at(
        other_predictors[scale_other_predictors],
        function(x) rescale(x)
      )
  }
  d
}

# Format data for stan modeling; called internally by fit_model()
prep_stan_dat <- function(d, predictors = "xscaled", arrange_by_group = TRUE) {
  
  # in case some missing now (e.g. cross-validation):
  d <- d %>% mutate(
    group = as.factor(as.character(group)),
    group_id = as.integer(group)
  )
  
  if (arrange_by_group) {
    d <- d %>% arrange(group_id, x)
  }
  d1 <- filter(d, !is.na(y1))
  dp <- filter(d, !is.na(yp))
  
  N0 <- nrow(d)
  N1 <- nrow(d1)
  Np <- nrow(dp)
  
  f0 <- as.formula(paste("y0 ~", paste(predictors, collapse = " + ")))
  f1 <- as.formula(paste("y1 ~", paste(predictors, collapse = " + ")))
  fp <- as.formula(paste("y ~", paste(predictors, collapse = " + ")))
  mm0 <- as.matrix(model.matrix(f0, data = d))
  mm1 <- as.matrix(model.matrix(f1, data = d1))
  mmp <- as.matrix(model.matrix(fp, data = dp))
  
  sd <- list(
    N0 = N0,
    N1 = N1,
    Np = Np,
    J01 = ncol(mm0),
    Jp = ncol(mmp),
    X0_ij = mm0,
    X1_ij = mm1,
    Xp_ij = mmp,
    y0_i = d$y0,
    y1_i = d1$y1,
    yp_i = dp$y,
    
    Ng = length(unique(as.character(d$group))),
    group_id0 = d$group_id,
    group_id1 = d1$group_id,
    group_idp = dp$group_id
  )
  
  list(stan_dat = sd, f = fp)
}

# Fit ZOIB model
fit_model <- function(d, model, predictors = "xscaled",
                      re = TRUE, # include group-level random effects? (T/F)
                      iter = 800L, chains = 4L,
                      adapt_delta = 0.8, log_lik = TRUE, ...) {
  
  pars <- c("b0_j", "b1_j", "bp_j", "phi")
  if (re) pars <- c(pars, c("sigma_z", "z0_g", "z1_g", "zp_g"))
  if (log_lik) pars <- c(pars, "log_lik")
  prep <- prep_stan_dat(d, predictors = predictors)
  m <- sampling(model,
                data = prep$stan_dat,
                iter = iter, chains = chains,
                pars = pars,
                control = list(adapt_delta = adapt_delta, max_treedepth = 20), ...
  )
  list(model = m, data = d, f = prep$f, stan_dat = prep$stan_dat)
}



# Make predictions, drawing the EXPECTED VALUE for y from the posterior predictive distribution
make_predictions_exp <- function(d, f, model,
                                 use_new_data = TRUE, # create new dataset for prediction? (T/F)
                                 Npred = 200L, # number of predicted values if above is TRUE
                                 re = TRUE, est_range = FALSE) {
  if (use_new_data) {
    d_pred <- expand.grid(
      xscaled = seq(min(d$xscaled), max(d$xscaled), length.out = Npred),
      group_id = unique(d$group_id), y = 1
    )
  } else {
    d_pred <- d
  }
  
  draws <- as_draws_rvars(as.matrix(model))
  
  # model matrix
  mm_pred <- as.matrix(model.matrix(f, data = d_pred))
  
  # multiply intercept
  pred_p_mu <- mm_pred[,1] * draws$bp_j[1]
  pred_0 <- mm_pred[,1] * draws$b0_j[1]
  pred_1 <- mm_pred[,1] * draws$b1_j[1]
  
  # multiply and add covariate(s)
  for( i in 2:ncol(mm_pred) ){
    pred_p_mu <- pred_p_mu + mm_pred[,i] * draws$bp_j[i]
    pred_0 <- pred_0 + mm_pred[,i] * draws$b0_j[i]
    pred_1 <- pred_1 + mm_pred[,i] * draws$b1_j[i]
  }
  
  # random effects
  if (re) {
    zp_g <- draws$zp_g[d_pred$group_id]
    z0_g <- draws$z0_g[d_pred$group_id]
    z1_g <- draws$z1_g[d_pred$group_id]
    
    # convert predictions w/ re to response scale, add to dataframe
    pred_p_mu_re <- as_draws_matrix(pred_p_mu + zp_g) %>%
      plogis() %>% as_draws_rvars()
    pred_0_re <- as_draws_matrix(pred_0 + z0_g) %>%
      plogis() %>% as_draws_rvars()
    pred_1_re <- as_draws_matrix(pred_1 + z1_g) %>%
      plogis() %>% as_draws_rvars()
    
    d_pred <- d_pred %>%
      mutate(pred_p_mu_re = pred_p_mu_re$x,
             pred_0_re = pred_0_re$x,
             pred_1_re = pred_1_re$x) %>%
      mutate(y_pred_re = (1 - pred_0_re) * (pred_1_re + (1 - pred_1_re) * pred_p_mu_re) )
  }
  
  # convert predictions w/o re to response scale, add to dataframe
  pred_p_mu <- as_draws_matrix(pred_p_mu) %>%
    plogis() %>% as_draws_rvars()
  pred_0 <- as_draws_matrix(pred_0) %>%
    plogis() %>% as_draws_rvars()
  pred_1 <- as_draws_matrix(pred_1) %>%
    plogis() %>% as_draws_rvars()
  
  d_pred <- d_pred %>%
    mutate(pred_p_mu = pred_p_mu$x,
           pred_0 = pred_0$x,
           pred_1 = pred_1$x) %>%
    mutate(y_pred = (1 - pred_0) * (pred_1 + (1 - pred_1) * pred_p_mu) )
  
  pred_df <- tibble(
    xscaled = d_pred$xscaled,
    x =  d_pred$xscaled * 2 * unique(d$xraw_sd) + unique(d$xraw_mean),
    group = d_pred$group_id,
    est = d_pred$y_pred
  )
  
  if (re) {
    pred_df$est_re <- d_pred$y_pred_re
  }
  
  # if condition is true, return median and 90% CI rather than full rvar object
  if (est_range) {
    pred_df <- pred_df %>%
      mutate(est_lwr = t(quantile(est, 0.05))[,1],
             est_upr = t(quantile(est, 0.95))[,1],
             est = t(quantile(est, 0.5))[,1])
    if (re) {
      pred_df <- pred_df %>%
        mutate(est_re_lwr = t(quantile(est_re, 0.05))[,1],
               est_re_upr = t(quantile(est_re, 0.95))[,1],
               est_re = t(quantile(est_re, 0.5))[,1])
    }
  }
  pred_df
}


# Makes predictions, drawing new observations for y from the posterior predictive distribution
# By definition, these draws have higher variance than draws of the expected value of the posterior predictive distribution
# because they include aleatoric uncertainty (i.e., sampling uncertainty);
# see https://mc-stan.org/docs/stan-users-guide/posterior-prediction.html#sampling-from-the-posterior-predictive-distribution
make_predictions_yrep <- function(d, f, model,
                                  use_new_data = TRUE, # create new dataset for prediction? (T/F)
                                  Npred = 200L, # number of predicted values if above is TRUE
                                  re = TRUE, est_range = FALSE) {
  if (use_new_data) {
    d_pred <- expand.grid(
      xscaled = seq(min(d$xscaled), max(d$xscaled), length.out = Npred),
      group_id = unique(d$group_id), y = 1
    )
  } else {
    d_pred <- d
  }
  
  draws <- as_draws_rvars(as.matrix(model))
  
  # model matrix
  mm_pred <- as.matrix(model.matrix(f, data = d_pred))
  
  # multiply intercept
  pred_p_mu <- mm_pred[,1] * draws$bp_j[1]
  pred_0 <- mm_pred[,1] * draws$b0_j[1]
  pred_1 <- mm_pred[,1] * draws$b1_j[1]
  
  # multiply and add covariate(s)
  for( i in 2:ncol(mm_pred) ){
    pred_p_mu <- pred_p_mu + mm_pred[,i] * draws$bp_j[i]
    pred_0 <- pred_0 + mm_pred[,i] * draws$b0_j[i]
    pred_1 <- pred_1 + mm_pred[,i] * draws$b1_j[i]
  }
  
  # random effects
  if (re) {
    zp_g <- draws$zp_g[d_pred$group_id]
    z0_g <- draws$z0_g[d_pred$group_id]
    z1_g <- draws$z1_g[d_pred$group_id]
    
    # convert predictions w/ re to response scale, add to dataframe
    pred_p_mu_re <- as_draws_matrix(pred_p_mu + zp_g) %>%
      plogis() %>% as_draws_rvars()
    pred_0_re <- as_draws_matrix(pred_0 + z0_g) %>%
      plogis() %>% as_draws_rvars()
    pred_1_re <- as_draws_matrix(pred_1 + z1_g) %>%
      plogis() %>% as_draws_rvars()
    
    d_pred <- d_pred %>%
      mutate(pred_p_mu_re = pred_p_mu_re$x, # mu for beta
             pred_0_re = pred_0_re$x, # p0
             pred_1_re = pred_1_re$x, # p1
             phi = draws$phi,
             shape1_re = pred_p_mu_re * phi,
             shape2_re = (1 - pred_p_mu_re) * phi,
             y_pred_p_re = rvar_rng(rbeta, n(),
                                    shape1 = shape1_re, shape2 = shape2_re),
             y_pred_0_re = rvar_rng(rbinom, n(),
                                    size = 1, prob = pred_0_re),
             y_pred_1_re = rvar_rng(rbinom, n(),
                                    size = 1, prob = pred_1_re)) %>%
      mutate(y_pred_re = rvar_ifelse(y_pred_0_re == 1, 0,
                                     rvar_ifelse(y_pred_1_re == 1, 1,
                                                 y_pred_p_re)))
  }
  
  # convert predictions w/o re to response scale, add to dataframe
  pred_p_mu <- as_draws_matrix(pred_p_mu) %>%
    plogis() %>% as_draws_rvars()
  pred_0 <- as_draws_matrix(pred_0) %>%
    plogis() %>% as_draws_rvars()
  pred_1 <- as_draws_matrix(pred_1) %>%
    plogis() %>% as_draws_rvars()
  
  d_pred <- d_pred %>%
    mutate(pred_p_mu = pred_p_mu$x, # mu for beta
           pred_0 = pred_0$x, # p0
           pred_1 = pred_1$x, # p1
           phi = draws$phi,
           shape1 = pred_p_mu * phi,
           shape2 = (1 - pred_p_mu) * phi,
           y_pred_p = rvar_rng(rbeta, n(),
                               shape1 = shape1, shape2 = shape2),
           y_pred_0 = rvar_rng(rbinom, n(),
                               size = 1, prob = pred_0),
           y_pred_1 = rvar_rng(rbinom, n(),
                               size = 1, prob = pred_1)) %>%
    mutate(y_pred = rvar_ifelse(y_pred_0 == 1, 0,
                                rvar_ifelse(y_pred_1 == 1, 1,
                                            y_pred_p)))
  
  pred_df <- tibble(
    xscaled = d_pred$xscaled,
    x = d_pred$xscaled * 2 * unique(d$xraw_sd) + unique(d$xraw_mean),
    group = d_pred$group_id,
    est = d_pred$y_pred
  )
  
  if (re) {
    pred_df$est_re <- d_pred$y_pred_re
  }
  
  # if condition is true, return median and 90% CI rather than full rvar object
  if (est_range) {
    pred_df <- pred_df %>%
      mutate(est_lwr = t(quantile(est, 0.05))[,1],
             est_upr = t(quantile(est, 0.95))[,1],
             est = t(quantile(est, 0.5))[,1])
    if (re) {
      pred_df <- pred_df %>%
        mutate(est_re_lwr = t(quantile(est_re, 0.05))[,1],
               est_re_upr = t(quantile(est_re, 0.95))[,1],
               est_re = t(quantile(est_re, 0.5))[,1])
    }
  }
  pred_df
}


# Plot posterior median predictions without fire-level random effects
plot_post_med <- function(d, m, 
                          type = "exp", # type of predictions; either expected value ("exp") or posterior prediction ("yrep")
                          xlab = "x", ylab = "Proportion", title = "", 
                          return_data = FALSE) # return dataframe rather than plot? (T/F)
{
  
  range_xscaled <- range(d$xscaled)
  newdata <- d
  cn <- colnames(m$stan_dat$Xp_ij)
  cn <- cn[!grepl("\\(", cn)]
  cn <- cn[!grepl("\\:", cn)]
  cn <- cn[!grepl("xscaled", cn)]
  newdata[, cn] <- 0
  newdata[, "xscaled"] <- seq(range_xscaled[1],
                              range_xscaled[2],
                              length.out = nrow(newdata) )
  
  if(type == "exp"){
    p <- make_predictions_exp(
      d = newdata, f = m$f,
      model = m$model, re = FALSE, use_new_data = FALSE)
    
  } else if(type == "yrep"){
    p <- make_predictions_yrep(
      d = newdata, f = m$f,
      model = m$model, re = FALSE, use_new_data = FALSE)
    
  } else {
    warning("'type' must equal either 'exp' or 'yrep'")
  }
  
  if (!return_data) {
    g <- ggplot() +
      geom_point(data = d, mapping = aes(x, y), alpha = 0.1) +
      geom_line(data = p, mapping = aes(x, median(est)),
                lwd = 1.5, color="steelblue") +
      geom_ribbon(data = p, mapping = aes(x, 
                                          ymin = t(quantile(est, 0.05)), 
                                          ymax = t(quantile(est, 0.95))),
                  alpha = 0.2, fill="steelblue") +
      ylim(0, 1) +
      ylab(ylab) +
      xlab(xlab) +
      coord_cartesian(expand = TRUE) +
      ggtitle(title) +
      theme_bw() +
      theme(panel.grid = element_blank())
    return(g)
  } else {
    return(as_tibble(p))
  }
}

# Plot effects of primary predictor and one interacting covariate
plot_interaction <- function(d, m,
                             int_var, # name of interaction variable
                             int_var_cont, # is interaction variable continuous (TRUE) or categorical (FALSE)?
                             .quant = c(0.9, 0.1), # quantiles of interaction variable, if continuous
                             type = "exp", # type of predictions; either expected value ("exp") or posterior prediction ("yrep")
                             title = "", xlab = "x", ylab = "Proportion",
                             return_data = FALSE) # return dataframe rather than plot? (T/F)
{
  
  int_var <- tolower(int_var)
  newdata <- d
  cn <- colnames(m$stan_dat$Xp_ij)
  cn <- cn[!grepl("\\(", cn)]
  cn <- cn[!grepl("\\:", cn)]
  cn <- cn[!grepl("xscaled", cn)]
  newdata[, cn] <- 0
  newdata[, "xscaled"] <- seq(min(newdata[, "xscaled"]),
                              max(newdata[, "xscaled"]),
                              length.out = nrow(newdata) )
  
  if (int_var_cont) {
    newdata[, int_var] <- quantile(d[, int_var, drop = TRUE], .quant[[1]])
    newdata$level <- as.character(.quant[[1]])
    newdata1 <- newdata
    
    for (i in seq(2, length(.quant))) {
      newdata2 <- newdata1
      newdata2[, int_var] <- quantile(d[, int_var, drop = TRUE], .quant[[i]])
      newdata2$level <- as.character(.quant[[i]])
      newdata <- bind_rows(newdata, newdata2)
    }
    
  } else {
    int_var_vals <- unique(d[int_var])
    newdata[, int_var] <- int_var_vals[1,]
    newdata$level <- as.character(int_var_vals[1,])
    newdata1 <- newdata
    
    for (i in seq(2, nrow(int_var_vals))) {
      newdata2 <- newdata1
      newdata2[, int_var] <- int_var_vals[i,]
      newdata2$level <- as.character(int_var_vals[i,])
      newdata <- bind_rows(newdata, newdata2)
    }
  }
  
  if(type == "exp"){
    p <- make_predictions_exp(
      d = newdata, f = m$f,
      model = m$model, re = FALSE, use_new_data = FALSE)
    
  } else if(type == "yrep"){
    p <- make_predictions_yrep(
      d = newdata, f = m$f,
      model = m$model, re = FALSE, use_new_data = FALSE)
    
  } else {
    warning("'type' must equal either 'exp' or 'yrep'")
  }
  
  p$level <- newdata$level
  
  if (!return_data) {
    g <- ggplot(p, aes(x,
                       y = median(est),
                       ymin = t(quantile(est, 0.05)),
                       ymax = t(quantile(est, 0.95)),
                       group = level,
                       fill = level
    )) +
      geom_point(data = d, mapping = aes(x, y), color = "gray", alpha = 0.2, inherit.aes = FALSE) +
      geom_line(lwd = 1.5, aes(color = level)) +
      geom_ribbon(alpha = 0.3) +
      scale_fill_viridis_d(option = 'mako', direction = -1, begin = 0.2, end = 0.8) +
      scale_color_viridis_d(option = 'mako', direction = -1, begin = 0.2, end = 0.8) +
      guides(color = guide_legend(reverse = T),
             fill = guide_legend(reverse = T)) +
      ylim(0, 1) +
      ylab(ylab) +
      xlab(xlab) +
      labs(
        color = paste0(ifelse(int_var_cont, "Quantile:\n", "Level:\n"), int_var),
        fill = paste0(ifelse(int_var_cont, "Quantile:\n", "Level:\n"), int_var)
      ) +
      coord_cartesian(expand = TRUE) +
      ggtitle(title) +
      theme_bw() +
      theme(panel.grid = element_blank())
    return(g)
  } else {
    return(as_tibble(p))
  }
}

