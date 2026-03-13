# Fit ZOIB models to use in producing atlas
# Code developed by M. Buonanduci

# Load packages
library(here)
library(tidyverse)
library(rstan)

# Load custom functions
source(here("00_zoib_functions.R"))

# Use multiple cores
options(mc.cores = 4)

# Load Stan ZOIB model
sm_zoib <- stan_model(here("zoib.stan"),
                      save_dso = TRUE, auto_write = TRUE)

# Specify number of iterations and chains
iter <- 2000L
chains <- 4L

# Load input data
input_data <- here("input_data.csv") %>% read_csv()


# Basal area mortality ----------

# Format data
other_predictors <- c("postfire2", "TT_cwd", "RAP_tcov")
scale_other_predictors <- c(FALSE, TRUE, TRUE)

# Get dataframe with scaled predictors for model-fitting
d <- get_dat(input_data,
             response = "prop_killedBA", main_predictor = "rdnbr_w_offset",
             other_predictors = other_predictors, 
             scale_other_predictors = scale_other_predictors)

# Get dataframe with unscaled predictors
# Calculate mean and sd to use later for burn severity map predictions
d_unscaled <- get_dat(input_data,
                      response = "prop_killedBA", main_predictor = "rdnbr_w_offset",
                      other_predictors = other_predictors, 
                      scale_other_predictors = rep(FALSE, length(other_predictors))) 

d_mean <- d_unscaled %>% summarise(across(all_of(tolower(c("x", other_predictors))), mean))
d_sd <- d_unscaled %>% summarise(across(all_of(tolower(c("x", other_predictors))), sd))


# Fit model
m <- fit_model(d, model = sm_zoib, re = TRUE,
               predictors = c("xscaled", paste("xscaled", tolower(other_predictors), sep = "*")),
               iter = iter, chains = chains)

# Add predictor means and SDs to list
m$data_mean = d_mean
m$data_sd = d_sd

# Add number of post-warmup draws
m$ndraws = dim( as.matrix(m$model) )[1]
  
# Write to file
saveRDS(m, here("m_west_BA.rds"))


# Stem mortality ----------

# Format data
other_predictors <- c("postfire2")
scale_other_predictors <- c(FALSE)

# Get dataframe with scaled predictors for model-fitting
d <- get_dat(input_data,
             response = "prop_killedSTEM", main_predictor = "rdnbr_w_offset",
             other_predictors = other_predictors, 
             scale_other_predictors = scale_other_predictors)

# Get dataframe with unscaled predictors
# Calculate mean and sd to use later for burn severity map predictions
d_unscaled <- get_dat(input_data,
                      response = "prop_killedSTEM", main_predictor = "rdnbr_w_offset",
                      other_predictors = other_predictors, 
                      scale_other_predictors = rep(FALSE, length(other_predictors))) 

d_mean <- d_unscaled %>% summarise(across(all_of(tolower(c("x", other_predictors))), mean))
d_sd <- d_unscaled %>% summarise(across(all_of(tolower(c("x", other_predictors))), sd))


# Fit model
m <- fit_model(d, model = sm_zoib, re = TRUE,
               predictors = c("xscaled", paste("xscaled", tolower(other_predictors), sep = "*")),
               iter = iter, chains = chains)

# Add predictor means and SDs to list
m$data_mean = d_mean
m$data_sd = d_sd

# Add number of post-warmup draws
m$ndraws = dim( as.matrix(m$model) )[1]

# Write to file
saveRDS(m, here("m_west_stems.rds"))


# Canopy cover change ----------

# Format data
other_predictors <- c("postfire2", "RAP_tcov")
scale_other_predictors <- c(FALSE, TRUE)

# Get dataframe with scaled predictors for model-fitting
d <- get_dat(input_data,
             response = "delt_livecc", main_predictor = "rdnbr_w_offset",
             other_predictors = other_predictors, 
             scale_other_predictors = scale_other_predictors)

# Get dataframe with unscaled predictors
# Calculate mean and sd to use later for burn severity map predictions
d_unscaled <- get_dat(input_data,
                      response = "delt_livecc", main_predictor = "rdnbr_w_offset",
                      other_predictors = other_predictors, 
                      scale_other_predictors = rep(FALSE, length(other_predictors))) 

d_mean <- d_unscaled %>% summarise(across(all_of(tolower(c("x", other_predictors))), mean))
d_sd <- d_unscaled %>% summarise(across(all_of(tolower(c("x", other_predictors))), sd))


# Fit model
m <- fit_model(d, model = sm_zoib, re = TRUE,
               predictors = c("xscaled", paste("xscaled", tolower(other_predictors), sep = "*")),
               iter = iter, chains = chains)

# Add predictor means and SDs to list
m$data_mean = d_mean
m$data_sd = d_sd

# Add number of post-warmup draws
m$ndraws = dim( as.matrix(m$model) )[1]

# Write to file
saveRDS(m, here("m_west_CC.rds"))


# Dead needle index ----------

# Format data
other_predictors <- c("TT_cwd", "LF_slp")
scale_other_predictors <- c(TRUE, TRUE)

# Get dataframe with scaled predictors for model-fitting
d <- get_dat(input_data,
             response = "prop_deadneedle", main_predictor = "rdnbr_w_offset",
             other_predictors = other_predictors, 
             scale_other_predictors = scale_other_predictors)

# Get dataframe with unscaled predictors
# Calculate mean and sd to use later for burn severity map predictions
d_unscaled <- get_dat(input_data,
                      response = "prop_deadneedle", main_predictor = "rdnbr_w_offset",
                      other_predictors = other_predictors, 
                      scale_other_predictors = rep(FALSE, length(other_predictors))) 

d_mean <- d_unscaled %>% summarise(across(all_of(tolower(c("x", other_predictors))), mean))
d_sd <- d_unscaled %>% summarise(across(all_of(tolower(c("x", other_predictors))), sd))


# Fit model
m <- fit_model(d, model = sm_zoib, re = TRUE,
               predictors = c("xscaled", paste("xscaled", tolower(other_predictors), sep = "*")),
               iter = iter, chains = chains)

# Add predictor means and SDs to list
m$data_mean = d_mean
m$data_sd = d_sd

# Add number of post-warmup draws
m$ndraws = dim( as.matrix(m$model) )[1]

# Write to file
saveRDS(m, here("m_west_needle.rds"))


# Bole char (circumference) ----------

# Format data
other_predictors <- c("TT_cwd")
scale_other_predictors <- c(TRUE)

# Get dataframe with scaled predictors for model-fitting
d <- get_dat(input_data,
             response = "prop_meanbolescorch", main_predictor = "rdnbr_w_offset",
             other_predictors = other_predictors, 
             scale_other_predictors = scale_other_predictors)

# Get dataframe with unscaled predictors
# Calculate mean and sd to use later for burn severity map predictions
d_unscaled <- get_dat(input_data,
                      response = "prop_meanbolescorch", main_predictor = "rdnbr_w_offset",
                      other_predictors = other_predictors, 
                      scale_other_predictors = rep(FALSE, length(other_predictors))) 

d_mean <- d_unscaled %>% summarise(across(all_of(tolower(c("x", other_predictors))), mean))
d_sd <- d_unscaled %>% summarise(across(all_of(tolower(c("x", other_predictors))), sd))


# Fit model
m <- fit_model(d, model = sm_zoib, re = TRUE,
               predictors = c("xscaled", paste("xscaled", tolower(other_predictors), sep = "*")),
               iter = iter, chains = chains)

# Add predictor means and SDs to list
m$data_mean = d_mean
m$data_sd = d_sd

# Add number of post-warmup draws
m$ndraws = dim( as.matrix(m$model) )[1]

# Write to file
saveRDS(m, here("m_west_bchar_circ.rds"))



# Bole char (height) ----------

# Format data
other_predictors <- c("LF_psr")
scale_other_predictors <- c(TRUE)

# Get dataframe with scaled predictors for model-fitting
d <- get_dat(input_data,
             response = "prop_meancharht", main_predictor = "rdnbr_w_offset",
             other_predictors = other_predictors, 
             scale_other_predictors = scale_other_predictors)

# Get dataframe with unscaled predictors
# Calculate mean and sd to use later for burn severity map predictions
d_unscaled <- get_dat(input_data,
                      response = "prop_meancharht", main_predictor = "rdnbr_w_offset",
                      other_predictors = other_predictors, 
                      scale_other_predictors = rep(FALSE, length(other_predictors))) 

d_mean <- d_unscaled %>% summarise(across(all_of(tolower(c("x", other_predictors))), mean))
d_sd <- d_unscaled %>% summarise(across(all_of(tolower(c("x", other_predictors))), sd))


# Fit model
m <- fit_model(d, model = sm_zoib, re = TRUE,
               predictors = c("xscaled", paste("xscaled", tolower(other_predictors), sep = "*")),
               iter = iter, chains = chains)

# Add predictor means and SDs to list
m$data_mean = d_mean
m$data_sd = d_sd

# Add number of post-warmup draws
m$ndraws = dim( as.matrix(m$model) )[1]

# Write to file
saveRDS(m, here("m_west_bchar_ht.rds"))


# Surface char ----------

# Format data
other_predictors <- c("postfire2", "TT_cwd", "LF_slp", "Latitude_WGS84")
scale_other_predictors <- c(FALSE, TRUE, TRUE, TRUE)

# Get dataframe with scaled predictors for model-fitting
d <- get_dat(input_data,
             response = "prop_surfcharmean", main_predictor = "rdnbr_w_offset",
             other_predictors = other_predictors, 
             scale_other_predictors = scale_other_predictors)

# Get dataframe with unscaled predictors
# Calculate mean and sd to use later for burn severity map predictions
d_unscaled <- get_dat(input_data,
                      response = "prop_surfcharmean", main_predictor = "rdnbr_w_offset",
                      other_predictors = other_predictors, 
                      scale_other_predictors = rep(FALSE, length(other_predictors))) 

d_mean <- d_unscaled %>% summarise(across(all_of(tolower(c("x", other_predictors))), mean))
d_sd <- d_unscaled %>% summarise(across(all_of(tolower(c("x", other_predictors))), sd))


# Fit model
m <- fit_model(d, model = sm_zoib, re = TRUE,
               predictors = c("xscaled", paste("xscaled", tolower(other_predictors), sep = "*")),
               iter = iter, chains = chains)

# Add predictor means and SDs to list
m$data_mean = d_mean
m$data_sd = d_sd

# Add number of post-warmup draws
m$ndraws = dim( as.matrix(m$model) )[1]

# Write to file
saveRDS(m, here("m_west_schar.rds"))

