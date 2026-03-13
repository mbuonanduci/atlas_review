# Use fitted models to generate atlas
# Code developed by M. Buonanduci

# Load packages and functions
library(here)
library(tidyverse)
library(posterior)
library(terra)
library(tidyterra)
library(sf)
library(parallel)
library(foreach)
library(doParallel)

# Load custom functions
source(here("00_zoib_functions.R"))


# Specify response variables --------
# Vector of: BA, stems, CC, needle, bchar_circ, bchar_ht, schar
resp_vars <- c("BA", "stems", "CC", "needle", "bchar_circ", "bchar_ht", "schar")


# Specify settings for parallel processing ---------

# Processing "chunk" size
csize <- 1000

# Number of cores
cl = makeCluster(min(10, detectCores()))
registerDoParallel(cl)


# Read in model fits --------
mods <- list()
for(r in 1:length(resp_vars)){
  mods[[r]] <- here(paste0("m_west_", resp_vars[r], ".rds")) %>% readRDS()
}


# Read in fires --------

# Read in fires
fire_perims <- here("data", "fire_perims_atlas.gpkg") %>% st_read()

# List of fires 
fire_IDs <- unique(fire_perims$fireid)

# List of years
years <- unique(fire_perims$fire_year) %>% as.numeric()


# Read in climate, forest structure, and topography covariates ----------

# Climate water deficit raster
# CRS is WGS84, ~220m resolution
cwd_rast <- here("source_data", "TopoTerra", "TopoTerra_1985-2015.tif") %>% rast()
cwd_rast <- cwd_rast$def

# Historical percent stand-replacing raster
# This raster has been reclassified from categorical to numeric, reflecting 'percent stand-replacing fire'
# CRS is EPSG:5070, 30m resolution
LF_psr_rast <- here("source_data", "LANDFIRE", "LF2016_PFS_200_CONUS", "Tif", "LC16_PFS_200_West_sr.tif") %>% rast()

# Slope raster
# CRS is EPSG:5070, 30m resolution
LF_slp_rast <- here("source_data", "LANDFIRE", "LF2020_SlpD_220_CONUS", "Tif", "LC20_SlpD_220_West.tif") %>% rast()
# Set NA flag
NAflag(LF_slp_rast) <- -9999


# Start processing loop ----------

# Loop through fire years ---------
for(y in 1:length(years)){
  
  # Read in pre-fire tree cover raster -------
  tcov_rast <- here("source_data", "RAP", paste0("RAP_TRE_", years[y]-1, ".tif")) %>% rast()
  
  # Loop through fire events ----------
  # Fires occurring in fire year y 
  fire_perims_y <- fire_perims %>%
    filter(fire_year == years[y])
  fire_IDs_y <- fire_perims_y$fireid
  
  for(f in 1:length(fire_IDs_y)){
    
    # Read in satellite index raster --------
    # CRS is EPSG:5070, 30m resolution
    # Mask using fire perimeter; only want to make predictions for cells within fire event
    satindex_rast <- here("data", "RdNBR_for_atlas", paste0( fire_IDs_y[f], "_rdnbr_w_offset.tif") ) %>% rast() %>%
      mask(fire_perims_y[f,])
    
    # Get raster cell coordinates, convert to sf object
    crds_sf <- st_as_sf(tibble(lon = crds(satindex_rast, na.rm=TRUE)[,1],
                               lat = crds(satindex_rast, na.rm=TRUE)[,2]),
                        coords = c("lon", "lat"), crs = st_crs("EPSG:5070"))
    
    # Get latitude of cells in WGS84
    crds_latitude <- crds_sf %>%
      st_transform(crs = st_crs("EPSG:4326")) %>%
      mutate(lat = st_coordinates(.)[,2])
    
    # Get flag for forest vs. non-forest areas 
    # Use pre-fire tree cover threshold of >=10% tree cover to designate "forest"
    # Extract pre-fire tree cover
    # Raster CRS and scale are the same
    tcov <- terra::extract(tcov_rast, vect(crds_sf), method = "simple") %>%
      mutate(forest_flag = ifelse(TRE<10, NA, 1))
    
    # Loop through burn severity models --------
    # Make predictions; store in list
    burnsev_pred <- list()
    for(r in 1:length(resp_vars)){
      
      # Get model
      m <- mods[[r]]
      
      # Get response variable
      resp <- resp_vars[r]
      
      # Create dataframe of predictors ----------
      # Get satellite index (main predictor, denoted 'x')
      # Set sampling interval = 0, indicating 1-year post-fire predictions)
      # Set group_id (fire id) = NA, as we don't include random effects
      # Add forest flag to carry through calculations
      predictors <- tibble(lon = crds(satindex_rast, na.rm=TRUE)[,1],
                           lat = crds(satindex_rast, na.rm=TRUE)[,2],
                           x = as.numeric(values(satindex_rast, na.rm = TRUE))) %>%
        mutate(xscaled = (x - m$data_mean$x) / (2 * m$data_sd$x)) %>%
        mutate(xraw_mean = m$data_mean$x,
               xraw_sd = m$data_sd$x) %>%
        mutate(postfire2 = 0) %>% # 1-year post-fire predictions
        mutate(group_id = NA) %>% # no fire event-level random effects
        mutate(forest_flag = tcov$forest_flag) # add forest flag
      
      # Add pre-fire tree cover covariate
      if(resp %in% c("BA", "CC")){
        tcov <- tcov %>%
          mutate(rap_tcov = (TRE - m$data_mean$rap_tcov) / (2 * m$data_sd$rap_tcov))
        predictors <- predictors %>%
          mutate(rap_tcov = tcov$rap_tcov)
      }
      
      # Add CWD covariate
      # Raster CRS and scale differ
      # Need to re-project points used for extraction
      if(resp %in% c("BA", "needle", "bchar_circ", "schar")){
        crds_sf_reproj <- st_transform(crds_sf, crs = st_crs(cwd_rast))
        cwd <- terra::extract(cwd_rast, vect(crds_sf_reproj), method = "simple") %>%
          mutate(tt_cwd = (def - m$data_mean$tt_cwd) / (2 * m$data_sd$tt_cwd))
        predictors <- predictors %>%
          mutate(tt_cwd = cwd$tt_cwd)
      }
      
      # Add historical percent stand-replacing covariate
      # Raster CRS and scale are the same
      if(resp == "bchar_ht"){
        LF_psr <- terra::extract(LF_psr_rast, vect(crds_sf), method = "simple") %>%
          mutate(lf_psr = (Value - m$data_mean$lf_psr) / (2 * m$data_sd$lf_psr))
        predictors <- predictors %>%
          mutate(lf_psr = LF_psr$lf_psr)
      }
      
      # Add slope covariate 
      # Raster CRS and scale are the same
      if(resp %in% c("needle", "schar")){
        LF_slp <- terra::extract(LF_slp_rast, vect(crds_sf), method = "simple") %>%
          mutate(lf_slp = (LC20_SlpD_220 - m$data_mean$lf_slp) / (2 * m$data_sd$lf_slp))
        predictors <- predictors %>%
          mutate(lf_slp = LF_slp$lf_slp)
      }
      
      # Add latitude covariate
      if(resp %in% c("schar")){
        crds_latitude <- crds_latitude %>%
          mutate(latitude_wgs84 = (lat - m$data_mean$latitude_wgs84) / (2 * m$data_sd$latitude_wgs84))
        predictors <- predictors %>%
          mutate(latitude_wgs84 = crds_latitude$latitude_wgs84)
      }
      
      # Drop potential NAs for predictors
      if("rap_tcov" %in% names(predictors)) predictors <- predictors %>% drop_na(rap_tcov)
      if("tt_cwd" %in% names(predictors)) predictors <- predictors %>% drop_na(tt_cwd)
      if("lf_psr" %in% names(predictors)) predictors <- predictors %>% drop_na(lf_psr)
      if("lf_slp" %in% names(predictors)) predictors <- predictors %>% drop_na(lf_slp)
      
      # Use fitted model to make predictions --------
      # Run in parallel
      
      # Create processing "chunks"
      nchunk <- ceiling(nrow(predictors) / csize)
      predictors <- predictors %>% mutate(chunk = rep(c(1:nchunk), length.out = nrow(predictors)))
      
      # Create output layer names
      resp_lwr = paste0(resp, "_lwr")
      resp_upr = paste0(resp, "_upr")
      resp_gt75 = paste0(resp, "_gt75")
      
      # Start parallel operation
      burnsev_pred_r <- foreach(i = 1:nchunk, .packages = c('tidyverse', 'posterior')) %dopar% {
        
        # Get predictors
        predictors_i <- filter(predictors, chunk == i) %>% mutate(y = 1)
        
        # Calculate expected value
        exp_i <- make_predictions_exp(d = predictors_i, f = m$f, model = m$model, 
                                      use_new_data = FALSE, re = FALSE, est_range = FALSE)
        
        # Calculate posterior predictive distribution
        yrep_i <- make_predictions_yrep(d = predictors_i, f = m$f, model = m$model, 
                                        use_new_data = FALSE, re = FALSE, est_range = FALSE)
        
        tibble(lon = predictors_i$lon,
               lat = predictors_i$lat,
               forest = predictors_i$forest_flag,
               rdnbr = predictors_i$x,
               !!sym(resp) := median(exp_i$est),
               !!sym(resp_lwr) := t(quantile(yrep_i$est, 0.05))[,1],
               !!sym(resp_upr) := t(quantile(yrep_i$est, 0.95))[,1],
               !!sym(resp_gt75) := sum(yrep_i$est > 0.75) / (m$ndraws))
      }
      
      burnsev_pred[[r]] <- bind_rows(burnsev_pred_r)
      
    }
    
    # Join together all results
    burnsev_pred_join <- tibble(lon = crds(satindex_rast, na.rm=TRUE)[,1],
                                lat = crds(satindex_rast, na.rm=TRUE)[,2])
    for(r in 1:length(resp_vars)){
      burnsev_pred_join <- burnsev_pred_join %>%
        left_join(burnsev_pred[[r]])
    }
    
    # Create raster object from outputs
    burnsev_pred_rast <- burnsev_pred_join %>%
      as_spatraster(xycols = 1:2, crs = crs(satindex_rast))
    
    # Write to file
    writeRaster(burnsev_pred_rast, 
                here("data", "atlas", paste0(fire_IDs_y[f], ".tif")), 
                overwrite=TRUE)
  }
}


