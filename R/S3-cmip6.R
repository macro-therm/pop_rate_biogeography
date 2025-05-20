#### script for assemblying the Future CMIP6 models (RCP 8.5, 2041-2060)
## data downloaded from https://geodata.ucdavis.edu/cmip6 for 2.5 resolution

library(tidyverse)
library(terra)


## a. cmip6_tmax ----
# cmip6_tmax
cmip6_models <-  c("ACCESS-CM2", "ACCESS-ESM1-5", "AWI-CM-1-1-MR", "BCC-CSM2-MR",
                   "CanESM5", "CanESM5-CanOE", "CMCC-ESM2", "CNRM-CM6-1",
                   "CNRM-CM6-1-HR", "CNRM-ESM2-1", "EC-Earth3-Veg", 
                   "EC-Earth3-Veg-LR","FIO-ESM-2-0", "GISS-E2-1-G",
                   "GISS-E2-1-H", "HadGEM3-GC31-LL", "INM-CM4-8", "INM-CM5-0",
                   "IPSL-CM6A-LR", "MIROC-ES2L", "MIROC6", "MPI-ESM1-2-HR",
                   "MPI-ESM1-2-LR", "MRI-ESM2-0", "UKESM1-0-LL")

cmip6_models_tmax_tbl <- tibble(ID = NULL,
                                month = NULL,
                                tmax = NULL,
                                model = NULL)

for(model_cmip6 in cmip6_models){
  
  print(paste0(100*which(cmip6_models == model_cmip6)/length(cmip6_models),"%"))
  model_i <- paste0("wc2.1_2.5m_tmax_", model_cmip6, "_ssp585_2041-2060.tif")
  raster_model_i <- terra::rast(paste0("D:/cmip6/", model_i))
  points_tmax_cmip6_i<- terra::extract(raster_model_i,
                                       points_vect_intrates, 
                                       id = points_vect_intrates$id_location) |> 
    as_tibble() |> 
    pivot_longer(cols = -ID,
                 names_to = "month",
                 values_to = "tmax") |> 
    mutate(month = as_factor(str_sub(month, -2)),
           model = as_factor(model_cmip6))
  cmip6_models_tmax_tbl <- bind_rows(cmip6_models_tmax_tbl, points_tmax_cmip6_i)
}

## b. cmip6_tmin ----
# cmip6_tmin
cmip6_models_tmin_tbl <- tibble(ID = NULL,
                                month = NULL,
                                tmin = NULL,
                                model = NULL)

for(model_cmip6 in cmip6_models){
  
  print(paste0(100*which(cmip6_models == model_cmip6)/length(cmip6_models),"%"))
  model_i <- paste0("wc2.1_2.5m_tmin_", model_cmip6, "_ssp585_2041-2060.tif")
  raster_model_i <- terra::rast(paste0("D:/cmip6/", model_i))
  points_tmin_cmip6_i<- terra::extract(raster_model_i,
                                       points_vect_intrates, 
                                       id = points_vect_intrates$id_location) |> 
    as_tibble() |> 
    pivot_longer(cols = -ID,
                 names_to = "month",
                 values_to = "tmin") |> 
    mutate(month = as_factor(str_sub(month, -2)),
           model = as_factor(model_cmip6))
  cmip6_models_tmin_tbl <- bind_rows(cmip6_models_tmin_tbl, points_tmin_cmip6_i)
}

## c. cmip6_tavg ----

cmip6_models_tavg <- inner_join(cmip6_models_tmin_tbl, cmip6_models_tmax_tbl) |> 
  mutate(tavg = map2_dbl(.x = tmin, 
                         .y = tmax,
                         .f = ~mean(c(.x, .y)))
         )
saveRDS(object = cmip6_models_tavg,
        file = here("data/data_source/cmip6_tavg_2041-2060_ssp585_res25.rds"))



