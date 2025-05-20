# Aim of the script: climatic data extraction


# 0. Packages and data loading --------------------------------------------
library(here)
library(tidyverse)
library(sf)
library(terra)
library(leaflet)
library(lme4)
library(nlme)
library(viridis)
library(brms)

# read fitted TPCs' thermal traits data set
therm_traits <- readRDS(here("data/data_sink/therm_traits_intrate.rds")) |> 
  mutate(across(c(reference, species, order, family, feed_guild,source_dataset),
                ~as_factor(.x))) |> 
  filter(ctmin > -5, #unrealistic fits
         ctmax <55) |> #unrealistic fits  
  mutate(therm_range = topt-ctmin) |> 
  group_by(reference, lon, lat)
## 1. climatic data extraction-----------------------------------------

### a) Present ----
#### 1. tmin ----

# choose data points
list_locations <- therm_traits |> 
  group_by(reference, species) |> 
  mutate(id_location = cur_group_id()) |> 
  ungroup() |> 
  select(id_location, lon, lat)

# convert to sf
points_intrates <-st_as_sf(x = list_locations,
                           coords = c("lon", "lat")) |> 
  st_set_crs(4326)

# convert to vector
points_vect_intrates <- vect(points_intrates)

#WorldClim Tmin
present_tmin_climate_wc <-  geodata::worldclim_global(var = "tmin",
                                                      res = 2.5,
                                                      path = tempdir())
## extract to spatial coords
points_tmin_present <- terra::extract(present_tmin_climate_wc,
                                      points_vect_intrates, 
                                      id = points_vect_intrates$id_location) |> 
  as_tibble() |> 
  pivot_longer(cols = -ID,
               names_to = "month",
               values_to = "tmin_avg") |> 
  mutate(month = str_sub(month, -2),
         month = as_factor(month))

#### 2. tmax ----

#WorldClim Tmax
present_tmax_climate_wc <-  geodata::worldclim_global(var = "tmax",
                                                      res = 2.5,
                                                      path = tempdir())
## extract to spatial coords
points_tmax_present <- terra::extract(present_tmax_climate_wc,
                                      points_vect_intrates, 
                                      id = points_vect_intrates$id_location) |> 
  as_tibble() |> 
  pivot_longer(cols = -ID,
               names_to = "month",
               values_to = "tmax_avg") |> 
  mutate(month = str_sub(month, -2),
         month = as_factor(month))

#### 3. tavg ----

points_tavg_present <- inner_join(points_tmin_present,
                                  points_tmax_present) |> 
  mutate(tavg = map2_dbl(.x = tmin_avg,
                         .y = tmax_avg,
                         .f = ~mean(c(.x, .y))),
         model = as_factor("WorldClim_1970-2000_v21")) |> 
  rename(tmin = tmin_avg,
         tmax = tmax_avg) |> 
  pivot_longer(cols = 3:5, 
               names_to = "var",
               values_to = "temp_value") |> 
  mutate(time_scenario = as_factor("Present")) |>
  filter(var == "tavg")

### b) Future ----
#WorldClim's Future data set (CMIP6) 
points_tavg_future <- readRDS(here("data/data_source/cmip6_tavg_2041-2060_ssp585_res25.rds")) |> 
  pivot_longer(cols = c("tmin", "tmax", "tavg"), 
               names_to = "var",
               values_to = "temp_value") |> 
  mutate(time_scenario = as_factor("Future"))


### c) Joined (Pres-Future) ----
points_tavg_worldclim <- points_tavg_present |> 
  bind_rows(points_tavg_future)
  
points_tavg_warmonth <- points_tavg_present |> 
  group_by(ID, time_scenario, model) |> 
  slice_max(temp_value) |> 
  ungroup() |> 
  group_by(ID, time_scenario) |> 
  summarise(temp_warmest_month = mean(temp_value))

points_tavg_warm_four <- points_tavg_worldclim |> 
  group_by(ID, time_scenario, model) |> 
  slice_max(temp_value, n = 4) |> 
  ungroup() |> 
  group_by(ID, time_scenario) |> 
  summarise(temp_warmest_four = mean(temp_value))

## 2. Modelling ----
### a) Warming Tolerance ----
#### a.1. data wrangling ----
##### as defined by Deutsch and what Rezende & Bozinovic (2019) call "Thermal Safety Margin": CTmax - T_exp_warmonth
points_tavg_warmonth <- points_tavg_present |> 
  group_by(ID, time_scenario, model) |> 
  slice_max(temp_value) |> 
  ungroup() |> 
  group_by(ID, time_scenario) |> 
  summarise(temp_warmest_month = mean(temp_value))

##obtain warming tolerances
wt_traits_nogauss <- points_tavg_warmonth |> 
  rename(id_location = ID) |> 
  inner_join(list_locations) |> 
  inner_join(therm_traits_nogaussians) |> 
  ungroup() |> 
  mutate(wt_match = map2_dbl(.x = ctmax,
                              .y = temp_warmest_month,
                              .f = ~abs(.x-.y)))
## set priors
prior_wt_lat <- c(prior(normal(0, 1), class = Intercept), #WT
                   prior(normal(0.1, 0.1), class = b),
                   prior(cauchy(0, 0.5), class = sd))

## warmest month
set.seed(2023)

wt_lat_bayes <- brm(wt_match ~ abs(lat) + (1|reference), #if we account for weights (se(vi)), Rhat = 3.00
                     data = wt_traits_nogauss,
                     prior = prior_wt_lat,
                     iter = 5000,
                     control = list(adapt_delta = 0.91))
summary(wt_lat_bayes)

## warmest quarter
wt_traits_four <- points_tavg_warm_four |> 
  rename(id_location = ID) |> 
  inner_join(list_locations) |> 
  inner_join(therm_traits_nogaussians) |> 
  ungroup() |> 
  mutate(wt_match_quarter = map2_dbl(.x = topt,
                                      .y = temp_warmest_four,
                                      .f = ~abs(.x-.y))) |> 
  select(id_location, time_scenario, lon, lat, reference, species, wt_match_quarter)

set.seed(2023)

wt_lat_bayes_quarter <- brm(wt_match_quarter ~ abs(lat) + (1|reference), #if we account for weights (se(vi)), Rhat = 3.00
                    data = wt_traits_four,
                    prior = prior_wt_lat,
                    iter = 5000,
                    control = list(adapt_delta = 0.91))
summary(wt_lat_bayes_quarter)

### b) Thermal Safety Margins ----
#### b.1. data wrangling ----
##### as defined by Deutsch et al. (2008): Topt - T_hab
####  T_hab either as mean annual temperature (Deutsch2008), mean of warmest quarter
####  and mean of warmest month (Rezende 2019)



tsm_traits <- points_tavg_warmonth |> 
  rename(id_location = ID) |> 
  inner_join(list_locations) |> 
  inner_join(therm_traits) |> 
  ungroup() |> 
  mutate(tsm_match = map2_dbl(.x = ctmax,
                                      .y = temp_warmest_month,
                                      .f = ~.x-.y)) |> 
  select(id_location, time_scenario, lon, lat, reference, species, tsm_match)

tsm_traits_four <- points_tavg_warm_four |> 
  rename(id_location = ID) |> 
  inner_join(list_locations) |> 
  inner_join(therm_traits) |> 
  ungroup() |> 
  mutate(tsm_match_quarter = map2_dbl(.x = ctmax,
                                     .y = temp_warmest_four,
                                     .f = ~.x-.y)) |> 
  select(id_location, time_scenario, lon, lat, reference, species, tsm_match_quarter)

tsm_traits_both <- inner_join(tsm_traits, tsm_traits_four)

tsm_traits_present <- tsm_traits_both |> 
  filter(time_scenario == "Present")


tsm_lat_lm <- lm(tsm_match ~ abs(lat),
                 tsm_traits_present)
performance::check_distribution(tsm_traits_present$tsm_match)
performance::check_model(tsm_lat_lm)

tsm_lat_lm_quarter <- lm(tsm_match_quarter ~ abs(lat),
                         tsm_traits_present)
## Bayesian for TSM warmest

library(brms)
prior_tsm_lat <- c(prior(normal(0, 1), class = Intercept), #TSM
                  prior(normal(0.1, 0.1), class = b),
                  prior(cauchy(0, 0.5), class = sd))


## random-intercept model
set.seed(2023)

tsm_lat_bayes <- brm(tsm_match ~ abs(lat) + (1|reference), #if we account for weights (se(vi)), Rhat = 3.00
                    data = tsm_traits_present,
                    prior = prior_tsm_lat,
                    iter = 5000,
                    control = list(adapt_delta = 0.91))
summary(tsm_lat_bayes)
plot(tsm_lat_bayes)
pp_check(tsm_lat_bayes)
ranef(tsm_lat_bayes)

post_samples <- posterior_samples(tsm_lat_bayes, c("^b", "^sd"))
names(post_samples) <- c("intercept","tsm_lat_est", "tau")


ggplot(aes(x = tsm_lat_est), data = post_samples) +
  geom_density(fill = "lightblue",                # set the color
               color = "lightblue", alpha = 0.7) +  
  geom_point(y = 0,                               # add point at mean
             x = mean(post_samples$tsm_lat_est)) +
  labs(x = expression(italic(TSM)[warmest]),
       y = element_blank()) +
  theme_minimal()


ggplot(aes(x = tau), data = post_samples) +
  geom_density(fill = "lightgreen",               # set the color
               color = "lightgreen", alpha = 0.7) +  
  geom_point(y = 0, 
             x = mean(post_samples$tau)) +        # add point at mean
  labs(x = expression(tau),
       y = element_blank()) +
  theme_minimal()

## bayesian for tsm quarter
## random-intercept model
set.seed(2023)

tsm_lat_bayes_quarter <- brm(tsm_match_quarter ~ abs(lat) + (1|reference), #if we account for weights (se(vi)), Rhat = 3.00
                             data = tsm_traits_present,
                             prior = prior_tsm_lat,
                             iter = 5000,
                             control = list(adapt_delta = 0.91))
summary(tsm_lat_bayes_quarter)
plot(tsm_lat_bayes_quarter)
pp_check(tsm_lat_bayes_quarter)
ranef(tsm_lat_bayes_quarter)

post_samples <- posterior_samples(tsm_lat_bayes_quarter, c("^b", "^sd"))
names(post_samples) <- c("intercept","tsm_lat_est", "tau")


ggplot(aes(x = tsm_lat_est), data = post_samples) +
  geom_density(fill = "lightblue",                # set the color
               color = "lightblue", alpha = 0.7) +  
  geom_point(y = 0,                               # add point at mean
             x = mean(post_samples$tsm_lat_est)) +
  labs(x = expression(italic(TSM)[warmest]),
       y = element_blank()) +
  theme_minimal()


ggplot(aes(x = tau), data = post_samples) +
  geom_density(fill = "lightgreen",               # set the color
               color = "lightgreen", alpha = 0.7) +  
  geom_point(y = 0, 
             x = mean(post_samples$tau)) +        # add point at mean
  labs(x = expression(tau),
       y = element_blank()) +
  theme_minimal()





## 3. Visualizations ----
#### 3.1. tsm ~ latitude ----
summary(tsm_lat_nogauss_lmer)

sim_and_plot_linears(model_object = tsm_lat_nogauss_lmer,
                     var_x = abs(tsm_traits_present_nogauss$lat),
                     var_y = tsm_traits_present_nogauss$tsm_match,
                     n_sims = 1000,
                     your_title = "Thermal Safety Margins ~ Latitude",
                     your_subtitle = NULL,
                     lab_y = expression(CT[max]~-~T[warmest~month]~(ºC)),
                     lab_x = "Absolute Latitude (º)")
ggsave(filename = here("data/data_sink/figs/tsm_lat_nogauss_lmer.png"),
       width = 16, height = 16, units = "cm")

## tsm warmest month
tsm_lat <-ggplot(tsm_traits_present, aes(x = lat, y = tsm_match))+
  geom_point(color = "#355070", size = 1.5, alpha = 0.56)+
  geom_smooth(color = "#b56576", fill = "#eaac8b")+
  labs(x = "Latitude(º)",
       y = expression(T[opt]~-~T[warmest~month]~(ºC)),
       title = "Thermal Safety Margins")+
  ggthemes::theme_few()+
  scale_x_continuous(breaks = seq(-50, 50, by = 10))+
  theme(plot.title = element_text(face = "bold"))


ggsave(filename = here("data/data_sink/figs/tsm_lat_loess.png"),
       width = 16, height = 16, units = "cm")
ggsave(filename = here("data/data_sink/figs/tsm_lat_loess.svg"),
       height = 3, width = 3)

## tsm warmest four months
ggplot(tsm_traits_present, aes(x = lat, y = tsm_match_quarter))+
  geom_point(color = "#b56576", size = 3, alpha = 0.6)+
  geom_smooth(color = "#355070", fill = "#eaac8b")+
  labs(x = "Latitude(º)",
       y = expression(T[opt]~-~T[warmest~month]~(ºC)),
  )+
  ggthemes::theme_clean()+
  scale_x_continuous(breaks = seq(-50, 50, by = 10))+
  coord_flip()

ggsave(filename = here("data/data_sink/figs/tsm_lat_loess_four.png"),
       width = 16, height = 16, units = "cm")
ggsave(filename = here("data/data_sink/figs/tsm_lat_loess_four.svg"),
       height = 3, width = 3)

#### 3.2. wt ~ latitude ----
## wt warmest month
wt_lat_loess <- ggplot(wt_traits_nogauss, aes(x = lat, y = wt_match))+
  geom_point(color = "#355070", size = 1.5, alpha = 0.56)+
  geom_smooth(color = "#b56576", fill = "#eaac8b")+
  labs(x = "Latitude(º)",
       y = expression(CT[max]~-~T[warmest~month]~(ºC)),
       title = "Warming Tolerances"
  )+
  ggthemes::theme_few()+
  scale_x_continuous(breaks = seq(-50, 50, by = 10))+
  theme(plot.title = element_text(face = "bold"))

ggsave(filename = here("data/data_sink/figs/wt_lat_loess.png"),
       width = 16, height = 16, units = "cm")
ggsave(filename = here("data/data_sink/figs/wt_lat_loess.svg"),
       height = 3, width = 3)


#### a.3. tsm ~ topt ----

sim_and_plot_linears(model_object = tsm_topt_nogauss_lmer,
                     var_x = abs(tsm_traits_present_nogauss$topt),
                     var_y = tsm_traits_present_nogauss$tsm_match,
                     n_sims = 1000,
                     your_title = "Thermal Safety Margins ~ T optimum",
                     your_subtitle = NULL,
                     lab_y = expression(CT[max]~-~Temp[warmest~month]~(ºC)),
                     lab_x = expression(T[opt]~(ºC)))
ggsave(filename = here("data/data_sink/figs/tsm_topt_nogauss_lmer.png"),
       width = 16, height = 16, units = "cm")

#### a.4. tsm ~ rmax ----
sim_and_plot_linears(model_object = tsm_rmax_nogauss_lmer,
                     var_x = abs(tsm_traits_present_nogauss$rmax),
                     var_y = tsm_traits_present_nogauss$tsm_match,
                     n_sims = 1000,
                     your_title = "Thermal Safety Margins ~ Maximum rate",
                     your_subtitle = NULL,
                     lab_y = expression(CT[max]~-~Temp[warmest~month]~(ºC)),
                     lab_x = expression(log(italic(r)[max])))
ggsave(filename = here("data/data_sink/figs/tsm_rmax_nogauss_lmer.png"),
       width = 16, height = 16, units = "cm")

#### a.5. wt + tsm abslat ----
wt_traits_comb <- wt_traits_nogauss |> 
  select(id_location, lon, lat, reference, species, wt_match)
tsm_traits_comb <- tsm_traits_present |> 
  select(id_location, lon, lat, reference, species, tsm_match)


therm_match_join <- left_join(tsm_traits_comb, wt_traits_comb) |> 
  pivot_longer(cols = c("tsm_match", "wt_match"),
               names_to = "matching_index",
               values_to = "therm_match")

match_lat <-ggplot(data = therm_match_join, aes(x = abs(lat), y = therm_match))+
  geom_point(size = 1.5, alpha = 0.56, aes(color = matching_index))+
  geom_smooth(method = "lm", aes(color = matching_index,
                                 fill = matching_index))+
  labs(x = "Absolute sLatitude(º)",
       y = "Temperature(ºC)",
       title = "Thermal Matching",
       color = NULL,
       fill = NULL)+
  ggthemes::theme_few()+
  scale_color_discrete(type =  c("#D19187", "#8D516D"),
                       labels = c("Thermal Safety Margins",
                                  "Warming Tolerances"))+
  scale_fill_discrete(type = c("#D19187", "#8D516D"),
                      labels = c("Thermal Safety Margins",
                                 "Warming Tolerances"))+
  
  scale_x_continuous(breaks = seq(-50, 50, by = 10))+
  theme(plot.title = element_text(face = "bold"),
        legend.position = "top")
match_lat



### a.6. Figure 1 ----

figure1_plotgrid <- cowplot::plot_grid(ctmin_to_lat_ref_plot,
                                       ctmax_to_lat_ref_plot,
                                       breadth_lat_plot,
                                       match_lat,
                                       tsm_lat,
                                       wt_lat_loess,
                                       nrow = 2,
                                       labels = "AUTO")
ggsave(figure2_plotgrid,
       filename = here("data/data_sink/figs/traits_lat/figure1_plotgrid.png"),
       width = 30, height = 20, units = "cm")
ggsave()