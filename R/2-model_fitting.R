
# Aim of the script: fit tpcs to each study and apply model selection


# 0. Packages and data loading --------------------------------------------
library(rTPC)
library(nls.multstart)
library(here)
library(readxl)
library(tidyverse)
source(here("analyses/S1-functions.R"))

# setwd() <- choose your directory's' highest level folder
## load data
int_rate_pest <- read_delim("data/data_source/int_rate_pest_dataset.csv") |> 
  as_tibble() |>
  mutate(source_dataset = as_factor(source_dataset)) |> 
  mutate(across(c(reference, species, order, family, feed_guild, source_dataset),
                ~as_factor(.x)))
# 1. Model fitting --------------------------------------------------------

## we'll use our function `fit_tpc()` and store its output (a tibble) as a new nested column in `int_rate_pest`

## this function fits nonlinear-regression models (see packages `rTPC` `nls.multstart`)
fit_tpcs <- int_rate_pest |>
  group_by(reference, species) |>
  summarise(nest(fit_tpc(temperature = temperature, # <- nest for processing
                         rate = int_rate),
                         .key = "tpcs_params"))

#unnest, disentangle for obtaining the fitted models' table
fit_tpcs_unnest <- fit_tpcs |> 
  unnest(cols = c(tpcs_params))
saveRDS(fit_tpcs_unnest, file = here("data/data_sink/fit_tpcs_unnest.rds"))
  
# 2. Model selection & Plots  --------------------------------------------------------

fitted_tpcs_joined <- readRDS(here("data/data_sink/fit_tpcs_unnest.rds")) |> 
  mutate(reference = as_factor(reference),
         species = as_factor(species)) |> 
  ungroup() |>
  inner_join(int_rate_pest,
             by = c("reference", "species")) |> 
  mutate(int_rate = map_dbl(.x = int_rate,
                            .f = ~if_else(.x < 0, # <- to help model fitting, negative values are assumed to be zero.
                                          0,
                                          .x)))

## plot all fitted TPCs and their uncertainties

for(study_i in distinct_studies){
  print(study_i)
  
  species_distinct_study <- fitted_tpcs_joined |> 
    filter(reference == study_i) |>
    distinct(species) |> 
    pull(species)
  for(species_i in species_distinct_study) {
    study_species_i <- fitted_tpcs_joined |> 
      filter(species == species_i &
               reference == study_i)
    study_sp_models <- study_species_i |> 
      distinct(model_name) |> 
      pull(model_name)
    
    for(model_sp_i in study_sp_models) {
      study_sp_model_i <- study_species_i |> 
        filter(model_name == model_sp_i)
      plot_tpc_species <- plot_tpcs(temperature = study_sp_model_i |> distinct(temperature) |> pull(temperature),
                                    rate = study_sp_model_i |> group_by(temperature) |> summarise(int_rate = mean(int_rate)) |> pull(int_rate),
                                    fitted_parameters = study_sp_model_i |> select(param_name, start_vals, param_est, param_se, 
                                                                                   model_name, model_AIC, model_BIC, model_fit),
                                    study = study_i,
                                    species = species_i,
                                    aic_model = study_sp_model_i |> distinct(model_AIC) |> pull(model_AIC),
                                    facets = FALSE,
                                    uncertainty = TRUE,
                                    uncertainty_sims = 100)  
      ggsave(plot_tpc_species,
             filename =  here(paste0("data/data_sink/figs/tpcs/", study_i, str_sub(word(species_i, 1), 1, 1), str_sub(word(species_i, 2), 1, 1),"_", model_sp_i, ".png")),
             device =  "png", 
             width = 15,
             height = 15,
             units = "cm")
    }
  }
  
}
  


# 3. Exploration of fitted_tpcs --------------------------------------------------------

fitted_tpcs_intrate <- readRDS(here("data/data_sink/tpcs_fit.rds"))

## succesfully fitted
fitted_tpcs_intrate |> 
  group_by(reference, species) |> 
  count()

## which models perform better?
distinct_models <- fitted_tpcs_intrate |> 
  group_by(reference, species, model_name) |> 
  distinct(model_name) 
distinct_studies_models <- distinct_models |> 
  ungroup() |> 
  count(model_name)

top_scorers <- fitted_tpcs_intrate %>%
  group_by(reference, species) %>%
  arrange(model_AIC) %>%
  mutate(min_AIC_diff = model_AIC - lag(model_AIC),
         min_AIC_diff = ifelse(is.na(min_AIC_diff), 
                               Inf, 
                               min_AIC_diff)) %>%
  filter(min_AIC_diff > 2) %>%
  slice_head(n = 1) %>%
  select(model_name)

distinct_scorers <- top_scorers |> 
  ungroup() |> 
  count(model_name)

  
ggplot(data = distinct_scorers, aes(x = fct_reorder(model_name, n),
                                          y = n))+
  geom_bar(aes(fill = model_name), stat = "identity")+
  scale_fill_manual(values = MoMAColors::MoMAPalettes$Klein[[1]])+
  theme_bw()+
  coord_flip()+
  theme(legend.position = "none")+
  labs(x = "Selected model by lowest AIC",
       y = "number of records")
ggsave(here("data/data_sink/figs/alltpcs_by_aic.png"),
       height = 13,
       width = 13,
       units = "cm")


ggplot(data = distinct_studies_models, aes(x = fct_reorder(model_name, n),
                                           y = n))+
  geom_bar(aes(fill = model_name), stat = "identity")+
  scale_fill_manual(values = MoMAColors::MoMAPalettes$Klein[[1]])+
  theme_bw()+
  coord_flip()+
  theme(legend.position = "none")+
  labs(x = "Selected model by lowest AIC",
       y = "number of records")
ggsave(here("data/data_sink/figs/alltpcs_by_aic.png"),
       height = 13,
       width = 13,
       units = "cm")

# 4. Selected tpcs ----

## read csv with manual decisions on tpc selection

model_selection_tbl <- read_delim(here("data/data_source/model_selection_tpcs.csv")) |> 
  filter(decision %in% c("pass", "choose", "choose_condition"))


selected_tpcs_intrate <- fitted_tpcs_intrate |> 
  inner_join(model_selection_tbl, by = c("reference", "species", "model_name")) |> 
  ungroup() |> 
  select(-model_AIC.y) |> 
  rename(model_AIC = model_AIC.x)
  

##### a) n effect sizes, studies and species ----------------------------------------------------

# number of effect sizes
effect_sizes_id <- selected_tpcs_intrate |> 
  group_by(reference) |> 
  count(species) |> 
  select(reference, species)

n_ef_size <-  selected_tpcs_intrate |> 
  distinct(species) |> 
  nrow() 
# 113 effect sizes

# number of studies
n_studies <- selected_tpcs_intrate |> 
  ungroup() |> 
  distinct(reference) |> 
  nrow()
# 109 studies

# number of distinct species
n_species <- selected_tpcs_intrate |> 
  distinct(species) |> 
  nrow()
 #113 species 

##### b) therm_traits_intrate ----------------------------------------------------

# let's extract thermal traits from `nls_multstart` models using `rTPC`
therm_traits_intrate_raw <- selected_tpcs_intrate |> 
  group_by(reference, species) |> 
  slice(1) |> 
  select(reference, species, order, family, model_fit, decision, model_name) |> 
  mutate(model_fit = map(.x = model_fit, 
                         .f = ~model_fit[[1]]),
         thermal_traits = map(.x = model_fit,
                              .f = ~rTPC::calc_params(.x))
         )
  
therm_traits_intrate <- therm_traits_intrate_raw |> 
  unnest(thermal_traits) |> 
  inner_join(int_rate_pest) |> 
  select(-temperature, -int_rate ) |> 
  group_by(reference, species) |> 
  summarise(across(everything(),
                   ~ifelse(is.numeric(.x),
                           median(.x),
                           unique(.x)
                           )
                   )
            ) |> 
  ungroup() |> 
  rename(vi = int_rate_se)

## how many are suitable for strict meta-analysis
n_suit_meta <- therm_traits_intrate |> 
  filter(!is.na(vi)) |> 
  n_distinct()
# 60 effect sizes

##### c) export therm_traits_intrate ----------------------------------------------------
saveRDS(therm_traits_intrate,
        file = here("data/data_sink/therm_traits_intrate.rds"))

##### d) plot all tpcs ----------------------------------------------------
rtpc_formulas_2join <- rTPC_choice_tbl |> 
  rename(model_name = rTPC_models)
int_rate_pest_id <- int_rate_pest |> 
  select(reference, species, temperature, int_rate)

selected_tpcs_joined <- inner_join(selected_tpcs_intrate, int_rate_pest_id, 
                                   by = c("reference", "species")) |> 
  inner_join(rtpc_formulas_2join)
seq_temps <- seq(min(therm_traits$ctmin), max(therm_traits$ctmax), by = 0.01)
list_refs <- unique(therm_traits$reference)
ggplot2fill <- ggplot()+
  theme_few()+
  labs(title = "Fitted TPCs",
       x = "Temperature (ºC)",
       y = expression(intrinsic~rate~of~increase~(italic(r)[m])))

##plot TPCs
for(study in unique(selected_tpcs_joined$reference)){t
  reference_tbl <- selected_tpcs_joined |> 
    filter(reference == study)
  for(species_i in unique(reference_tbl$species)){
    reference_species_tbl <- reference_tbl |> 
      filter(species == species_i)
  predict2fill <- tibble(temperature = NULL,
                         rate = NULL)
  params_i <- reference_species_tbl |>
    distinct(param_est) |> 
    pull(param_est)
   formula_i <- reference_species_tbl |> 
     distinct(nls_formula_params) |> 
      pull(nls_formula_params)
    ##predict based on parameters
    explore_preds <- tibble(temperature = seq(-5,
                                              60,
                                              .01),
                            preds = NULL)
    fit_vals_tbl <- explore_preds |>
      dplyr::select(temperature) |>
      mutate(formula = formula_i) |>
      mutate(preds = purrr::map_dbl(.x = temperature,
                                    .f = reformulate(unique(formula_i)))) |>
      filter(preds >= 0) |>
      dplyr::select(-formula)
    predict2fill <- predict2fill |>
      bind_rows(fit_vals_tbl)
  ggplot_row <- ggplot()+
      geom_line(data = predict2fill,
                aes(x = temperature, 
                    y = preds),
                color = "darkcyan",
                linewidth = 0.25,
                alpha = .6)+
      theme_void()+
      theme(legend.position = "none")
  
  ggplot2fill <- ggplot2fill+ggplot_row[[2]]
  }
  }
ggplot2fill
ggsave(filename = here("data/data_sink/figs/all_tpcs_fitted.png"),
       width = 12, height = 12, units = "cm"
       )
 ## same but excluding gaussians
rtpc_formulas_2join_v2 <- rTPC_choice_tbl |> 
  rename(model_name = rTPC_models) |> 
  filter(model_name != "gaussian_1987")
selected_tpcs_joined_v2 <- inner_join(selected_tpcs_intrate, int_rate_pest, 
                                   by = c("reference", "species", "order", "family")) |> 
  inner_join(rtpc_formulas_2join_v2)
seq_temps <- seq(min(therm_traits$ctmin), max(therm_traits$ctmax), by = 0.01)
list_refs <- unique(therm_traits$reference)
ggplot2fill_v2 <- ggplot()+
  theme_few()+
  labs(title = "Fitted TPCs",
       x = "Temperature (ºC)",
       y = expression(intrinsic~rate~of~increase~(italic(r)[m])))
for(study in unique(selected_tpcs_joined_v2$reference)){t
  reference_tbl <- selected_tpcs_joined_v2 |> 
    filter(reference == study)
  for(species_i in unique(reference_tbl$species)){
    reference_species_tbl <- reference_tbl |> 
      filter(species == species_i)
    predict2fill <- tibble(temperature = NULL,
                           rate = NULL)
    params_i <- reference_species_tbl |>
      distinct(param_est) |> 
      pull(param_est)
    formula_i <- reference_species_tbl |> 
      distinct(nls_formula_params) |> 
      pull(nls_formula_params)
    ##predict based on parameters
    explore_preds <- tibble(temperature = seq(-5,
                                              60,
                                              .01),
                            preds = NULL)
    fit_vals_tbl <- explore_preds |>
      dplyr::select(temperature) |>
      mutate(formula = formula_i) |>
      mutate(preds = purrr::map_dbl(.x = temperature,
                                    .f = reformulate(unique(formula_i)))) |>
      filter(preds >= 0) |>
      dplyr::select(-formula)
    predict2fill <- predict2fill |>
      bind_rows(fit_vals_tbl)
    ggplot_row <- ggplot()+
      geom_line(data = predict2fill,
                aes(x = temperature, 
                    y = preds),
                color = "darkcyan",
                linewidth = 0.25,
                alpha = .6)+
      theme_void()+
      theme(legend.position = "none")
    
    ggplot2fill_v2 <- ggplot2fill_v2+ggplot_row[[2]]
  }
}
ggplot2fill_v2
ggsave(filename = here("data/data_sink/figs/all_tpcs_fitted_nogauss.png"),
       width = 12, height = 12, units = "cm"
)

## Arrhenius space:
##plot TPCs
for(study in unique(selected_tpcs_joined$reference)){
  reference_tbl <- selected_tpcs_joined |> 
    filter(reference == study)
  for(species_i in unique(reference_tbl$species)){
    reference_species_tbl <- reference_tbl |> 
      filter(species == species_i)
    predict2fill <- tibble(temperature = NULL,
                           rate = NULL)
    params_i <- reference_species_tbl |>
      distinct(param_est) |> 
      pull(param_est)
    formula_i <- reference_species_tbl |> 
      distinct(nls_formula_params) |> 
      pull(nls_formula_params)
    ##predict based on parameters
    explore_preds <- tibble(temperature = seq(-5,
                                              60,
                                              .01),
                            preds = NULL)
    fit_vals_tbl <- explore_preds |>
      dplyr::select(temperature) |>
      mutate(arrhenius_temp = 1/(8.617e-05*(temperature + 273.15))) |> 
      mutate(formula = formula_i) |>
      mutate(preds = purrr::map_dbl(.x = arrhenius_temp,
                                    .f = reformulate(unique(formula_i)))) |>
      dplyr::select(-formula) |> 
      mutate(preds = log(preds))
    predict2fill_arrhenius <- predict2fill |>
      bind_rows(fit_vals_tbl) 
    
    points_arrhenius <- reference_species_tbl |> 
      select(temperature, int_rate) |> 
      group_by(temperature) |> 
      summarise(int_rate = mean(int_rate)) |> 
      mutate(arrhenius_temp = 1/(8.617e-05*(temperature + 273.15)),
             log_int_rate = log(int_rate))
    
    
    ggplot_row <- ggplot()+
      
      geom_line(data = predict2fill_arrhenius,
                aes(x = arrhenius_temp, 
                    y = preds),
                color = "darkcyan",
                linewidth = 0.25,
                alpha = .6)+
      geom_point(data = points_arrhenius, aes(x = arrhenius_temp,
                                              y = log_int_rate))+
      theme_void()+
      theme(legend.position = "none")
    
    ggplot2fill <- ggplot2fill+ggplot_row[[2]]
  }
}
ggplot2fill
ggsave(filename = here("data/data_sink/figs/all_tpcs_fitted_arrhenius.png"),
       width = 12, height = 12, units = "cm"
)


