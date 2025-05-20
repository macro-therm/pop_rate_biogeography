### scripts with functions

## 1. Model fitting with `rTPC` and `nls.multstart` ----
library(rTPC)
library(nls.multstart)
library(broom)
library(pracma)

#first define model_names list
rTPC_models <- c(# "delong_2017",  <- no error generation
                 "deutsch_2008",
                #  "flinn_1991", <- delete it to consistently large std. errors
                 "gaussian_1987",
                 "pawar_2018", 
                # "ratkowsky_1983", <- delete it to unrealistic fitting
                 "rezende_2019",
                 "sharpeschoolfull_1981",
               #  "sharpeschoolhigh_1981", <- similar to pawar
               #  "spain_1982", <- bad fitting
                 "weibull_1995")

nls_ms_formula <- c(# "delong_2017(temperature, c, eb, ef, tm, ehc)", <- no error generation
                          "deutsch_2008(temperature, rmax, topt, ctmax, a)",
                        #  "flinn_1991(temperature, a, b, c)", <- delete it to consistently large std. errors
                          "gaussian_1987(temperature, rmax, topt, a)",
                          "pawar_2018(temperature, r_tref, e, eh, topt, tref = 25)",
                       #   "ratkowsky_1983(temperature, tmin, tmax, a, b)", <- delete it to unrealistic fitting
                          "rezende_2019(temperature, q10, a, b, c)",
                          "sharpeschoolfull_1981(temperature, r_tref, e, el, tl, eh, th, tref = 25)",
                      #   "sharpeschoolhigh_1981(temperature, r_tref, e, eh, th, tref = 25)", <- similar to pawar
                       #   "spain_1982(temperature, a, b, c, r0)", <- bad fitting
                          "weibull_1995(temperature, a, topt, b, c)")
nls_formula_params <- c(# "delong_2017(.x, params_i[1], params_i[2], params_i[3],params_i[4], params_i[5])", <- no error generation
                        "deutsch_2008(.x, params_i[1], params_i[2], params_i[3], params_i[4])",
                       # "flinn_1991(.x, params_i[1], params_i[2], params_i[3])", <- delete it to consistently large std. errors
                        "gaussian_1987(.x, params_i[1], params_i[2], params_i[3])",
                        "pawar_2018(.x, params_i[1], params_i[2], params_i[3], params_i[4], tref = 25)",
                       # "ratkowsky_1983(.x, params_i[1], params_i[2], params_i[3], params_i[4])", <- delete it to unrealistic fitting
                        "rezende_2019(.x, params_i[1], params_i[2], params_i[3], params_i[4])",
                        "sharpeschoolfull_1981(.x, params_i[1], params_i[2], params_i[3], params_i[4], params_i[5], params_i[6], tref = 25)",
                       # "sharpeschoolhigh_1981(.x, params_i[1], params_i[2], params_i[3], params_i[4], tref = 25)", <- similar to pawar
                       # "spain_1982(.x, params_i[1], params_i[2], params_i[3], params_i[4])", <-- bad convergence
                        "weibull_1995(.x, params_i[1], params_i[2], params_i[3],params_i[4])")
rTPC_n_params <- c(# 5, 
                   4, 
                   # 3, <- delete it to consistently large std. errors
                   3, 
                   4, 
                   # 4, <- delete it to unrealistic fitting 
                   4, 
                   6, 
                  # 4, 
                  # 4, 
                   4)
rTPC_choice_tbl <- tibble(rTPC_models,
                          nls_ms_formula,
                          rTPC_n_params,
                          nls_formula_params)


extract_param_names <- function(nls_object){
  parameter_est <- coef(nls_object)
  param_names <- names(parameter_est)
  return(param_names)
}



fit_tpc <- function(temperature, rate) {
  models_2fit <- rTPC_choice_tbl |> 
    filter(rTPC_n_params <= n_distinct(temperature)) |> 
    pull(rTPC_models)
  list_fit_models <- vector("list", length = length(models_2fit))
  list_param <- dplyr::tibble(param_name = NULL,
                              start_vals = NULL,
                              param_est = NULL,
                              param_se = NULL,
                              model_name = NULL,
                              model_AIC = NULL,
                              model_BIC = NULL,
                              model_fit = NULL)

    for (i in models_2fit) {
    message(paste0("fitting model ", i)) # to let people know that the function is working and R is not crashing
    
      possible_error <- tryCatch(expr = {start_vals_i <- rTPC::get_start_vals(x = temperature,
                                                                        y = rate,
                                                                        model_name = i)
      model_i <- rTPC_choice_tbl |>
        filter(rTPC_models == i)
      therm_perf_df <- dplyr::tibble(temperature, rate)
      fit_nls <- nls_multstart(formula =reformulate(response = "rate",
                                                    termlabels = unique(model_i$nls_ms_formula)),
                               data = therm_perf_df,
                               iter = 500, 
                               start_lower = start_vals_i-10,
                               start_upper = start_vals_i+10,
                               lower = get_lower_lims(therm_perf_df$temperature,
                                                      therm_perf_df$rate,
                                                      model_name = i),
                               upper = get_upper_lims(therm_perf_df$temperature,
                                                      therm_perf_df$rate,
                                                      model_name = i),
                               supp_errors = "Y")
      list_fit_models[[which(rTPC_models == i)]] <- fit_nls
      sum_fit_nls <- summary(fit_nls)
      list_param_tbl <- dplyr::tibble(param_name = extract_param_names(fit_nls),
                                      start_vals = tidyr::replace_na(start_vals_i, 0),
                                      param_est = sum_fit_nls$parameters[1:model_i$rTPC_n_params, 1],
                                      param_se = sum_fit_nls$parameters[1:model_i$rTPC_n_params, 2],
                                      model_name = i,
                                      model_AIC = AIC(fit_nls),
                                      model_BIC = BIC(fit_nls),
                                      model_fit = list(fit_nls))
      }, # <- inside tryCatch
      error = function(e) e)
      if(inherits(possible_error, "error")) {
        fit_nls <- NULL
      }
      if(is.null(fit_nls)) {list_param <- list_param}
      else {list_param <- list_param |>
        dplyr::bind_rows(list_param_tbl)}
    } # <- loop ends
  return(list_param)
}



## 2. TPC uncertainty simulation ----

sim_tpc_gridparams <- function(grid_parameters,
                               temperature, 
                               model_2sim){
  model_eq <- rTPC_choice_tbl |>
    filter(rTPC_models == model_2sim) |> 
    pull(nls_formula_params)
  params_i <- grid_parameters
  
  tpc_sim_i <- purrr::map(.x =  temperature,
                          .f = reformulate(termlabels = model_eq)
  )
  tpc_sim_tbl <- tibble(temperature,
                        pred_rate = tpc_sim_i) |>
    mutate(pred_rate = unlist(pred_rate))
  return(tpc_sim_tbl)
}

round_roots_2sim <- function(n_simulations, n_parameters){
  n_reps <- round(nthroot(n_simulations,
                          n_parameters))
  if(n_reps < n_simulations) {
    nreps <- n_reps + 1
  } else {nreps}
  return(nreps)
}

sim_tpcs_uncertainty <- function(fitted_parameters,
                                 model_2sim,
                                 temperature,
                                 rate,
                                 study,
                                 species,
                                 aic_model,
                                 number_simulations = 100,
                                 return_object = "plot") {
  therm_perf_data <- tibble(temperature, rate)
  model_ext_list <- fitted_parameters |>
    filter(model_name == model_2sim) |>
    slice(1) |>
    pull(model_fit)
  model_extracted <- model_ext_list[[1]]
  sum_model<- summary(model_extracted)
  n_parameters_model <- nrow(sum_model$coefficients)
  params_model <- sum_model$parameters[,1]
  
  if(number_simulations == 1000) {
    n_sims <- case_when(n_parameters_model == 2 ~  round_roots_2sim(1000, 2),
                        n_parameters_model == 3 ~  round_roots_2sim(1000, 3),
                        n_parameters_model == 4 ~  round_roots_2sim(1000, 4),
                        n_parameters_model == 5 ~  round_roots_2sim(1000, 5),
                        n_parameters_model == 6 ~  round_roots_2sim(1000, 6),
                        n_parameters_model == 7 ~  round_roots_2sim(1000, 7)
                        )
    
  } else {
    n_sims <- case_when(n_parameters_model == 2 ~  round_roots_2sim(100, 2),
                        n_parameters_model == 3 ~  round_roots_2sim(100, 3),
                        n_parameters_model == 4 ~  round_roots_2sim(100, 4),
                        n_parameters_model == 5 ~  round_roots_2sim(100, 5),
                        n_parameters_model == 6 ~  round_roots_2sim(100, 6),
                        n_parameters_model == 7 ~  round_roots_2sim(100, 7)
                        )
  } ## avoid extremely time-costly simulations
 
  allparams_sim_list <- list(simulations = NULL)
  set.seed(2023)
  for(parameter_index in 1:n_parameters_model){
    param_i_est <- sum_model$parameters[parameter_index, 1]
    param_i_se <- sum_model$parameters[parameter_index, 2]
    param_i_tbl <- tibble(estimate = param_i_est,
                          se = param_i_se)
    params_tbl_i <- rnorm(n = n_sims,
                          mean = param_i_tbl$estimate,
                          sd = param_i_tbl$se
    )
    allparams_sim_list[[parameter_index]] <- params_tbl_i
  }
  params_grid_raw <- expand.grid(allparams_sim_list)
  set.seed(2023)
  params_grid <- params_grid_raw |>
    slice_sample(n = number_simulations)
  params_names <- extract_param_names(model_extracted)
  colnames(params_grid) <- params_names
  all_tpcs_sim <- tibble(temperature = NULL,
                         dev_rate = NULL,
                         n_sim = NULL)
  
  
  for(n_sim_curve in c(1:nrow(params_grid))){
    print(paste0("Simulating Uncertainty TPC:","   ", n_sim_curve,"/",nrow(params_grid)))
    params_grid_i <- params_grid |>
      slice(n_sim_curve)
    curve_sim_i <- sim_tpc_gridparams(grid_parameters = params_grid_i,
                                      temperature = seq(min(temperature)-20,
                                                        max(temperature) +20,
                                                        .1),
                                      model_2sim = model_2sim)
    curve_sim_isim <- curve_sim_i |>
      mutate(n_sim = n_sim_curve,
             color = "gray78",
             linewidth = 0.5,
             alpha = 0.15)
    all_tpcs_sim <- all_tpcs_sim |>
      bind_rows(curve_sim_isim)
  }
  central_tpc <- sim_tpc_gridparams(grid_parameters = params_model,
                                    temperature = seq(min(temperature -10),
                                                      max(temperature +5),
                                                      0.1),
                                    model_2sim = model_2sim) |>
    mutate(color = "goldenrod",
           linewidth = 1.5,
           alpha = 1)
  simulated_tpcs <- bind_rows(all_tpcs_sim,
                              central_tpc) |> 
    filter(pred_rate >= 0,
           pred_rate < 1.5*max(rate))
    
  
  n_sim <- simulated_tpcs |> distinct(n_sim) |> pull(n_sim) |> length()
  
  plot_all_curves <- ggplot()+
    xlim(c(min(temperature) - 9,
           max(temperature) + 9))+
    geom_line(data = simulated_tpcs,
              aes(x = temperature,
                  y = pred_rate,
                  color = as_factor(color),
                  linewidth = linewidth,
                  group = n_sim,
                  alpha = alpha))+
    geom_point(data = therm_perf_data,
               aes(x = temperature, y = rate),
               color = "gray32",
               size = 3)+
    scale_color_identity()+
    scale_alpha_identity()+
    scale_linewidth_identity()+
    labs(title = species,
         subtitle = paste(model_2sim, "AIC =", round(aic_model,2)),
         tag = study,
         x = "Temperature (ºC)",
         y = expression(intrinsic~rate~of~increase~(italic(r)[m])))+
    ggtitle(parse(text = paste0("bolditalic('", species, "')")))+
  ggthemes::theme_clean()
if(return_object == "tbl") {
  return(all_tpcs_sim)
} else {
    return(plot_all_curves)
  }
}


# 3. Plot TPCs to help model selection ------------------------------------

plot_tpcs <- function(temperature, 
                      rate, 
                      fitted_parameters, 
                      study,
                      species,
                      aic_model,
                      facets = TRUE, 
                      uncertainty = FALSE, 
                      uncertainty_sims = 100) {
  
  predict2fill <- tibble(temperature = NULL,
                         rate = NULL,
                         model_name = NULL,
                         model_AIC = NULL)
  model_names2plot <- fitted_parameters |>
    distinct(model_name) |>
    pull(model_name)
  therm_perf_df <- tibble(temperature, rate)
  
  for(i in model_names2plot){
    fitted_parameters_i <- fitted_parameters |>
      filter(model_name == i)
    model_AIC_i <-fitted_parameters_i |>
      pull(model_AIC)
    params_i <- fitted_parameters_i |>
      pull(param_est)
    formula_i <- rTPC_choice_tbl |> 
      filter(rTPC_models == i) |>
      pull(nls_formula_params)
    ##predict based on parameters
    explore_preds <- tibble(temperature = seq(min(therm_perf_df$temperature)-10,
                                              max(therm_perf_df$temperature) +5,
                                              .1),
                            model_name = i,
                            model_AIC = model_AIC_i[1],
                            preds = NULL,
                            n_params = length(params_i))
    fit_vals_tbl <- explore_preds |>
      dplyr::select(temperature, model_name, model_AIC, n_params) |>
      mutate(formula = formula_i) |>
      mutate(preds = purrr::map_dbl(.x = temperature,
                                    .f = reformulate(unique(formula_i)))) |>
      filter(preds >= 0) |>
      dplyr::select(-formula)
    predict2fill <- predict2fill |>
      bind_rows(fit_vals_tbl)
  }
  
  aic_text <-  predict2fill  |>
    group_by(model_name)  |>
    summarise(aic = mean(model_AIC),
              n_params = paste(mean(n_params), "parameters"))  |>
    arrange(aic)
  aic_order <- aic_text  |>
    pull(model_name)
  aic_values <- aic_text |>
    mutate(aic =   paste("AIC =",
                         round(aic, 2)),
           temperature = min(therm_perf_df$temperature),
           preds = 1.5*max(therm_perf_df$rate))
  if(facets == FALSE && uncertainty == TRUE){
    plot_tpc_uncertainty <- sim_tpcs_uncertainty(fitted_parameters = fitted_parameters,
                                                 model_2sim = model_names2plot,
                                                 temperature,
                                                 rate,
                                                 study = study,
                                                 aic_model,
                                                 species = species,
                                                 number_simulations = uncertainty_sims)
    return(plot_tpc_uncertainty)
  } else {
    ggplot_models <- ggplot()+
      geom_line(data = predict2fill |>
                  filter(preds < (1.5*max(therm_perf_df$rate))),
                aes(x = temperature, 
                    y = preds,
                    color = model_name, 
                    group = model_name),
                linewidth = 1.3)+
      geom_point(data = therm_perf_df, aes(x = temperature,
                                           y = rate),
                 color = "gray32",
                 size = 2)+
      facet_wrap(~factor(model_name, levels = aic_order))+
      ggthemes::theme_few()+
      theme(legend.position = "none")+
      labs(title = study,
           x = "Temperature (ºC)",
           y = expression(intrinsic~rate~of~increase~(italic(r)[m])))+
      geom_label(data = aic_values,
                 aes(label = aic,
                     x = temperature,
                     y = preds,
                     fill = model_name),
                 color = "white",
                 size = 3)
    return(ggplot_models)
  }
  
}

# 4. Simulate linear models ----
linear_fun <- function(slope, intercept, x){
  y = slope*x + intercept
  return(y)
}
line_plot <- function(plot, simulation_tbl){
  plot +
    geom_abline(aes(slope = slope,
                    intercept = intercept))
}

sim_and_plot_linears <- function(model_object,
                                 var_x,
                                 var_y,
                                 n_sims,
                                 your_title,
                                 your_subtitle,
                                 lab_x,
                                 lab_y,
                                 color_points,
                                 color_central,
                                 color_uncertainty) {
  
  modelled_relation <- tibble(indep_var = var_x, 
                              dep_var = var_y)
  sum_model <- summary(model_object)
  model_params <- tibble(intercept = sum_model$coefficients[1,1],
                         slope =  sum_model$coefficients[2,1],
                         intercept_se = sum_model$coefficients[1,2],
                         slope_se =  sum_model$coefficients[2,2],
                         sample_size = nrow(modelled_relation)) |> 
    mutate(intercept_sd = intercept_se * sqrt(sample_size),
           slope_sd = slope_se * sqrt(sample_size))
  set.seed(2023)
  sim_slope_model <- rnorm(n_sims, 
                           mean = model_params |> select(slope) |> as_vector(),
                           sd = model_params |> select(slope_se) |> as_vector())
  
  sim_intercept_model <- rnorm(n_sims, 
                               mean = model_params |> select(intercept) |> as_vector(),
                               sd = model_params |> select(intercept_se) |> as_vector())
  
  sim_model_lm <- tibble(slope = sim_slope_model,
                         intercept = sim_intercept_model)
  
  n_draws <- n_sims
  alpha_level <- .15
  col_draw <- color_uncertainty
  col_median <-  color_central
  
  model_plot <- ggplot(modelled_relation,
                       aes(x = indep_var, y = dep_var))+
    labs(title= your_title,
         subtitle= your_subtitle,
         x = lab_x,
         y = lab_y)+
    ggthemes::theme_few()+ 
    geom_abline(aes(slope = slope,
                    intercept =intercept),
                data = slice_sample(sim_model_lm, n = n_draws),
                color = col_draw,
                alpha = alpha_level)+
    geom_abline(aes(slope = model_params$slope,
                    intercept = model_params$intercept),
                color = color_central,
                linewidth = 1.4)+
    geom_point(alpha=0.82,
               color= color_points)+
    theme(plot.title = element_text(face = "bold"))
  
  model_plot
  
}

plot_preds_hier_model <- function(data2fit, model_fit) {
  ### example
  families_3points <- therm_traits_order |> 
    group_by(family) |> 
    count() |> 
    filter(n >2) |> 
    pull(family)
  therm_traits_family <- therm_traits |> 
    filter(family %in% families_3points)
  model_fit <- lmerTest::lmer(ctmax ~ abs(lat) + (1|family),
                         data = therm_traits)
  coef(model_fit)
  
  
  ###
  groups_coef <- coef(model_fit)[[1]] |> 
    rownames_to_column(var = "grouping_var") |> 
    rename(intercept = `(Intercept)`,
           slope = `abs(lat)`)
  preds_groups <- tibble(group = NULL,
                         var_x = NULL,
                         var_y = NULL)
  for(group_i in unique(groups_coef$grouping_var)) {
    intercept_i <- groups_coef |> 
      filter(grouping_var == group_i) |> 
      pull(intercept)
    slope_i <- groups_coef |> 
      filter(grouping_var == group_i) |> 
      pull(slope)
    preds_tbl <- tibble(group = group_i,
                        var_x = seq(min(data2pred$var_x),
                                    max(data2pred$var_x),
                                    length.out = 100)) |> 
      mutate(var_y = map_dbl(.x = var_x,
                             .f = ~linear_fun(slope = slope_i,
                                              intercept = intercept_i,
                                              x = .x)
                             )
             )
    preds_groups <- bind_rows(preds_groups, preds_tbl)
  }
 ggplot(preds_groups, aes(x = var_x,
                          y = var_y))+
   geom_line(aes(color = group))
  ###
}


# 5. Predict future rates ----

predict_r <- function(fitted_model, temp){
sum_fitted_model <- summary(fitted_model)
model_name <-   stringr::str_extract(as.character(sum_fitted_model$formula),
                                    "\\w+(?=\\()")[[3]]

params_i <- sum_fitted_model$parameters[,1]
therm_perf_df <- tibble(temperature = temp, 
                        pred_rate = NULL)

  formula_i <- rTPC_choice_tbl |> 
    filter(rTPC_models == model_name) |>
    pull(nls_formula_params)
  ##predict based on parameters
  explore_preds <- tibble(temperature = temp,
                          model_name = model_name,
                          preds = NULL,
                          n_params = length(params_i))
  pred_rate <- explore_preds |>
    dplyr::select(temperature, model_name, n_params) |>
    mutate(formula = formula_i) |>
    mutate(preds = purrr::map_dbl(.x = temperature,
                                  .f = reformulate(unique(formula_i)))) |>
    pull(preds)
  return(pred_rate)
}

# 6. bootstrap params uncertainty ----

boots_uncertainty <- function(model_fit, temperature, rate, return_object = "tbl") {
  data2refit <- tibble(temperature, rate) 
  nlms_model_fit <- model_fit
  model_reformula <- deparse(formula(nlms_model_fit)[[3]])
  first_parenthesis_position <- str_locate(formula_as_character, "\\(")[1, "start"]
  model_name_refit <- str_sub(model_reformula, end = first_parenthesis_position - 2)
  lower_lims <- get_lower_lims(data2refit$temperature,
                               data2refit$rate,
                               model_name = model_name_refit)
  upper_lims <- get_upper_lims(data2refit$temperature,
                               data2refit$rate,
                               model_name = model_name_refit)
  
  fit_nlsLM <- minpack.lm::nlsLM(formula = formula(nlms_model_fit),
                                 data = data2refit,
                                 start = coef(nlms_model_fit),
                                 lower = lower_lims,
                                 upper = upper_lims,
                                 weights = rep(1, times = length(temperature)),
                                 control = minpack.lm::nls.lm.control(maxiter = 200))
  
  extra_params <- calc_params(fit_nlsLM) %>%
    pivot_longer(everything(), names_to =  'param', values_to = 'estimate') 
  set.seed(2020)
  ci_extra_params <- Boot(fit_nlsLM, 
                          f = function(x){unlist(calc_params(x))}, 
                          labels = names(calc_params(fit_nlsLM)),
                          R = 200, 
                          method = 'case') %>%
    confint(., method = 'bca') %>%
    as.data.frame() %>%
    rename(conf_lower = 1, conf_upper = 2) %>%
    rownames_to_column(., var = 'param') %>%
    mutate(method = 'case bootstrap')
  
  ci_extra_params <- left_join(ci_extra_params, extra_params)
  
  point_range_plot <-  ggplot(ci_extra_params, aes(param, estimate)) +
    geom_point(size = 4) +
    geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
    theme_bw() +
    facet_wrap(~param, scales = 'free') +
    scale_x_discrete('') +
    labs(title = 'Uncertainty of calculated parameters',
         subtitle = expression(Bootstrap~(see~Padfield~italic(et~al.)~2021)))
  return(ci_extra_params)
}

