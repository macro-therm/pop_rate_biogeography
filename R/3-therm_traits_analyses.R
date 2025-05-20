# Aim of the script: thermal traits (meta-)analyses


# 0. Packages and data loading --------------------------------------------
library(here)
library(tidyverse)
library(leaflet)
library(RColorBrewer)
library(sf)
library(nlme)
library(lme4)
library(emmeans)
library(ggthemes)
library(ggdist)
library(wesanderson)
source(here("analyses/S1-functions.R"))

##read traits
therm_traits <- readRDS(here("data/data_sink/therm_traits_intrate.rds")) |> 
  mutate(across(c(reference, species, order, family, feed_guild),
                ~as_factor(.x))) |> 
  filter(ctmin > -5, #unrealistic fits ( <- maybe we could select other model for these such as gaussian)
         ctmax <55) |> 
  mutate(therm_range = topt-ctmin) |> 
  mutate(log_rmax = log(rmax))


therm_traits_nogaussians <- therm_traits |> 
  filter(model_name != "gaussian_1987")

therm_traits_insects <- therm_traits |> 
  filter(order != "Acari")

# 1. Explorating thermal traits--------------------------------------------

## a) tpc model --------------------------------------------

### all
count_best_models <- therm_traits |> 
  count(model_name) |> 
  mutate(perc = map_dbl(.x = n,
                        .f = ~.x*100/sum(n))
         )
ggplot(count_best_models, aes(x = fct_reorder(as_factor(model_name), n), y = n))+
  geom_bar(stat = "identity", aes(fill = model_name))+
  scale_fill_manual(values = MoMAColors::MoMAPalettes$Klein[[1]])+
  labs(x = NULL, y = "Number of selections",
       title = "TPC models ranking")+
  coord_flip()+
  ggthemes::theme_clean()+
  theme(legend.position = "none")

### by order
count_best_models_order <- therm_traits |> 
  group_by(order) |> 
  count(model_name) |> 
  mutate(perc = map_dbl(.x = n,
                        .f = ~.x*100/sum(n))
         )

ggplot(count_best_models_order, aes(x = fct_reorder(as_factor(model_name), n), y = n))+
  geom_bar(stat = "identity", aes(fill = model_name))+
  scale_fill_manual(values = MoMAColors::MoMAPalettes$Klein[[1]])+
  labs(x = NULL, y = "Number of selections",
       title = "TPC models ranking")+
  coord_flip()+
  facet_wrap(~order)+
  ggthemes::theme_clean()+
  theme(legend.position = "none")
  
## b) taxonomic  --------------------------------------------
count_order <- therm_traits |> 
  count(order) |> 
  mutate(order = as_factor(order)) |> 
  mutate(perc = map_dbl(.x = n,
                        .f = ~.x*100/sum(n))
  )
         
         
ggplot(data = count_order, aes(x = fct_reorder(order, n),
                                           y = n))+
  geom_bar(aes(fill = order), stat = "identity")+
  scale_fill_manual(values = MoMAColors::MoMAPalettes$Klein[[1]])+
  theme_bw()+
  coord_flip()+
  labs(x = NULL,
       y = "number of effect sizes",
       title = "Taxonomic orders included")+
  ggthemes::theme_clean()+
  theme(legend.position = "none")
ggsave(here("data/data_sink/figs/count_order.png"), width = 15, height = 15, units = "cm")

## c) feeding guild  --------------------------------------------

count_feed_guild <- therm_traits |> 
  count(feed_guild) |> 
  mutate(feed_guild = as_factor(feed_guild)) |> 
  mutate(perc = map_dbl(.x = n,
                        .f = ~.x*100/sum(n))
  )

### all
ggplot(data = count_feed_guild, aes(x = fct_reorder(feed_guild, n),
                               y = n))+
  geom_bar(aes(fill = feed_guild), stat = "identity")+
  scale_fill_manual(values = MoMAColors::MoMAPalettes$Klein[[1]])+
  coord_flip()+
  labs(x = NULL,
       y = "number of effect sizes",
       title = "Feeding guilds included")+
  ggthemes::theme_clean()+
  theme(legend.position = "none")
ggsave(here("data/data_sink/figs/count_feed_guild.png"), width = 15, height = 15, units = "cm")

### feeding guild within orders
count_feed_guild_within_orders <- therm_traits |> 
  count(order,feed_guild) |> 
  mutate(feed_guild = as_factor(feed_guild)) 
  
ggplot(data = count_feed_guild_within_orders, aes(x = fct_reorder(feed_guild, n),
                                    y = n))+
  geom_bar(aes(fill = feed_guild), stat = "identity")+
  scale_fill_manual(values = MoMAColors::MoMAPalettes$Klein[[1]])+
  coord_flip()+
  facet_wrap(~order)+
  labs(x = NULL,
       y = "number of effect sizes",
       title = "Feeding guilds within orders")+
  theme_bw()+
  theme(legend.position = "none")

### orders within feeding_guilds
count_orders_within_feed_guild <- therm_traits |> 
  count(feed_guild, order) |> 
  mutate(feed_guild = as_factor(feed_guild),
         order = as_factor(order)) 

ggplot(data = count_orders_within_feed_guild, aes(x = fct_reorder(order, n),
                                                  y = n))+
  geom_bar(aes(fill = feed_guild), stat = "identity")+
  scale_fill_manual(values = MoMAColors::MoMAPalettes$Klein[[1]])+
  coord_flip()+
  facet_wrap(~feed_guild)+
  labs(x = NULL,
       y = "number of effect sizes",
       title = "Feeding guilds within orders")+
  theme_bw()+
  theme(legend.position = "none")


## combination 1
n_order <- therm_traits |>
  count(order, name = "n_order")

count_composed <- therm_traits |>
  group_by(order) |> 
  count(feed_guild) |> 
  inner_join(n_order)

ggplot(data = count_composed, aes(x = fct_reorder(order, n_order), y = n))+
  geom_bar(position = "stack", 
           stat = "identity",
           aes(fill = feed_guild))+
  scale_fill_manual(values = MoMAColors::MoMAPalettes$Klein[[1]])+
  ggthemes::theme_clean()+
  coord_flip()+
  labs(fill = "Feeding guild",
       title = "Order - Feeding guild composition",
       x = NULL,
       y = "Number of effect sizes")

## combination 2
n_fg <- therm_traits |>
  count(feed_guild, name = "n_fg")

count_composed2 <- therm_traits |>
  group_by(feed_guild) |> 
  count(order) |> 
  inner_join(n_fg)

ggplot(data = count_composed2, aes(x = fct_reorder(feed_guild, n_fg), y = n))+
  geom_bar(position = "stack", 
           stat = "identity",
           aes(fill = order))+
  scale_fill_manual(values = MoMAColors::MoMAPalettes$Klein[[1]])+
  ggthemes::theme_clean()+
  coord_flip()+
  labs(fill = "Order",
       title = "Feeding guild - Order composition",
       x = NULL,
       y = "Number of effect sizes")
  
## d) locations and latitude --------------------------------------------

##normal map
map_studies <- ggplot(data = therm_traits, aes(x = lon, y = lat))+
  borders("world", colour = "transparent", fill = "lightgrey")+
  geom_point(aes(color = order), 
             alpha = 0.25,
             position = position_jitter(width = 1,
                                        height = 1),
             size = 3)+
  ggthemes::theme_map()+
  theme(legend.position = "bottom")
map_studies
ggsave(here("data/data_sink/figs/map_studies.png"),
       height=15,
       width=25,
       units="cm",
       bg = "transparent")


##leaflet map

mytext <- paste(
  "<i>Species</i>: ", therm_traits$species,"<br/>", 
  "<i>Reference</i>: ", therm_traits$reference, "<br/>", 
  "<i>Order</i>: ", therm_traits$order, "<br/>",
  "<i>Family</i>: ", therm_traits$family, "<br/>",
  "<i>Feeding guild</i>: ", therm_traits$feed_guild, "<br/>",
  "<i>TPC model</i>: ", therm_traits$model_name, "<br/>",
  sep="") |>
  lapply(htmltools::HTML)

colors_orders <- MoMAColors::MoMAPalettes$Klein[[1]]

palette_shift <- colorFactor(palette = colors_orders,
                              domain = unique(therm_traits$order))

intrapests_orders_locations<- leaflet(data = therm_traits) |> 
  addTiles() |> 
  addCircleMarkers(lng = ~jitter(lon, factor = 1), 
                   lat = ~jitter(lat, factor = 1), 
                   stroke = FALSE,
                   fillColor = ~palette_shift(order), 
                   fillOpacity = 0.8, 
                   color = ~palette_shift(order), 
                   popup = ~reference,
                   label = mytext,
                   group = "reference",
                   labelOptions = labelOptions( 
                     style = list("font-weight" = "normal", 
                                  padding = "3px 8px"), 
                     textsize = "13px", 
                     direction = "auto")) |> 
  addProviderTiles('CartoDB.DarkMatterNoLabels') |>
  addLegend(pal= palette_shift, 
            values= ~order, 
            opacity = 0.9,
            title = "Taxonomic order",
            position = "bottomleft")
intrapests_orders_locations
htmlwidgets::saveWidget(intrapests_orders_locations,
                        here("data/data_sink/figs/map_locations_interactive.html"))

## custom map

map_source_locations <- ggplot(data = therm_traits, aes(x = lon, y = lat))+
  borders("world", colour = "transparent", fill = "gray")+
  geom_point(aes(color = order), 
             position = position_jitter(width = 1,
                                        height = 1),
             alpha = 0.8,
             size = 3)+
  scale_color_manual(values = colors_orders)+
  ggthemes::theme_clean()+
  theme(legend.position = "bottom",
        panel.background = element_rect(fill = "gray93"))+
  labs(y = "Latitude",
       x = "Longitude",
       color = "Taxonomic order")
map_source_locations
ggsave(here("data/data_sink/figs/map_studies.png"),
       height=15,
       width=25,
       units="cm",
       bg = "transparent")
## e) body sizes --------------------------------------------
therm_traits_bodysize <- therm_traits |> 
  mutate(log_body_length_mm = log(body_length_mm))

### order distribution
log_size <- expression(log[10] ~ italic("body size"))

ggplot(therm_traits_bodysize, aes(x = order, y = log_body_length_mm, 
                                  color = order,
                                  fill = order))+
  ggdist::stat_halfeye(adjust = .5, width = .9, .width = 0.01, 
                       justification = -0.3, point_colour = NA) + 
  geom_boxplot(width = .5, outlier.shape = NA,
               color = "gray23",
               fill = NA)+
  scale_fill_manual(values = MoMAColors::MoMAPalettes$Klein[[1]])+
  scale_color_manual(values = MoMAColors::MoMAPalettes$Klein[[1]])+
  ggdist::stat_dots(side = "left", dotsize = .6, 
                    justification = 1.3, 
                    binwidth = .1)+
  coord_flip()+
  theme_bw()+
  labs(x = log_size,
       y = "Taxonomic order")+
  theme(legend.position = "none")
ggsave(here("data/data_sink/figs/bodysize_order.png"), width = 15, height = 15, units = "cm")

### feed_guild distribution
log_size <- expression(log[10] ~ italic("body size"))
ggplot(therm_traits_bodysize, aes(x = feed_guild, y = log_body_length_mm, 
                                  color = feed_guild,
                                  fill = feed_guild))+
  ggdist::stat_halfeye(adjust = .5, width = .9, .width = 0.01, 
                       justification = -0.3, point_colour = NA) + 
  geom_boxplot(width = .5, outlier.shape = NA,
               color = "gray23",
               fill = NA)+
  ggdist::stat_dots(side = "left", dotsize = .6, 
                    justification = 1.3, 
                    binwidth = .1)+
  scale_fill_manual(values = MoMAColors::MoMAPalettes$Klein[[1]])+
  scale_color_manual(values = MoMAColors::MoMAPalettes$Klein[[1]])+
  coord_flip()+
  theme_bw()+
  labs(x = log_size,
       y = "Feeding guild")+
  theme(legend.position = "none")
ggsave(here("data/data_sink/figs/bodysize_feed.png"), width = 15, height = 15, units = "cm")



### overall distribution
ggplot(therm_traits_bodysize, aes(x = log_body_length_mm))+

  geom_histogram(aes(y = after_stat(density)),
                 fill = "goldenrod3",
                 color = "darkslateblue"
                 )+
  ggthemes::theme_few()                
           
summary(therm_traits_bodysize$log_body_length_mm) 

## f) normality --------------------------------------------

therm_traits_vertical <- therm_traits |> 
  pivot_longer(cols = c(ctmin, ctmax, topt, rmax, e, eh, q10, thermal_safety_margin,
                        thermal_tolerance, therm_range, breadth, skewness),
               names_to = "parameter",
               values_to = "estimate")

plot_param <-  ggplot(therm_traits_vertical, aes(x = estimate))+
  geom_histogram(aes(y = after_stat(density),
                     fill = parameter),
                 color = "gray59",
                 bins = 15)+
  facet_wrap(~parameter,scales = "free")+
  scale_fill_manual(values = c(MoMAColors::MoMAPalettes$Klein[[1]], "slateblue"))+
  theme_few()+
  theme(legend.position  = "none")
plot_param
ggsave(filename = here("data/data_sink/figs/params_histograms.png"),
       width = 22, height = 19, units = "cm")


# 2. Summary effects--------------------------------------------

parameters_intercept_output <- tibble(pred = NULL,
                                      ci_l = NULL,
                                      ci_h = NULL,
                                      predint_l = NULL,
                                      predint_h = NULL,
                                      parameter = NULL)

parameters_order_output <- tibble(pred = NULL,
                                  ci_l = NULL,
                                  ci_h = NULL,
                                  predint_l = NULL,
                                  predint_h = NULL,
                                  parameter = NULL)

parameters_feed_guild_output <- tibble(pred = NULL,
                                          ci_l = NULL,
                                          ci_h = NULL,
                                          predint_l = NULL,
                                          predint_h = NULL,
                                          parameter = NULL)


## a) ct_min --------------------------------------------
#### i. intercept ----
ctmin_intercept_lmer_ref <- lmer(ctmin ~ 1 + (1|reference),
                             na.action = na.omit,
                             data = therm_traits)
sum_ctmin_ref_lmer <- summary(ctmin_intercept_lmer_ref)

#### ii. order ----
ctmin_order_lmer_ref <- lmerTest::lmer(ctmin ~ order + (1|reference),
                             data = therm_traits,
                             na.action = na.omit)

summary(ctmin_order_lmer_ref)
anova(ctmin_order_lmer_ref)

#### iii. feeding guild ----

ctmin_feed_guild_lmer_ref <- lmerTest::lmer(ctmin ~ feed_guild + (1|reference),
                                       data = therm_traits,
                                       na.action = na.omit)
summary(ctmin_feed_guild_lmer_ref)
anova(ctmin_feed_guild_lmer_ref)

intervals_ctmin_order <- bernr::bolker_ci(ctmin_order_lmer_ref,
                                          newdat = data.frame(order = levels(therm_traits$order)),
                                          pred_int = TRUE)

ctmin_order_output <- intervals_ctmin_order|> 
  select(pred, ci_l, ci_h, predint_l, predint_h) |>  
  mutate(order = levels(therm_traits$order)) |>  
  as_tibble() |> 
  mutate(parameter = rep("ctmin", length(levels(therm_traits$order))))

parameters_order_output <- bind_rows(parameters_order_output, ctmin_order_output)
 
## b) ct_max  --------------------------------------------
#### i. intercept ----
ctmax_intercept_lmer_ref <- lmer(ctmax ~ 1 + (1|reference),
                                 na.action = na.omit,
                                 data = therm_traits)
sum_ctmax_ref_lmer <- summary(ctmax_intercept_lmer_ref)
#### ii. order ----
ctmax_order_lmer_ref <- lmerTest::lmer(ctmax ~ order + (1|reference),
                                       data = therm_traits_nogaussians,
                                       na.action = na.omit)

summary(ctmax_order_lmer_ref)
anova(ctmax_order_lmer_ref)

#### iii. feeding guild ----
ctmax_feed_guild_lmer_ref <- lmerTest::lmer(ctmax ~ feed_guild + (1|reference),
                                            data = therm_traits_nogaussians,
                                            na.action = na.omit)

summary(ctmax_feed_guild_lmer_ref)
anova(ctmax_feed_guild_lmer_ref)
intervals_ctmax_order <- bernr::bolker_ci(ctmax_order_lmer_ref,
                                          newdat = data.frame(order = levels(therm_traits$order)),
                                          pred_int = TRUE)

ctmax_order_output <- intervals_ctmax_order|> 
  select(pred, ci_l, ci_h, predint_l, predint_h) |>  
  mutate(order = levels(therm_traits$order)) |>  
  as_tibble() |> 
  mutate(parameter = rep("ctmax", length(levels(therm_traits$order))))

parameters_order_output <- bind_rows(parameters_order_output, ctmin_order_output)
## c) t_opt  --------------------------------------------
#### i. intercept ----
topt_intercept_lmer_ref <- lmerTest::lmer(topt ~ 1 + (1|reference),
                                 na.action = na.omit,
                                 data = therm_traits)
sum_topt_ref_lmer <- summary(topt_intercept_lmer_ref)

#### ii. order ----
topt_order_lmer_ref <- lmerTest::lmer(topt ~ order + (1|reference),
                                       data = therm_traits,
                                       na.action = na.omit)

summary(topt_order_lmer_ref)
anova(topt_order_lmer_ref)

#### iii. feeding guild ----
topt_feed_guild_lmer_ref <- lmerTest::lmer(topt ~ feed_guild + (1|reference),
                                            data = therm_traits,
                                            na.action = na.omit)

summary(topt_feed_guild_lmer_ref)
anova(topt_feed_guild_lmer_ref)

intervals_topt_order <- bernr::bolker_ci(topt_order_lmer_ref,
                                          newdat = data.frame(order = levels(therm_traits$order)),
                                          pred_int = TRUE)

topt_order_output <- intervals_topt_order|> 
  select(pred, ci_l, ci_h, predint_l, predint_h) |>  
  mutate(order = levels(therm_traits$order)) |>  
  as_tibble() |> 
  mutate(parameter = rep("topt", length(levels(therm_traits$order))))

parameters_order_output <- bind_rows(parameters_order_output, topt_order_output)

## d) r_max --------------------------------------------
#### i. intercept ----
rmax_intercept_lmer_ref <- lmer(rmax ~ 1 + (1|reference),
                                na.action = na.omit,
                                data = therm_traits)
summary(rmax_intercept_lmer_ref)


#### ii. order ----
rmax_order_lmer_ref <- lmerTest::lmer(log(rmax) ~ order + (1|reference),
                                      data = therm_traits,
                                      na.action = na.omit)

summary(rmax_order_lmer_ref)
anova(rmax_order_lmer_ref)

#### iii. feeding guild ----
rmax_feed_guild_lmer_ref <- lmerTest::lmer(log(rmax) ~ feed_guild + (1|reference),
                                            data = therm_traits,
                                            na.action = na.omit)

summary(rmax_feed_guild_lmer_ref)
anova(rmax_feed_guild_lmer_ref)

## e) thermal_tolerance  --------------------------------------------
#### i. intercept ----
tolerance_intercept_lmer_ref <- lmerTest::lmer(thermal_tolerance ~ 1 + (1|reference),
                                          na.action = na.omit,
                                          data = therm_traits_nogaussians)
summary(tolerance_intercept_lmer_ref)

#### ii. order ----
tolerance_order_lmer_ref_nogauss <- lmerTest::lmer(thermal_tolerance ~ 1 + order + (1|reference),
                                               na.action = na.omit,
                                               data = therm_traits_nogaussians)
summary(tolerance_order_lmer_ref_nogauss)

tolerance_order_lmer_ref <- lmerTest::lmer(thermal_tolerance ~ 1 + order + (1|reference),
                                                       na.action = na.omit,
                                                       data = therm_traits_nogaussians)
summary(tolerance_order_lmer_ref)
anova(tolerance_order_lmer_ref)

#### iii. feeding guild ----
tolerance_feed_guild_lmer_ref_nogauss <- lmerTest::lmer(thermal_tolerance ~ 1 + feed_guild + (1|reference),
                                                       na.action = na.omit,
                                                       data = therm_traits_nogaussians)
summary(tolerance_intercept_lmer_ref_nogauss)

tolerance_feed_guild_lmer_ref <- lmerTest::lmer(thermal_tolerance ~ 1 + feed_guild + (1|reference),
                                               na.action = na.omit,
                                               data = therm_traits_nogaussians)
summary(tolerance_feed_guild_lmer_ref)
anova(tolerance_feed_guild_lmer_ref)

## f) therm_range  --------------------------------------------
#### i. intercept ----
therm_window_lmer_ref <- lmerTest::lmer(therm_range ~ 1 + (1|reference),
                                               na.action = na.omit,
                                               data = therm_traits)
summary(therm_window_lmer_ref)


#### ii. order ----
therm_window_lmer_ref_order <- lmerTest::lmer(therm_range ~ 1 + order + (1|reference),
                                        na.action = na.omit,
                                        data = therm_traits)
summary(therm_window_lmer_ref_order)
anova(therm_window_lmer_ref_order)

#### iii. feeding guild ----
therm_window_lmer_ref_feed_guild <- lmerTest::lmer(therm_range ~ 1 + feed_guild + (1|reference),
                                              na.action = na.omit,
                                              data = therm_traits)
summary(therm_window_lmer_ref_feed_guild)
anova(therm_window_lmer_ref_feed_guild)


## g) breadth  --------------------------------------------
#### i. intercept ----
breadth_intercept_lmer_ref <- lmer(breadth ~ 1 + (1|reference),
                                na.action = na.omit,
                                data = therm_traits)
summary(breadth_intercept_lmer_ref)

#### 1. order ----
breadth_order_lmer_ref <- lmerTest::lmer(breadth ~ order + (1|reference),
                                         data = therm_traits,
                                         na.action = na.omit)

summary(breadth_order_lmer_ref)
anova(breadth_order_lmer_ref)

### 2. feed guild ----
breadth_feed_guild_lmer_ref <- lmerTest::lmer(breadth ~ feed_guild + (1|reference),
                                            data = therm_traits,
                                            na.action = na.omit)
summary(breadth_feed_guild_lmer_ref)
anova(breadth_feed_guild_lmer_ref)

##k) plot thermal traits ----
## add summary effect and uncertainty

params_maintraits <- params_ctmin_intercept_lmer_ref |> 
  bind_rows(params_ctmax_intercept_lmer_ref, params_topt_intercept_lmer_ref)

unique(therm_traits$feed_guild)
traitpoints_feed_guild_freq <- therm_traits |> 
  select(reference, feed_guild, species, ctmax, ctmin, topt) 

ctmax_hist <- ggplot(therm_traits) +
  ggdist::stat_dist_halfeye(aes(x = ctmax),breaks = 30,
                            color = "#9b2226",
                            fill = "#bb3e03",
  ) +   
  ggdist::stat_dots(aes(x = ctmax),
                    side = "left",
                    dotsize = 1.2,
                    justification = 1.05,
                    binwidth = .5,
                    color = "#9b2226",
                    fill = "#bb3e03",
                    alpha = .5
  )+
  theme_clean()+
  theme(axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        strip.text.x = element_blank())+
  labs(x = "Temperature (ºC)",
       y = "")

ctmin_hist <- ggplot(therm_traits) +
  ggdist::stat_dist_halfeye(aes(x = ctmin), breaks = 30,
                            adjust = .5,
                            width = .5, ## set slab interval to show IQR and 95% data range
                            color = "#005f73",
                            fill = "#0a9396"
  ) +
  ggdist::stat_dots(aes(x = ctmin),
                    side = "left",
                    dotsize = 1.25,
                    justification = 1.05,
                    binwidth = .5,
                    color = "#005f73",
                    fill = "#0a9396",
                    alpha = .5
  )+
  theme_clean()+
  theme(axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        strip.text.x = element_blank())+
  labs(x = "Temperature (ºC)",
       y = " ")


topt_hist <- ggplot(therm_traits) +
  ggdist::stat_dist_halfeye(aes(x = topt),breaks = 30,
                            adjust = .5,
                            width = .5, ## set slab interval to show IQR and 95% data range
                            color = "#E1AE01",
                            fill = "#EBCC2A"
  ) + 
  ggdist::stat_dots(aes(x = topt),
                    side = "left",
                    dotsize = .7,
                    justification = 1.05,
                    binwidth = .5,
                    color = "#E1AF00",
                    fill = "#EBCC2A",
                    alpha = .5
  )+
  theme_clean()+
  theme(axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        strip.text.x = element_blank())+
  labs(x = "Temperature (ºC)",
       y = "")

cowplot::plot_grid(ctmin_hist, topt_hist, ctmax_hist, 
                   nrow = 1, labels = c("CTmin", "Topt", "CTmax"))
ggsave(filename = here("data/data_sink/figs/rainclouds_traits.png"),
       width = 25, height = 9, dpi = 300, units = "cm")



# 3. Macrophysiology analyses--------------------------------------------
## a) traits ~ latitude --------------------------------------------
### a.1. ctmax - lat ----
ctmax_lat_ref <- lmerTest::lmer(ctmax ~ abs(lat) + (1|reference),
                       data = therm_traits)
summary(ctmax_lat_ref) 

therm_traits_deutsch <- therm_traits |> filter(model_name == "deutsch_2008")
ctmax_lat_ref_nogauss <-  lmerTest::lmer(ctmax ~ abs(lat) + (1|reference),
                           data = therm_traits_nogaussians)

summary(ctmax_lat_ref_nogauss)

## all models
ctmax_to_lat_plot <- sim_and_plot_linears(ctmax_lat,
                                          var_x = abs(therm_traits$lat),
                                          var_y = therm_traits$ctmax,
                                          n_sims = 1000,
                                          your_title = expression(CT[max]),
                                          your_subtitle = NULL,
                                          lab_x = "Absolute Latitude (º)",
                                          lab_y = "Temperature (ºC)",
                                          color_points = "#bb3e03",
                                          color_central = "#9b2226",
                                          color_uncertainty = "#ee9b00")

ctmax_to_lat_ref_plot <- sim_and_plot_linears(ctmax_lat_ref_nogauss,
                                          var_x = abs(therm_traits_nogaussians$lat),
                                          var_y = therm_traits_nogaussians$ctmax,
                                          n_sims = 1000,
                                          your_title = "Critical Thermal Maxima",
                                          your_subtitle = "Gaussian TPCs excluded",
                                          lab_x = "Absolute Latitude (º)",
                                          lab_y = "Temperature (ºC)",
                                          color_points = "#bb3e03",
                                          color_central = "#9b2226",
                                          color_uncertainty = "#ee9b00")

ctmax_to_lat_ref_deutsch_plot <- sim_and_plot_linears(ctmax_lat_ref_deutsch,
                                              var_x = abs(therm_traits_deutsch$lat),
                                              var_y = therm_traits_deutsch$ctmax,
                                              n_sims = 1000,
                                              your_title = expression(CT[max]),
                                              your_subtitle = "Only deutsch_2008 models",
                                              lab_x = "Absolute Latitude (º)",
                                              lab_y = "ºC",
                                              color_points = "#bb3e03",
                                              color_central = "#9b2226",
                                              color_uncertainty = "#ee9b00")

ggsave(ctmax_to_lat_ref_deutsch_plot,
       filename = here("data/data_sink/figs/traits_lat/ctmax_to_lat_ref_deutsch_plot.png"),
       width = 16, height = 16, units = "cm")

## deutsch models

ctmax_lat_deutsch_plot <- sim_and_plot_linears(ctmax_lat_deutsch,
                                               var_x = abs(therm_traits_deutsch$lat),
                                               var_y = therm_traits_deutsch$ctmax,
                                               n_sims = 1000,
                                               your_title = expression(CT[max]),
                                               your_subtitle = "Model: deutsch_2008",
                                               lab_x = "Absolute Latitude (º)",
                                               lab_y = "ºC",
                                               color_points = "#bb3e03",
                                               color_central = "#9b2226",
                                               color_uncertainty = "#ee9b00")

ggsave(ctmax_lat_deutsch_plot,
       filename = here("data/data_sink/figs/traits_lat/ctmax_lat_deutsch_plot.png"),
       width = 16, height = 16, units = "cm")

###### by order ----
for(order_i in levels(therm_traits$order)){
  therm_traits_order <- therm_traits |> 
    filter(order == order_i)
  possible_error <- tryCatch(expr = {ctmax_lat_order_i <-  lmer(ctmax ~ abs(lat) + (1|species),
                                                                data = therm_traits_order)
  ctmax_to_lat_plot_order <- sim_and_plot_linears(ctmax_lat_order_i,
                                                  var_x = abs(therm_traits_order$lat),
                                                  var_y = therm_traits_order$ctmax,
                                                  n_sims = 1000,
                                                  your_title = expression(CT[max]),
                                                  your_subtitle = order_i,
                                                  lab_x = "Absolute Latitude (º)",
                                                  lab_y = "ºC",
                                                  color_points = "#bb3e03",
                                                  color_central = "#9b2226",
                                                  color_uncertainty = "#ee9b00")
  },  # <- inside tryCatch
  error = function(e) e)
  if(inherits(possible_error, "error")) {
    ctmax_lat_order_i <- NULL
    ggplot(therm_traits_order, aes(x = abs(therm_traits_order$lat),
                                   y = therm_traits_order$ctmax))+
      geom_point(color = "#bb3e03")+
      geom_smooth(method = "lm",
                  fill = "#ee9b00",
                  color = "#9b2226")+
      ggthemes::theme_clean()+
      labs(x = "Absolute Latitude",
           y = "ºC",
           title = expression(CT[max]),
           subtitle = "No lmer fitting: using `geom_smooth` instead for exploratory analysis")
  }
  ggsave(filename = here(paste0("data/data_sink/figs/traits_lat/ctmax_lat_", order_i,"_plot.png")),
         width = 16, height = 16, units = "cm")
}

###### by feed_guild ----
for(feed_guild_i in levels(therm_traits$feed_guild)){
  therm_traits_feed_guild <- therm_traits |> 
    filter(feed_guild == feed_guild_i)
  possible_error <- tryCatch(expr = {ctmax_lat_feed_guild_i <-  lmer(ctmax ~ abs(lat) + (1|species),
                                                                data = therm_traits_feed_guild)
  ctmax_to_lat_plot_feed_guild <- sim_and_plot_linears(ctmax_lat_feed_guild_i,
                                                  var_x = abs(therm_traits_feed_guild$lat),
                                                  var_y = therm_traits_feed_guild$ctmax,
                                                  n_sims = 1000,
                                                  your_title = expression(CT[max]),
                                                  your_subtitle = feed_guild_i,
                                                  lab_x = "Absolute Latitude (º)",
                                                  lab_y = "ºC",
                                                 ctmax)
  },  # <- inside tryCatch
  error = function(e) e)
  if(inherits(possible_error, "error")) {
    ctmax_lat_feed_guild_i <- NULL
    ggplot(therm_traits_feed_guild, aes(x = abs(therm_traits_feed_guild$lat),
                                   y = therm_traits_feed_guild$ctmax))+
      geom_point(color = "#bb3e03")+
      geom_smooth(method = "lm",
                  fill = "#ee9b00",
                  color = "#9b2226")+
      ggthemes::theme_clean()+
      labs(x = "Absolute Latitude",
           y = "ºC",
           title = expression(CT[max]),
           subtitle = "No lmer fitting: using `geom_smooth` instead for exploratory analysis")
  }
  ggsave(filename = here(paste0("data/data_sink/figs/traits_lat/ctmax_lat_", feed_guild_i,"_plot.png")),
         width = 16, height = 16, units = "cm")
}

###### by size categories ----

therm_traits_sizecats <- therm_traits |> 
  mutate(log_body_length_mm = log10(body_length_mm),
         body_length_cat = case_when(log_body_length_mm < 0 ~ as_factor("small"),
                                     log_body_length_mm >= 0 & log_body_length_mm < 1 ~ as_factor("medium"),
                                     log_body_length_mm >= 1 ~ as_factor("large"))) 


for(body_length_cat_i in levels(therm_traits_sizecats$body_length_cat)){
  therm_traits_body_length <- therm_traits_sizecats |> 
    filter(body_length_cat == body_length_cat_i)
  possible_error <- tryCatch(expr = {ctmax_lat_body_length_i <-  lmer(ctmax ~ abs(lat) + (1|species),
                                                                     data = therm_traits_body_length)
  ctmax_to_lat_plot_body_length <- sim_and_plot_linears(ctmax_lat_body_length_i,
                                                       var_x = abs(therm_traits_body_length$lat),
                                                       var_y = therm_traits_body_length$ctmax,
                                                       n_sims = 1000,
                                                       your_title = expression(CT[max]),
                                                       your_subtitle = body_length_cat_i,
                                                       lab_x = "Absolute Latitude (º)",
                                                       lab_y = "ºC",
                                                       color_points = "#bb3e03",
                                                       color_central = "#9b2226",
                                                       color_uncertainty = "#ee9b00")
  },  # <- inside tryCatch
  error = function(e) e)
  if(inherits(possible_error, "error")) {
    ctmax_lat_body_length_i <- NULL
    ggplot(therm_traits_body_length, aes(x = abs(therm_traits_body_length$lat),
                                        y = therm_traits_body_length$ctmax))+
      geom_point(color = "#bb3e03")+
      geom_smooth(method = "lm",
                  fill = "#ee9b00",
                  color = "#9b2226")+
      ggthemes::theme_clean()+
      labs(x = "Absolute Latitude",
           y = "ºC",
           title = expression(CT[max]),
           subtitle = "No lmer fitting: using `geom_smooth` instead for exploratory analysis")
  }
  ggsave(filename = here(paste0("data/data_sink/figs/traits_lat/ctmax_lat_", body_length_cat_i,"_plot.png")),
         width = 16, height = 16, units = "cm")
}


###a.2. ctmin - lat ----
ctmin_lat_ref <- lmerTest::lmer(ctmin ~ abs(lat) + (1|reference),
                                data = therm_traits)
summary(ctmin_lat_ref) 


therm_traits_deutsch <- therm_traits |> filter(model_name == "deutsch_2008")
ctmin_lat_ref_deutsch <-  lmerTest::lmer(ctmin ~ abs(lat) + (1|reference),
                   data = therm_traits_deutsch)
summary(ctmin_lat_ref_deutsch)
#all models
ctmin_to_lat_plot <- sim_and_plot_linears(ctmin_lat,
                                          var_x = abs(therm_traits$lat),
                                          var_y = therm_traits$ctmin,
                                          n_sims = 1000,
                                          your_title = expression(CT[min]),
                                          your_subtitle = NULL,
                                          lab_x = "Absolute Latitude (º)",
                                          lab_y = "ºC",
                                          color_points = "#0a9396",
                                          color_central = "#005f73",
                                          color_uncertainty = "#94d2bd")
ggsave(ctmin_to_lat_plot,
       filename = here("data/data_sink/figs/traits_lat/ctmin_to_lat_plot.png"),
         width = 16, height = 16, units = "cm")

ctmin_to_lat_ref_plot <- sim_and_plot_linears(ctmin_lat_ref,
                                          var_x = abs(therm_traits$lat),
                                          var_y = therm_traits$ctmin,
                                          n_sims = 1000,
                                          your_title = "Critical Thermal Minima",
                                          your_subtitle = NULL,
                                          lab_x = "Absolute Latitude (º)",
                                          lab_y = "Temperature (ºC)",
                                          color_points = "#0a9396",
                                          color_central = "#005f73",
                                          color_uncertainty = "#94d2bd")
ggsave(ctmin_to_lat_ref_plot,
       filename = here("data/data_sink/figs/traits_lat/ctmin_to_lat_ref_plot.png"),
       width = 16, height = 16, units = "cm")

#only deutsch
ctmin_to_lat_ref_deutsch_plot <- sim_and_plot_linears(ctmin_lat_ref_deutsch,
                                                  var_x = abs(therm_traits_deutsch$lat),
                                                  var_y = therm_traits_deutsch$ctmin,
                                                  n_sims = 1000,
                                                  your_title = expression(CT[min]),
                                                  your_subtitle = "Only deutsch_2008 models",
                                                  lab_x = "Absolute Latitude (º)",
                                                  lab_y = "Temperature (ºC)",
                                                  color_points = "#0a9396",
                                                  color_central = "#005f73",
                                                  color_uncertainty = "#94d2bd")
ggsave(ctmin_to_lat_deutsch_plot,
       filename = here("data/data_sink/figs/traits_lat/ctmin_to_lat_deutsch_plot.png"),
       width = 16, height = 16, units = "cm")



traits_lat_plotgrid <- cowplot::plot_grid(ctmin_to_lat_plot,
                                          #ctmax_lat_deutsch_plot,
                                          ctmax_to_lat_plot,
                                          # ctmin_to_lat_deutsch_plot,
                                          nrow = 1,
                                          labels = "AUTO")
ggsave(traits_lat_plotgrid,
       filename = here("data/data_sink/figs/traits_lat/traits_lat_plotgrid.png"),
       width = 16, height = 16, units = "cm")

traits_lat_plotgrid_refs <- cowplot::plot_grid(ctmin_to_lat_ref_plot,
                                          #ctmax_lat_deutsch_plot,
                                          ctmax_to_lat_ref_plot,
                                          # ctmin_to_lat_deutsch_plot,
                                          nrow = 1,
                                          labels = "AUTO")
ggsave(traits_lat_plotgrid_refs,
       filename = here("data/data_sink/figs/traits_lat/traits_lat_plotgrid_refs.png"),
       width = 16, height = 16, units = "cm")

traits_lat_plotgrid_deutsch_ref <- cowplot::plot_grid(ctmax_to_lat_ref_deutsch_plot, #ctmax_lat_deutsch_plot,
                                                      ctmin_to_lat_ref_deutsch_plot,# ctmin_to_lat_deutsch_plot,
                                                      nrow = 1,
                                                      labels = "AUTO")
ggsave(traits_lat_plotgrid_deutsch_ref,
       filename = here("data/data_sink/figs/traits_lat/traits_lat_plotgrid_deutsch_ref.png"),
       width = 16, height = 16, units = "cm")




###### by order ----
for(order_i in levels(therm_traits$order)){
  therm_traits_order <- therm_traits |> 
    filter(order == order_i)
  possible_error <- tryCatch(expr = {ctmin_lat_order_i <-  lmer(ctmin ~ abs(lat) + (1|species),
                                                                data = therm_traits_order)
  ctmin_to_lat_plot_order <- sim_and_plot_linears(ctmin_lat_order_i,
                                                  var_x = abs(therm_traits_order$lat),
                                                  var_y = therm_traits_order$ctmin,
                                                  n_sims = 1000,
                                                  your_title = expression(CT[min]),
                                                  your_subtitle = order_i,
                                                  lab_x = "Absolute Latitude (º)",
                                                  lab_y = "ºC",
                                                  color_points = "#0a9396",
                                                  color_central = "#005f73",
                                                  color_uncertainty = "#94d2bd")
  },  # <- inside tryCatch
  error = function(e) e)
  if(inherits(possible_error, "error")) {
    ctmin_lat_order_i <- NULL
    ggplot(therm_traits_order, aes(x = abs(therm_traits_order$lat),
                                   y = therm_traits_order$ctmin))+
      geom_point(color = "#0a9396")+
      geom_smooth(method = "lm",
                  fill = "#94d2bd",
                  color = "#005f73")+
      ggthemes::theme_clean()+
      labs(x = "Absolute Latitude",
           y = "ºC",
           title = expression(CT[min]),
           subtitle = "No lmer fitting: using `geom_smooth` instead for exploratory analysis")
  }
  ggsave(filename = here(paste0("data/data_sink/figs/traits_lat/ctmin_lat_", order_i,"_plot.png")),
         width = 16, height = 16, units = "cm")
}

###### by feed_guild ----
for(feed_guild_i in levels(therm_traits$feed_guild)){
  therm_traits_feed_guild <- therm_traits |> 
    filter(feed_guild == feed_guild_i)
  possible_error <- tryCatch(expr = {ctmin_lat_feed_guild_i <-  lmer(ctmin ~ abs(lat) + (1|species),
                                                                     data = therm_traits_feed_guild)
  ctmin_to_lat_plot_feed_guild <- sim_and_plot_linears(ctmin_lat_feed_guild_i,
                                                       var_x = abs(therm_traits_feed_guild$lat),
                                                       var_y = therm_traits_feed_guild$ctmin,
                                                       n_sims = 1000,
                                                       your_title = expression(CT[min]),
                                                       your_subtitle = feed_guild_i,
                                                       lab_x = "Absolute Latitude (º)",
                                                       lab_y = "ºC",
                                                       color_points = "#0a9396",
                                                       color_central = "#005f73",
                                                       color_uncertainty = "#94d2bd")
  },  # <- inside tryCatch
  error = function(e) e)
  if(inherits(possible_error, "error")) {
    ctmin_lat_feed_guild_i <- NULL
    ggplot(therm_traits_feed_guild, aes(x = abs(therm_traits_feed_guild$lat),
                                        y = therm_traits_feed_guild$ctmin))+
      geom_point(color = "#0a9396")+
      geom_smooth(method = "lm",
                  fill = "#94d2bd",
                  color = "#005f73")+
      ggthemes::theme_clean()+
      labs(x = "Absolute Latitude",
           y = "ºC",
           title = expression(CT[min]),
           subtitle = "No lmer fitting: using `geom_smooth` instead for exploratory analysis")
  }
  ggsave(filename = here(paste0("data/data_sink/figs/traits_lat/ctmin_lat_", feed_guild_i,"_plot.png")),
         width = 16, height = 16, units = "cm")
}

###### by sizecats ----

therm_traits_sizecats <- therm_traits |> 
  mutate(log_body_length_mm = log10(body_length_mm),
         body_length_cat = case_when(log_body_length_mm < 0 ~ as_factor("small"),
                                     log_body_length_mm >= 0 & log_body_length_mm < 1 ~ as_factor("medium"),
                                     log_body_length_mm >= 1 ~ as_factor("large"))) 


for(body_length_cat_i in levels(therm_traits_sizecats$body_length_cat)){
  therm_traits_body_length <- therm_traits_sizecats |> 
    filter(body_length_cat == body_length_cat_i)
  possible_error <- tryCatch(expr = {ctmin_lat_body_length_i <-  lmer(ctmin ~ abs(lat) + (1|species),
                                                                      data = therm_traits_body_length)
  ctmin_to_lat_plot_body_length <- sim_and_plot_linears(ctmin_lat_body_length_i,
                                                        var_x = abs(therm_traits_body_length$lat),
                                                        var_y = therm_traits_body_length$ctmin,
                                                        n_sims = 1000,
                                                        your_title = expression(CT[min]),
                                                        your_subtitle = body_length_cat_i,
                                                        lab_x = "Absolute Latitude (º)",
                                                        lab_y = "ºC",
                                                        color_points = "#0a9396",
                                                        color_central = "#005f73",
                                                        color_uncertainty = "#94d2bd")
  },  # <- inside tryCatch
  error = function(e) e)
  if(inherits(possible_error, "error")) {
    ctmin_lat_body_length_i <- NULL
    ggplot(therm_traits_body_length, aes(x = abs(therm_traits_body_length$lat),
                                         y = therm_traits_body_length$ctmin))+
      geom_point(color = "#0a9396")+
      geom_smooth(method = "lm",
                  fill = "#94d2bd",
                  color = "#005f73")+
      ggthemes::theme_clean()+
      labs(x = "Absolute Latitude",
           y = "ºC",
           title = expression(CT[min]),
           subtitle = "No lmer fitting: using `geom_smooth` instead for exploratory analysis")
  }
  ggsave(filename = here(paste0("data/data_sink/figs/traits_lat/ctmin_lat_", body_length_cat_i,"_plot.png")),
         width = 16, height = 16, units = "cm")
}



###a.3. breadth - lat ----

breadth_lat_ref <- lmerTest::lmer(breadth ~ abs(lat) + (1|reference),
                                data = therm_traits)
summary(breadth_lat_ref) 


breadth_lat_plot <- sim_and_plot_linears(breadth_lat_ref,
                                         var_x = abs(therm_traits$lat),
                                         var_y = therm_traits$breadth,
                                         n_sims = 1000,
                                         your_title = "Thermal Breadth",
                                         your_subtitle = NULL,
                                         lab_x = "Absolute Latitude (º)",
                                         lab_y = "ºC",
                                         color_points = "#b56576",
                                         color_central = "#355070",
                                         color_uncertainty = "#eaac8b")
breadth_lat_plot
ggsave(filename = here("data/data_sink/figs/breadth_lat_plot.png"),
       width = 16, height = 16, units = "cm")

###a.4. thermal_tolerance - lat ----
thermal_tolerance_lat <-  lmerTest::lmer(thermal_tolerance ~ abs(lat) + (1|reference),
                               na.action = na.omit,
                               data = therm_traits)
summary(thermal_tolerance_lat)

thermal_tolerance_lat_nogauss <- lmerTest::lmer(thermal_tolerance ~abs(lat) + (1|reference),
                                                data = therm_traits_nogaussians)
summary(thermal_tolerance_lat_nogauss)

thermal_tolerance_lat_deutsch <- lmerTest::lmer(thermal_tolerance ~abs(lat) + (1|reference),
                                                data = therm_traits |> filter(model_name == "deutsch_2008"))
summary(thermal_tolerance_lat_deutsch)

thermal_tolerance_lat_plot <- sim_and_plot_linears(thermal_tolerance_lat,
                                         var_x = abs(therm_traits$lat),
                                         var_y = therm_traits$thermal_tolerance,
                                         n_sims = 1000,
                                         your_title = expression(T[tolerance]),
                                         your_subtitle = NULL,
                                         lab_x = "Absolute Latitude (º)",
                                         lab_y = "ºC",
                                         color_points = "#b56576",
                                         color_central = "#355070",
                                         color_uncertainty = "#eaac8b")
thermal_tolerance_lat_plot
ggsave(filename = here("data/data_sink/figs/thermal_tolerance_lat_plot.png"),
       width = 16, height = 16, units = "cm")

thermal_tolerance_lat_nogauss_plot <- sim_and_plot_linears(thermal_tolerance_lat,
                                                   var_x = abs(therm_traits_nogaussians$lat),
                                                   var_y = therm_traits_nogaussians$thermal_tolerance,
                                                   n_sims = 1000,
                                                   your_title = expression(T[tolerance]),
                                                   your_subtitle = "Gaussian models excluded",
                                                   lab_x = "Absolute Latitude (º)",
                                                   lab_y = "ºC",
                                                   color_points = "#b56576",
                                                   color_central = "#355070",
                                                   color_uncertainty = "#eaac8b")
thermal_tolerance_lat_nogauss_plot
ggsave(filename = here("data/data_sink/figs/thermal_tolerance_lat_nogauss_plot.png"),
       width = 16, height = 16, units = "cm")

##hemisphere partitioning
therm_traits_south <- therm_traits |> 
  filter(lat < 0)

therm_traits_north <- therm_traits |> 
  filter(lat > 0)

therm_traits_hemisphere <- therm_traits |> 
  mutate(hemisphere = case_when(lat < 0 ~ "Southern Hemisphere",
                                lat > 0 ~ "Northern Hemisphere"))

therm_traits_south_lmer <- lmerTest::lmer(thermal_tolerance ~ lat + (1|reference),
                                          na.action = na.omit,
                                          data = therm_traits_south |> filter(model_name != "gaussian_1987"))
summary(therm_traits_south_lmer)

therm_traits_north_lmer <- lmerTest::lmer(thermal_tolerance ~ lat + (1|reference),
                                          na.action = na.omit,
                                          data = therm_traits_north |> filter(model_name != "gaussian_1987"))
summary(therm_traits_north_lmer)



tolerance_lat_nogauss_north_plot<- sim_and_plot_linears(therm_traits_north_lmer,
                                                        var_x = therm_traits_north$lat,
                                                        var_y = therm_traits_north$thermal_tolerance,
                                                        n_sims = 1000,
                                                        your_title = expression(T[tolerance]),
                                                        your_subtitle = "Northern hemisphere; no gaussians",
                                                        lab_x = "º Latitude",
                                                        lab_y = "Temperature (ºC)",
                                                        color_points = "#5B6D83",
                                                        color_central = "#245C7D",
                                                        color_uncertainty = "#D28F86")

tolerance_lat_nogauss_south_plot<- sim_and_plot_linears(therm_traits_south_lmer,
                                                        var_x = therm_traits_south$lat,
                                                        var_y = therm_traits_south$thermal_tolerance,
                                                        n_sims = 1000,
                                                        your_title = expression(T[tolerance]),
                                                        your_subtitle = "Southern hemisphere; no gaussians",
                                                        lab_x = "º Latitude",
                                                        lab_y = "Temperature (ºC)",
                                                        color_points = "#5B6D83",
                                                        color_central = "#245C7D",
                                                        color_uncertainty = "#D28F86")

traits_nogauss_hemispheres <- cowplot::plot_grid(tolerance_lat_nogauss_south_plot, 
                                                 tolerance_lat_nogauss_north_plot)

ggsave(here("data/data_sink/figs/traits_nogauss_hemispheres.svg"),
       width = 8,
       height = 5)

###a.5. therm_range - lat ----
thermal_range_lat <-  lmerTest::lmer(therm_range ~ abs(lat) + (1|reference),
                               na.action = na.omit,
                               data = therm_traits)
summary(thermal_range_lat)
thermal_range_lat_plot <- sim_and_plot_linears(thermal_range_lat,
                                                   var_x = abs(therm_traits$lat),
                                                   var_y = therm_traits$therm_range,
                                                   n_sims = 1000,
                                                   your_title = expression(T[range]),
                                                   your_subtitle = NULL,
                                                   lab_x = "Absolute Latitude (º)",
                                                   lab_y = "ºC",
                                               color_points = "#b56576",
                                               color_central = "#355070",
                                               color_uncertainty = "#eaac8b")
thermal_range_lat_plot
ggsave(filename = here("data/data_sink/figs/thermal_range_lat_plot.png"),
       width = 16, height = 16, units = "cm")

###a.6. thermal_safety_margin - lat ----
thermal_safety_margin_lat <-  lmer(thermal_safety_margin ~ abs(lat) + (1|reference),
                                   na.action = na.omit,
                                   data = therm_traits)
summary(thermal_safety_margin_lat)
thermal_safety_margin_lat_plot <- sim_and_plot_linears(thermal_safety_margin_lat,
                                                       var_x = abs(therm_traits$lat),
                                                       var_y = therm_traits$thermal_safety_margin,
                                                       n_sims = 1000,
                                                       your_title = expression(T[safety~margin]),
                                                       your_subtitle = NULL,
                                                       lab_x = "Absolute Latitude (º)",
                                                       lab_y = "ºC",
                                                       color_points = "#b56576",
                                                       color_central = "#355070",
                                                       color_uncertainty = "#eaac8b")
thermal_safety_margin_lat_plot
ggsave(filename = here("data/data_sink/figs/thermal_safety_margin_lat_plot.png"),
       width = 16, height = 16, units = "cm")

###a.7. topt - lat ----
topt_lat <-  lmerTest::lmer(topt ~ abs(lat) + (1|species),
                             data = therm_traits)

topt_lat_ref <- lmerTest::lmer(topt ~ abs(lat) + (1|reference),
                                data = therm_traits)
summary(topt_lat)
summary(topt_lat_ref) 
topt_lat_ref_plot <- sim_and_plot_linears(topt_lat_ref,
                                          var_x = abs(therm_traits$lat),
                                          var_y = therm_traits$topt,
                                          n_sims = 1000,
                                          your_title = expression(T[opt]),
                                          your_subtitle = NULL,
                                          lab_x = "Absolute Latitude (º)",
                                          lab_y = "ºC",
                                          color_points = "#b56576",
                                          color_central = "#355070",
                                          color_uncertainty = "#eaac8b")
topt_lat_ref_plot
ggsave(filename = here("data/data_sink/figs/topt_lat_ref_plot.png"),
       width = 16, height = 16, units = "cm")

rmax_lat_ggplot <- ggplot(therm_traits, aes(x = abs(lat),
                                             y = log(rmax)))+
  geom_point(color = "#b56576")+
  geom_smooth(method = "lm", color = "#355070", fill = "#eaac8b")+
  labs(x = "Latitude (abs. º)",
       y = expression(italic(r)[max]))+
  theme_few()




###a.8. rmax - lat ----
rmax_lat <-  lmerTest::lmer(log(rmax) ~ abs(lat) + (1|species),
                            data = therm_traits)

rmax_lat_ref <- lmerTest::lmer(log(rmax) ~ abs(lat) + (1|reference),
                               data = therm_traits)
summary(rmax_lat)
summary(rmax_lat_ref) 

rmax_lat_ref_body <- lmerTest::lmer(log(rmax) ~ abs(lat) + log10(body_length_mm) + (1|reference),
                               data = therm_traits)
summary(rmax_lat_ref_body)


rmax_lat_plot <- sim_and_plot_linears(rmax_lat_ref_body,
                                      var_x = abs(therm_traits$lat),
                                      var_y = log(therm_traits$rmax),
                                      n_sims = 1000,
                                      your_title = expression(italic(r)[max]),
                                      your_subtitle = NULL,
                                      lab_x = "º Absolute Latitude ",
                                      lab_y = "Maximum performance (logarithm)",
                                      color_points = "#b56576",
                                      color_central = "#355070",
                                      color_uncertainty = "#eaac8b")
rmax_lat_plot
ggsave(filename = here("data/data_sink/figs/rmax_lat_plot.png"),
       width = 16, height = 16, units = "cm")

### a.9. Activation Energy - lat ----

### a) lmer
hist(log(therm_traits$e, 30))

e_lat_ref <- lmerTest::lmer(log(e) ~ abs(lat) + (1|reference),
                               data = therm_traits)
summary(e_lat_ref)

therm_traits_deutsch <- therm_traits |> 
  filter(model_name == "deutsch_2008")

e_lat_ref_nogauss <- lmerTest::lmer(log(e) ~ abs(lat) + (1|reference),
                            data = therm_traits_nogaussians)
summary(e_lat_ref_nogauss)

# try it without selection using only pawar and ssi
therm_traits_mechanistic <- readRDS(here("data/data_sink/therm_traits_mechanistic.rds"))

mechanistic_traits_selected <- therm_traits |> 
  filter(model_name %in% c("rezende_2019", "pawar_2018", "sharpeschoolfull_1981"))

mechanistic_non_selected <- therm_traits_mechanistic |> 
  group_by(reference, species, model_name) |> 
  mutate(model_fit = map(.x = model_fit, 
                         .f = ~model_fit[[1]]),
         model_AIC = map_dbl(.x = model_fit,
                             .f = ~AIC(.x))
  ) |> 
  ungroup(model_name) |> 
  arrange(model_AIC) |> 
  mutate(min_AIC_diff = model_AIC - lag(model_AIC),
         min_AIC_diff = ifelse(is.na(min_AIC_diff), 
                               Inf, 
                               min_AIC_diff))  |> 
  filter(min_AIC_diff > 2) |> 
  slice_head(n = 1) 

e_lat_ref_mechan <- lmerTest::lmer(log(e) ~ abs(lat) + (1|reference),
                                   data = mechanistic_non_selected)
summary(e_lat_ref_mechan) # no variation

## b) cold-adaptation hypothesis --------------------------------------------
### b.1. Variability ~ order ----
therm_traits_deutsch <- therm_traits |> 
  filter(model_name == "deutsch_2008")
trait_clouds_order <- ggplot(therm_traits_deutsch, aes(x = order)) + 
  ggdist::stat_halfeye(aes(y = ctmin),
                       adjust = .5,
                       width = .5, ## set slab interval to show IQR and 95% data range
                       color = "#0081a7",
                       fill = "#00afb9"
  ) + 
  ggdist::stat_dots(aes(y = ctmin),
                    side = "left",
                    dotsize = 2,
                    justification = 1.05,
                    binwidth = .5,
                    color = "#0081a7",
                    fill= "#0081a7",
                    alpha = .5
  )+
  coord_cartesian(xlim = c(1, NA))+
  ggdist::stat_halfeye(aes(y = ctmax),
                       adjust = .5,
                       width = .5, ## set slab interval to show IQR and 95% data range
                       color = "#9b2226",
                       fill = "#f07167"
  ) + 
  ggdist::stat_dots(aes(y = ctmax),
                    side = "left",
                    dotsize = 2,
                    justification = 1.05,
                    binwidth = .5,
                    color = "#9b2226",
                    fill= "#9b2226",
                    alpha = .5
  )+
  theme_few()+
  labs(x = "Order",
       y = "Thermal traits (ºC)")
trait_clouds_order
ggsave(filename = here("data/data_sink/figs/trait_clouds_order.png"),
       width = 16, height = 16, units = "cm")

# how about excluding gaussians, since they may overestimate the right part
# of the curve
trait_clouds_order_nogauss <- ggplot(therm_traits |> filter(model_name != "gaussian_1987"), aes(x = order)) + 
  ggdist::stat_halfeye(aes(y = ctmin),
                       adjust = .5,
                       width = .5, ## set slab interval to show IQR and 95% data range
                       color = "#0081a7",
                       fill = "#00afb9"
  ) + 
  ggdist::stat_dots(aes(y = ctmin),
                    side = "left",
                    dotsize = 2,
                    justification = 1.05,
                    binwidth = .5,
                    color = "#0081a7",
                    fill= "#0081a7",
                    alpha = .5
  )+
  coord_cartesian(xlim = c(1, NA))+
  ggdist::stat_halfeye(aes(y = ctmax),
                       adjust = .5,
                       width = .5, ## set slab interval to show IQR and 95% data range
                       color = "#9b2226",
                       fill = "#f07167"
  ) + 
  ggdist::stat_dots(aes(y = ctmax),
                    side = "left",
                    dotsize = 2,
                    justification = 1.05,
                    binwidth = .5,
                    color = "#9b2226",
                    fill= "#9b2226",
                    alpha = .5
  )+
  theme_few()+
  labs(title = "Thermal traits distribution",
       subtitle = "Gaussian curves excluded",
       x = "Order",
       y = "Temperature (ºC)")
trait_clouds_order_nogauss
ggsave(filename = here("data/data_sink/figs/trait_clouds_order_deutsch.png"),
       width = 16, height = 16, units = "cm")

trait_clouds_deutsch <- ggplot(therm_traits |> filter(model_name != "gaussian_1987"), aes(x = 1.5)) + 
  ggdist::stat_halfeye(aes(y = ctmin),
                       adjust = .5,
                       width = .5, ## set slab interval to show IQR and 95% data range
                       color = "#0081a7",
                       fill = "#00afb9"
  ) + 
  ggdist::stat_dots(aes(y = ctmin),
                    side = "left",
                    dotsize = 1.2,
                    justification = 1.05,
                    binwidth = .5,
                    color = "#0081a7",
                    fill= "#0081a7",
                    alpha = .5
  )+
  coord_cartesian(xlim = c(1, NA))+
  ggdist::stat_halfeye(aes(y = ctmax),
                       adjust = .5,
                       width = .5, ## set slab interval to show IQR and 95% data range
                       color = "#9b2226",
                       fill = "#f07167"
  ) + 
  ggdist::stat_dots(aes(y = ctmax),
                    side = "left",
                    dotsize = 1.2,
                    justification = 1.05,
                    binwidth = .5,
                    color = "#9b2226",
                    fill= "#9b2226",
                    alpha = .5
  )+
  theme_clean()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  labs(x = element_blank(),
       y = "Thermal traits (ºC)")
trait_clouds_deutsch
ggsave(filename = here("data/data_sink/figs/trait_clouds_nogauss.png"),
       width = 11, height = 11, units = "cm")
## and using only deutsch
trait_clouds_deutsch <- ggplot(therm_traits |> filter(model_name == "deutsch_2008"), aes(x = 1.5)) + 
  ggdist::stat_halfeye(aes(y = ctmin),
                       adjust = .5,
                       width = .5, ## set slab interval to show IQR and 95% data range
                       color = "#0081a7",
                       fill = "#00afb9"
  ) + 
  ggdist::stat_dots(aes(y = ctmin),
                    side = "left",
                    dotsize = 1.2,
                    justification = 1.05,
                    binwidth = .5,
                    color = "#0081a7",
                    fill= "#0081a7",
                    alpha = .5
  )+
  coord_cartesian(xlim = c(1, NA))+
  ggdist::stat_halfeye(aes(y = ctmax),
                       adjust = .5,
                       width = .5, ## set slab interval to show IQR and 95% data range
                       color = "#9b2226",
                       fill = "#f07167"
  ) + 
  ggdist::stat_dots(aes(y = ctmax),
                    side = "left",
                    dotsize = 1.2,
                    justification = 1.05,
                    binwidth = .5,
                    color = "#9b2226",
                    fill= "#9b2226",
                    alpha = .5
  )+
  theme_clean()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  labs(x = element_blank(),
       y = "Thermal traits (ºC)")
trait_clouds_deutsch
ggsave(filename = here("data/data_sink/figs/trait_clouds_deutsch.png"),
       width = 11, height = 11, units = "cm")
### b.2. F-test (see Herrando-Perez, 2019) ----
variances_ctmin <- therm_traits |> 
  select(ctmin, vi, reference) |> 
  mutate(parameter = rep("ctmin")) |> 
  rename(estimate = ctmin)
variances_ctmax <- therm_traits |> 
  select(ctmax, vi, reference) |> 
  mutate(parameter = rep("ctmax")) |> 
  rename(estimate = ctmax)
variances_traits = variances_ctmin |> 
  bind_rows(variances_ctmax) |> 
  glimpse()
#test variance differences
vartest_estimates <- var.test(estimate ~ parameter,
                              variances_traits, 
                              alternative = "greater") 
print(vartest_estimates) #CTmin estimate variance is significantly greater than CTmax estimate variance
### (2) Levene's test (see Hoffmann, 2013) ----
car::leveneTest(y = estimate ~ parameter,
                data = variances_traits) 
# their variance differ 

## c) hotter is better --------------------------------------------

#  a visual examination may help guide the analysis. 
#  See `data/data_sink/figs/all_tpcs_fitted.png`
#  the pattern seems to slightly emerge amidst the noise (as temperature rises, the peaks of the 
#  curves rise too)

# let's test this with modelling

###### c.1. log(rmax) ~ topt ----

## log(rmax) ~topt
logrmax_topt <- lmerTest::lmer(log(rmax) ~ topt + (1|reference),
                     data = therm_traits,
                     na.action = na.omit)
summary(logrmax_topt)

# by TPC models
logrmax_topt_deutsch <- lmerTest::lmer(log(rmax) ~ topt + (1|reference),
                               data = therm_traits |> 
                                 filter(model_name == "deutsch_2008"),
                               na.action = na.omit)
summary(logrmax_topt_deutsch)

logrmax_topt_nogauss <- lmerTest::lmer(log(rmax) ~ topt + (1|reference),
                                        data = therm_traits |> 
                                          filter(model_name != "gaussian_1987"),
                                        na.action = na.omit)
summary(logrmax_topt_nogauss)



#by families 
therm_traits |> count(family) |> arrange(-n)

therm_traits_family <- therm_traits |> 
  mutate(log_rmax = log(rmax)) |> 
  filter(family %in% c("Aphididae", "Tetranychidae", "Noctuidae", "Tephritidae", 
                       "Pyralidae", "Trichogrammatidae", "Thripidae", "Pseudococcidae",
                       "Crambidae", "Gelechiidae", "Liposcelididae", "Agromyzidae")) |> 
  mutate(cat_bs = case_when(log10(body_length_mm) < 0 ~ "small",
                            log10(body_length_mm) >= 0 & 
                              log10(body_length_mm) < 1 ~ "medium",
                            log10(body_length_mm) >=1 ~ "large")
        )
ggplot(therm_traits_family, aes(x = topt, y = log_rmax,
                                 color = family, fill = family))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(.~cat_bs)+
  theme_few()
ggsave(filename = here("data/data_sink/figs/therm_traits_family.png"),
       width = 16, height = 16, units = "cm")

#by orders
ggplot(therm_traits, aes(x = topt, y = log_rmax,
                         color = order, fill = order))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(.~order)+
  theme_few()
ggsave(filename = here("data/data_sink/figs/therm_traits_order.png"),
       width = 16, height = 16, units = "cm")

#see hemiptera's families
therm_traits_hemiptera <- therm_traits |> filter(order == "Hemiptera") #|> count(family)
#aphids are the only family with a different pattern.
range(therm_traits_hemiptera |> filter(family == "Aphididae") |> pull(body_length_mm))


rmax_topt_plot <- sim_and_plot_linears(logrmax_topt,
                                      var_x = therm_traits$topt,
                                      var_y = log(therm_traits$rmax),
                                      n_sims = 1000,
                                      your_title = NULL,
                                      your_subtitle = NULL,
                                      lab_x = expression(T[opt]~(ºC)),
                                      lab_y = expression(ln(italic(r)[max])),
                                      color_points = "#ce6a85",
                                      color_central = "#5c374c",
                                      color_uncertainty = "#faa275")+
  xlim(10, 40)+
  ylim(-5, 2.5)
rmax_topt_plot
ggsave(filename = here("data/data_sink/figs/rmax_topt_plot.png"),
       width = 12, height = 12, units = "cm")

rmax_topt_ggplot <- ggplot(therm_traits, aes(x = topt,
                                             y = log(rmax)))+
  geom_point(color = "#ce6a85")+
  geom_smooth(method = "lm", color = "#5c374c", fill = "#faa275")+
  labs(x =expression(T[opt]~(ºC)),
       y = expression(italic(r)[max]))+
  theme_few()
ggsave(filename = here("data/data_sink/figs/rmax_topt_ggplot.png"),
       width = 16, height = 16, units = "cm")

###### c.2. log(rmax) body size ----

## plot by categorical sizes
therm_traits_sizecats <- therm_traits |> 
  mutate(size_cat = case_when(log10(body_length_mm) < 0 ~ "small",
                              log10(body_length_mm) >= 0 &
                                log10(body_length_mm) < 1 ~ "medium",
                              log10(body_length_mm) >= 1 ~ "large")
         )

hotter_better_sizes_topt<- ggplot(therm_traits_sizecats,
                                  aes(x = topt, 
                                      y = log(rmax),
                                      fill = size_cat,
                                      color = size_cat))+
  geom_point(aes(shape = size_cat,
                 color = size_cat), alpha = .45)+
  geom_smooth(method = "lm")+
  scale_fill_manual(values = RColorBrewer::brewer.pal(3,"Dark2"))+
  scale_color_manual(values = RColorBrewer::brewer.pal(3,"Dark2"))+
  facet_wrap(.~size_cat)+
  ggthemes::theme_clean()+
  theme(legend.position = "none")+
  labs(x = "Topt",
       y = "Maximum performance")
hotter_better_sizes_topt # only evident for small organisms the HiB pattern

ggsave(filename = here("data/data_sink/figs/hotter_better_sizecats2.png"),
       width = 16, height = 16, units = "cm")




#correct in model:         
rmax_topt_bscor1 <- lmerTest::lmer(log(rmax)/body_length_mm ~ topt + (1|reference),
                                  data = therm_traits) ## nonsignificant
summary(rmax_topt_bscor1)
rmax_topt_bscor_plot <- sim_and_plot_linears(rmax_topt_bscor,
                                       var_x = therm_traits$topt,
                                       var_y = therm_traits$rmax_bsize,
                                       n_sims = 1000,
                                       your_title = NULL,
                                       your_subtitle = "Weights: Body length (mm)",
                                       lab_x = expression(T[opt]~(ºC)),
                                       lab_y = expression(log(italic(r)[max])),
                                       color_points = "#ce6a85",
                                       color_central = "#5c374c",
                                       color_uncertainty = "#faa275")
rmax_topt_bscor_plot

rmax_topt_bscor2 <- lmerTest::lmer(log(rmax) ~ topt + log10(body_length_mm) + (1|reference),
                                   data = therm_traits) ## significant
summary(rmax_topt_bscor2)

## latitude?
rmax_topt_bscor2_lat <- lmerTest::lmer(log(rmax) ~ topt + log10(body_length_mm) + abs(lat) + topt*abs(lat) + (1|reference),
                                   data = therm_traits) ## significant
summary(rmax_topt_bscor2_lat) # significant for topt/abs only for

rmax_topt_bscor_plot <- sim_and_plot_linears(rmax_topt_bscor2,
                                             var_x = therm_traits$topt,
                                             var_y = log(therm_traits$rmax),
                                             n_sims = 1000,
                                             your_title = NULL,
                                             your_subtitle = NULL,
                                             lab_x = expression(T[opt]~(ºC)),
                                             lab_y = expression(log(italic(r)[max])),
                                             color_points = "#ce6a85",
                                             color_central = "#5c374c",
                                             color_uncertainty = "#faa275")
ggsave(filename = here("data/data_sink/figs/hotter_better_sizecorrected.png"),
       width = 16, height = 16, units = "cm")

rmax_topt_bscor_models <- lmerTest::lmer(log(rmax) ~ topt + log10(body_length_mm) + (1|reference),
                                   data = therm_traits) ## significant

ggplot(therm_traits, aes(x = topt, y = log(rmax)))+
  geom_point(color = "#ce6a85")+
  geom_smooth(method = "lm", color = "#5c374c", fill = "#faa275")+
  labs(x = expression(T[opt]~(ºC)),
       y = expression(ln~(italic(r)[max])))+
  facet_wrap(~model_name)+
  theme_few()
ggsave(here("data/data_sink/figs/hib_model.png"), width = 15, height = 15, units = "cm")





## d) generalist-specialist trade-off --------------------------------------------
# we expect that specialists will have higher r_max at a cost of narrower thermal breadths
# thus, if the hypothesis is true, r_max should be negatively correlated to thermal breadth

### d.1. all models --------------------------------------------
logrmax_bs_breadth <- lmerTest::lmer(log(rmax) ~ breadth + (1|reference),
                           data = therm_traits)
summary(logrmax_bs_breadth)
logrmax_bs_breadth_plot <- sim_and_plot_linears(logrmax_bs_breadth,
                                                 var_x = therm_traits$breadth,
                                                 var_y = log(therm_traits$rmax),
                                                 n_sims = 1000,
                                                 your_title = NULL,
                                                 your_subtitle = NULL,
                                                 lab_x = expression(italic(T)[breadth~Q[80]]~(ºC)),
                                                 lab_y = expression(ln(italic(r)[max])),
                                                 color_points = "#822faf",
                                                 color_central = "#47126b",
                                                 color_uncertainty = "#d55d92")

logrmax_bs_therm_range <- lmerTest::lmer(log(rmax) ~ therm_range + (1|reference),
                                     data = therm_traits)
summary(logrmax_bs_therm_range)
logrmax_bs_therm_range_plot <- sim_and_plot_linears(logrmax_bs_therm_range,
                                                var_x = therm_traits$therm_range,
                                                var_y = log(therm_traits$rmax),
                                                n_sims = 1000,
                                                your_title = NULL,
                                                your_subtitle = NULL,
                                                lab_x = expression(italic(T)[range]~(ºC)),
                                                lab_y = NULL,
                                                color_points = "#b56576",
                                                color_central = "#355070",
                                                color_uncertainty = "#eaac8b")
logrmax_bs_thermal_tolerance <- lmerTest::lmer(log(rmax) ~ thermal_tolerance + (1|reference),
                                         data = therm_traits)
summary(logrmax_bs_thermal_tolerance)
logrmax_bs_thermal_tolerance_plot <- sim_and_plot_linears(logrmax_bs_thermal_tolerance,
                                                    var_x = therm_traits$thermal_tolerance,
                                                    var_y = log(therm_traits$rmax),
                                                    n_sims = 1000,
                                                    your_title = NULL,
                                                    your_subtitle = NULL,
                                                    lab_x = expression(italic(T)[tolerance]~(ºC)),
                                                    lab_y = NULL,
                                                    color_points = "#e26d5c",
                                                    color_uncertainty = "#ffe1a8")


grid_spec_gen <- cowplot::plot_grid(logrmax_bs_breadth_plot,
                                    logrmax_bs_therm_range_plot,
                                    logrmax_bs_thermal_tolerance_plot,
                                    labels = "AUTO",
                                    nrow = 1)
grid_spec_gen
ggsave(filename = here("data/data_sink/figs/grid_spec_gen.png"),
       width = 20, height = 15, units = "cm")
### d.2. no gauss --------------------------------------------

therm_traits_nogauss <- therm_traits |> 
  filter(model_name != "gaussian_1987")

logrmax_bs_breadth_nogauss <- lmerTest::lmer(log(rmax) ~ breadth + (1|reference),
                                     data = therm_traits_nogauss)
summary(logrmax_bs_breadth_nogauss)
logrmax_bs_breadth_nogauss_plot <- sim_and_plot_linears(logrmax_bs_breadth_nogauss,
                                                var_x = therm_traits_nogauss$breadth,
                                                var_y = log(therm_traits_nogauss$rmax),
                                                n_sims = 1000,
                                                your_title = NULL,
                                                your_subtitle = NULL,
                                                lab_x = expression(T[breadth~Q[80]]~(ºC)),
                                                lab_y = expression(log(italic(r)[max])),
                                                color_points = "#822faf",
                                                color_central = "#47126b",
                                                color_uncertainty = "#d55d92")

logrmax_bs_therm_range_nogauss <- lmerTest::lmer(log(rmax) ~ therm_range + (1|reference),
                                         data = therm_traits_nogauss)
summary(logrmax_bs_therm_range_nogauss)
logrmax_bs_therm_range_nogauss_plot <- sim_and_plot_linears(logrmax_bs_therm_range_nogauss,
                                                    var_x = therm_traits_nogauss$therm_range,
                                                    var_y = log(therm_traits_nogauss$rmax),
                                                    n_sims = 1000,
                                                    your_title = NULL,
                                                    your_subtitle = NULL,
                                                    lab_x = expression(T[range]~(ºC)),
                                                    lab_y = NULL,
                                                    color_points = "#b56576",
                                                    color_central = "#355070",
                                                    color_uncertainty = "#eaac8b")

logrmax_bs_thermal_tolerance_nogauss <- lmerTest::lmer(log(rmax) ~ thermal_tolerance + (1|reference),
                                               data = therm_traits_nogauss)
summary(logrmax_bs_thermal_tolerance_nogauss)
logrmax_bs_thermal_tolerance_nogauss_plot <- sim_and_plot_linears(logrmax_bs_thermal_tolerance_nogauss,
                                                          var_x = therm_traits_nogauss$thermal_tolerance,
                                                          var_y = log(therm_traits_nogauss$rmax),
                                                          n_sims = 1000,
                                                          your_title = NULL,
                                                          your_subtitle = NULL,
                                                          lab_x = expression(T[tolerance]~(ºC)),
                                                          lab_y = NULL,
                                                          color_points = "#e26d5c",
                                                          color_central = "#472d30",
                                                          color_uncertainty = "#ffe1a8")


grid_spec_gen_nogauss <- cowplot::plot_grid(logrmax_bs_breadth_nogauss_plot,
                                    logrmax_bs_therm_range_nogauss_plot,
                                    logrmax_bs_thermal_tolerance_nogauss_plot,
                                    labels = "AUTO",
                                    nrow = 1)

ggsave(filename = here("data/data_sink/figs/grid_spec_gen_nogauss.png"),
       width = 20, height = 15, units = "cm")

### d.3. thermal_tolerane ~ topt --------------------------------------------
tolerance_topt_lmer_nogauss <- lmerTest::lmer(thermal_tolerance ~ topt + (1|reference),
                                      data = therm_traits_nogaussians)
summary(tolerance_topt_lmer_nogauss)



