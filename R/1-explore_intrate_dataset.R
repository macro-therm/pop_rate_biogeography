# Aim of the script: load `int_rate_pest`data set and explore it


# 0. Packages and data loading --------------------------------------------

library(tidyverse)
library(here)
library(readxl)
library(ape)

# setwd("~/GitHub/int_rate_pest_tempcopy" ) <- your directory path goes here

## load the data set
int_rate_pest <- read_delim("data/data_source/int_rate_pest_dataset.csv") |> 
  as_tibble() |> 
  mutate(source_dataset = as_factor(source_dataset)) |> 
  group_by(reference, species) |> 
  mutate(id = cur_group_id())

# 1. Explore variables ----------------------------------------------------

####  a) Intrinsic rate of increase --------------------------------------

## overview
summary(int_rate_pest$int_rate)

## histograms
ggplot(data = int_rate_pest, aes(x = int_rate))+
  geom_histogram(fill = "cyan4", binwidth = .01)+
  scale_fill_brewer(palette = 2)+
  labs(x = expression(intrinsic~rate~of~increase~(italic(r)[m])),
       y = "Frequency")+
  geom_vline(xintercept = 0.0,
             color = "gray23", 
             linetype = "dashed")+
  ggthemes::theme_few()

## density plots
ggplot(data = int_rate_pest, aes(x = int_rate))+
  ggdist::stat_dist_halfeye(fill = "cyan4", color = "gray78")+
  scale_fill_brewer(palette = 2)+
  labs(x = expression(intrinsic~rate~of~increase~(italic(r)[m])),
       y = "Frequency")+
  geom_vline(xintercept = 0.0,
             color = "gray23", 
             linetype = "dashed")+
  ggthemes::theme_few()
ggsave(here("data/data_sink/figs/int_rate_distribution.png"),
       height = 13,
       width = 13,
       units = "cm")

## calculate how many studies (in %) reported standard errors and sample sizes:
count_r_se <- int_rate_pest |>
  select(int_rate_se, reference) |> 
  drop_na() |> 
  distinct(reference) |> 
  pull(reference) |> 
  length()

count_total <- int_rate_pest |> 
  ungroup() |> 
  distinct(reference) |> 
  pull(reference) |> 
  length()

count_r_se_n <- int_rate_pest |> 
  select(int_rate_se, reference, sample_size) |> 
  drop_na() |> 
  distinct(reference) |> 
  pull(reference) |> 
  length()

#paste(round(count_r_se*100/count_total, 2), "% of studies", "(i.e.", count_r_se, "), reported","r_m", "Standard Error")
#paste(round(count_r_se_n*100/count_total, 2), "% of studies", "(i.e.", count_r_se_n, "), reported both SE and Sample size, being suitable for meta-analysis")


####  b) Temperatures --------------------------------------
## overview
summary(int_rate_pest$temperature)

ggplot(data = int_rate_pest, aes(x = temperature))+
  geom_histogram(fill = "firebrick3", color = "gray78")+
  scale_fill_brewer(palette = 2)+
  labs(x = expression(intrinsic~rate~of~increase~(italic(r)[m])),
       y = "Frequency")+
  ggthemes::theme_few()

## avoid mujica for this characterisation due to different methodology
temperature_treatments <- int_rate_pest |> 
  group_by(reference, species) |> 
  summarise(n_temps = n_distinct(temperature)) |> 
  filter(reference != "mujica2016") # <- ILCyM simulations
summary(temperature_treatments)

ggplot(data = temperature_treatments, aes(x = n_temps))+
  geom_histogram(fill = "goldenrod3", color = "gray23", binwidth = 1)+
  ggthemes::theme_few()+
  labs(x = "Temperature (ºC)",
       y = "Frequency")
ggsave(here("data/data_sink/figs/temperature_treatments.png"),
       height = 13,
       width = 13,
       units = "cm")

####  c) Taxonomy --------------------------------------
# obtain distinct classifications

distinct_studies <- int_rate_pest |> 
  group_by(reference, species, family, order) |> 
  distinct(reference) 
distinct_studies_order <- distinct_studies |> 
  ungroup() |> 
  count(order)

distinct_species_order <- distinct_studies |> 
  ungroup() |> 
  group_by(order) |>
  summarise(species = unique(species)) |> 
  count()

ggplot(data = distinct_studies_order, aes(x = fct_reorder(order, n),
                                          y = n))+
  geom_bar(aes(fill = order), stat = "identity")+
  scale_fill_manual(values = MoMAColors::MoMAPalettes$Klein[[1]])+
  theme_bw()+
  theme(legend.position = "none")+
  scale_y_continuous(position = "right")+
  coord_flip()+
  labs(x = "N", y = "Order")
ggsave(here("data/data_sink/figs/distinct_studies_order.png"), width = 15, height = 15, units = "cm")


##### i. temperature distribution (see Pawar et al. 2016) ----
boltzmann_ct <- 8.617e-05 # <- constant

# calculate thermal range, thermal range location and thermal range spread
t_dist <- int_rate_pest |>
  group_by(id, reference, species) |> 
  summarise(temperature = unique(temperature)) |> 
  ungroup(reference, species) |> 
  mutate(temp_kelvin = temperature + 273.15,
         temp_mean = mean(temp_kelvin),
         range = max(temp_kelvin) - min(temp_kelvin),
         range_location = 1/(boltzmann_ct*temp_mean)) |> 
  ungroup() |> 
  mutate(range_spread_i = (temp_kelvin - temp_mean)^2) |> 
  group_by(id) |> 
  mutate(range_spread = sum(range_spread_i))

## histograms
range_hist <- ggplot(data = t_dist, aes(x = range))+
  geom_histogram(fill = "#E9C86B", color = "gray")+
  labs(y = NULL)+
  ggthemes::theme_few()

range_location_hist <- ggplot(data = t_dist, aes(x = range_location))+
  geom_histogram(fill = "#F9B69C", color = "gray")+
  labs(y = NULL)+
  ggthemes::theme_few()

temp_hist <- ggplot(data = t_dist, aes(x = temperature))+
  geom_histogram(fill = "#387D89", color = "gray")+
  labs(y = NULL)+
  ggthemes::theme_few()

range_spread_hist <- ggplot(data = t_dist, aes(x = range_spread))+
  geom_histogram(fill = "#E75139", color = "gray")+
  labs(y = NULL)+
  ggthemes::theme_few()

cowplot::plot_grid(range_hist, range_location_hist, range_spread_hist, temp_hist)
ggsave(here("data/data_sink/figs/temp_distribution.png"), width = 15, height = 15, units = "cm")

####  d) Ecological feeding guilds --------------------------------------

## distinct classifications of ecological feeding guilds
distinct_studies_feed_guild <- int_rate_pest |> 
  group_by(reference, feed_guild) |> 
  distinct(reference) 
count_feed_guild <- distinct_studies_feed_guild |> 
  ungroup() |> 
  count(feed_guild)

ggplot(data = count_feed_guild, aes(x = fct_reorder(feed_guild, n),
                                          y = n))+
  geom_bar(aes(fill = feed_guild), stat = "identity")+
  theme_bw()+
  scale_fill_manual(values = MoMAColors::MoMAPalettes$Klein[[1]])+
  theme(legend.position = "none")+
  scale_y_continuous(position = "right")+
  coord_flip()+
  labs(y = "N", x = "Feeding guild")
ggsave(here("data/data_sink/figs/count_feed_guild.png"), width = 15, height = 15, units = "cm")

