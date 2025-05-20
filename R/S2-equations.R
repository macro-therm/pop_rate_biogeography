### Model equations to fit ----

## Since some models which are more complex (e.g., SSI) require more data treatments 

## Deutsch et al. (2008) ----
# also in rTPC
deutsch2008 <- function(temp, t_opt, sigma_p, ct_max){
  est <- ifelse(temp <= t_opt,
                exp(-((temp - t_opt) / (2 * sigma_p ))^2),
                1 - ((temp - t_opt) / (t_opt - ct_max))^2)
  return(est)
}

## Rezende et al. (2019) ----

rezende2019 <- function(temp, q_10, c, d, t_th){
  est <- if (temp < t_th) {
    (c * exp(temp * log(q_10) / 10))
  } else if (temp >= t_th) {
    (c * exp(temp * log(q_10) / 10)) * (1 - d * (temp - t_th)^2)
  } else {NA}
  return(est)
}

# also in rTPC::rezende_2019()

## Angilletta (2006) ----

#inherited from Angilletta 2006
# rTPC::modifiedgaussian_2006()

## Lynch (1987) ----
# rTPC::gaussian_1987()

##Ratkowsky 2005 ----
# rTPC::ratkowsky_2005()


##more from rTPC: delong, pawar, sharpeschoolhich, spain, flinn

## workflow:
 # 1. Write a function to fit all these models (variation from mappestRisk)
 # 2. Fit all these models to each study (no simulations)
 # 3. Meta-analytic models (lme and metafor)

seq_temps <- seq(0, 50, by = 0.01)
seq_rezende <- purrr::map_dbl(.x = seq_temps,
                              .f = ~rezende2019(temp = .x,
                                                c = 0.008,
                                                q_10 = 3.789,
                                                d = 0.0051,
                                                t_th = 20.96))
df_rezende_ex <- tibble(temperature = seq_temps,
                        preds = seq_rezende) |> 
  filter(preds >= 0)
ggplot(data = df_rezende_ex, aes(x = temperature, y = preds))+
  geom_line()

example_rezende <- rezende2019(temp =)












seq_temps <- seq(0, 50, by = 0.01)
seq_deutsch <- purrr::map_dbl(.x = seq_temps,
                              .f = ~deutsch2008(temp = .x,
                                                t_opt = 22.82,
                                                sigma_p = 3.6731,
                                                ct_max = 27.74))
df_deutsch_ex <- tibble(temperature = seq_temps,
                        preds = seq_deutsch) |> 
  filter(preds >= 0)
ggplot(data = df_deutsch_ex, aes(x = temperature, y = preds))+
  geom_line()
