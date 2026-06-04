library(tidyverse)
library(lme4)

if (startsWith(getwd(), "/home/lakrids")) {
  path_prefix <- "/home/lakrids/GenomeDK"
} else {path_prefix <- "/faststorage/project"}

df <- read_csv(paste0(path_prefix, "/megaFauna/sa_megafauna/data/bioclim/modeling_data.csv"))

df <- df %>% rename(temperature = mean_temp, npp = mean_npp, precipitation = mean_prec) %>%
  arrange(species, time_kya) %>%
  group_by(species) %>%
  mutate(lag_temperature = lag(temperature), temperature_s = scale(temperature)) %>%
  ungroup() 

df <- df %>%
  mutate(window_size = window_end - time_kya)
df %>% filter(startsWith(species, "Panthera"))

df %>%
  ggplot(aes(x = log10(time_kya), y = log(Ne), colour = species)) +
  geom_step() +
  labs(x = "Time (kya)", y = "Ne", colour = "Species")

#### Core mixed-effects model: Ne ~ combination of different climate + human predictors

# Here are some example, but of course more combinations are possible:
# linear temperature effect
m_linT <- lmer(
  log(Ne) ~ temperature_s + p_human +
    (1 + temperature_s + p_human | species),
  data = df
)
summary(m_linT)
coef(m_linT)$species
m_simple <- lmer(log(Ne) ~ temperature + p_human + (1 | species), 
                 data = df, weights = window_size)
AIC(m_linT, m_simple)

# temperature and precipitation effect
m_full <- lmer(
  log(Ne) ~ temperature + precipitation + p_human +
    (1 + temperature + precipitation + p_human | species),
  data = df
)

# including lag
m_lag <- lmer(
  log(Ne) ~ temperature + lag_temperature + p_human +
    (1 + temperature + lag_temperature + p_human | species),
  data = df
)
# etc...

# compare models
AIC(m_linT, m_full, m_lag)

# body mass as moderator
# temperature * log(adult_mass) expands to:
#   temperature                    - effect of temperature when mass = 0 (i.e. log mass = 0, mass = 1kg)
#   log(adult_mass)                - effect of mass when temperature = 0
#   temperature:log(adult_mass)    - does the temperature slope differ by body mass?
# same logic applies to p_human * log(adult_mass)
m_mass <- lmer(
  log(Ne) ~ temperature * log(adult_mass) +  # climate sensitivity moderated by mass
    p_human * log(adult_mass) +      # human sensitivity moderated by mass
    (1 + temperature + p_human | species),    # each species gets its own intercept
  # and its own slopes for temperature
  # and human pressure
  data = df
)
summary(m_mass)
# generation time as moderator
# same structure as m_mass but using generation time instead
# generation time is arguably more directly linked to demographic recovery rate
# than mass, though the two are correlated
m_gentime <- lmer(
  log(Ne) ~ temperature * log(generation_time) +  # climate sensitivity moderated by generation time
    p_human * log(generation_time) +      # human sensitivity moderated by generation time
    (1 + temperature + p_human | species),
  data = df
)

# both traits as moderators
# this is the most complete version - asks whether mass and generation time
# each explain independent variation in climate and human sensitivity
# after accounting for the other trait
m_both <- lmer(
  log(Ne) ~ temperature * log(adult_mass) +        # climate x mass interaction
    p_human * log(adult_mass) +            # human x mass interaction
    temperature * log(generation_time) +   # climate x generation time interaction
    p_human * log(generation_time) +       # human x generation time interaction
    (1 + temperature + p_human | species),
  data = df
)

# compare models by AIC - lower is better
# if m_gentime beats m_mass, generation time is a better moderator than mass
# if m_both beats both, the two traits carry independent information
AIC(m_mass, m_gentime, m_both)

# summary shows fixed effect estimates and their significance
# key things to look at:
#   temperature:log(adult_mass)     - is climate sensitivity greater in larger species?
#   p_human:log(adult_mass)         - is human impact greater in larger species?
#   temperature:log(generation_time)- is climate sensitivity greater in slower species?
#   p_human:log(generation_time)    - is human impact greater in slower species?
summary(m_both)
