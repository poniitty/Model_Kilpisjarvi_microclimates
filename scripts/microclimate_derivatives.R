library(tidyverse) # Handles everything more nicely than base R
library(lubridate) # Handles time objects nicely

d <- read_csv("data/imputed_tomst_data.csv")

# Own helper function to round to two decimals
round2 <- function(x) round(x,2)
# Own function to change infinite values to NA
infmutate <- function(x) ifelse(is.infinite(x),NA,x)

d <- d %>% mutate(hyddate = date + days(122)) %>% 
  mutate(year = year(hyddate),
         month = month(date),
         doy = yday(date),
         hyddoy = yday(hyddate)) %>% 
  relocate(hyddate:hyddoy, .after = date)

########################################################################
# Annual means

am <- d %>% 
  group_by(site) %>% 
  summarise(T1_mean = mean(T1_mean),
            T3_mean = mean(T3_mean))

########################################################################
# MONTHLY means

mm <- d %>% 
  group_by(site, month) %>% 
  summarise(T1_mean = mean(T1_mean),
            T3_mean = mean(T3_mean)) %>% 
  pivot_wider(id_cols = site, names_from = month, values_from = T1_mean:T3_mean)

########################################################################
# Degree days

tdds <- d %>% 
  mutate(T1_mean = ifelse(T1_mean < 0, 0, T1_mean),
         T3_mean = ifelse(T3_mean < 0, 0, T3_mean)) %>% 
  group_by(site, year) %>% 
  summarise(T1_TDD = sum(T1_mean),
            T3_TDD = sum(T3_mean)) %>% 
  group_by(site) %>% 
  summarise(T1_TDD = mean(T1_TDD),
            T3_TDD = mean(T3_TDD))

gdd3s <- d %>% 
  mutate(T1_mean = ifelse(T1_mean < 3, 0, T1_mean),
         T3_mean = ifelse(T3_mean < 3, 0, T3_mean)) %>% 
  group_by(site, year) %>% 
  summarise(T1_GDD3 = sum(T1_mean),
            T3_GDD3 = sum(T3_mean)) %>% 
  group_by(site) %>% 
  summarise(T1_GDD3 = mean(T1_GDD3),
            T3_GDD3 = mean(T3_GDD3))

gdd5s <- d %>% 
  mutate(T1_mean = ifelse(T1_mean < 5, 0, T1_mean),
         T3_mean = ifelse(T3_mean < 5, 0, T3_mean)) %>% 
  group_by(site, year) %>% 
  summarise(T1_GDD5 = sum(T1_mean),
            T3_GDD5 = sum(T3_mean)) %>% 
  group_by(site) %>% 
  summarise(T1_GDD5 = mean(T1_GDD5),
            T3_GDD5 = mean(T3_GDD5))

fdds <- d %>% 
  mutate(T1_mean = ifelse(T1_mean > 0, 0, T1_mean),
         T3_mean = ifelse(T3_mean > 0, 0, T3_mean)) %>% 
  group_by(site, year) %>% 
  summarise(T1_FDD = sum(T1_mean),
            T3_FDD = sum(T3_mean)) %>% 
  group_by(site) %>% 
  summarise(T1_FDD = mean(T1_FDD),
            T3_FDD = mean(T3_FDD))


# thermal sum days

tdds_days <- d %>% 
  mutate(T1_mean = ifelse(T1_mean < 0, 0, 1),
         T3_mean = ifelse(T3_mean < 0, 0, 1)) %>% 
  group_by(site, year) %>% 
  summarise(T1_TDD_days = sum(T1_mean),
            T3_TDD_days = sum(T3_mean)) %>% 
  group_by(site) %>% 
  summarise(T1_TDD_days = mean(T1_TDD_days),
            T3_TDD_days = mean(T3_TDD_days))

gdd3s_days <- d %>% 
  mutate(T1_mean = ifelse(T1_mean < 3, 0, 1),
         T3_mean = ifelse(T3_mean < 3, 0, 1)) %>% 
  group_by(site, year) %>% 
  summarise(T1_GDD3_days = sum(T1_mean),
            T3_GDD3_days = sum(T3_mean)) %>% 
  group_by(site) %>% 
  summarise(T1_GDD3_days = mean(T1_GDD3_days),
            T3_GDD3_days = mean(T3_GDD3_days))

gdd5s_days <- d %>% 
  mutate(T1_mean = ifelse(T1_mean < 5, 0, 1),
         T3_mean = ifelse(T3_mean < 5, 0, 1)) %>% 
  group_by(site, year) %>% 
  summarise(T1_GDD5_days = sum(T1_mean),
            T3_GDD5_days = sum(T3_mean)) %>% 
  group_by(site) %>% 
  summarise(T1_GDD5_days = mean(T1_GDD5_days),
            T3_GDD5_days = mean(T3_GDD5_days))

fdds_days <- d %>% 
  mutate(T1_mean = ifelse(T1_mean > 0, 0, 1),
         T3_mean = ifelse(T3_mean > 0, 0, 1)) %>% 
  group_by(site, year) %>% 
  summarise(T1_FDD_days = sum(T1_mean),
            T3_FDD_days = sum(T3_mean)) %>% 
  group_by(site) %>% 
  summarise(T1_FDD_days = mean(T1_FDD_days),
            T3_FDD_days = mean(T3_FDD_days))

# Min & max

mins <- d %>% 
  group_by(site, year) %>% 
  summarise(T1_min = quantile(T1_min, 0.005, na.rm = T),
            T3_min = quantile(T3_min, 0.005, na.rm = T)) %>% 
  group_by(site) %>% 
  summarise(T1_min = mean(T1_min),
            T3_min = mean(T3_min))

maxs <- d %>% 
  group_by(site, year) %>% 
  summarise(T1_max = quantile(T1_max, 0.995, na.rm = T),
            T3_max = quantile(T3_max, 0.995, na.rm = T)) %>% 
  group_by(site) %>% 
  summarise(T1_max = mean(T1_max),
            T3_max = mean(T3_max))

# Summer frost

frost <- d %>% 
  filter(month %in% c(5:8),
         doy < 227,
         T1_mean >= 1) %>% 
  mutate(T3_min = ifelse(T3_min > 0, 0, T3_min)) %>% 
  group_by(site, year) %>% 
  summarise(T3_frost = sum(T3_min)) %>% 
  group_by(site) %>% 
  summarise(T3_frost = mean(T3_frost)) %>%
  arrange(T3_frost)

# Start of season

sos1 <- d %>% 
  group_by(site, year) %>% 
  filter(T1_mean >= 1 & T3_mean >= 1) %>% 
  summarise(sos1 = min(doy),
            eos1 = max(doy)) %>% 
  group_by(site) %>% 
  summarise(sos1 = mean(sos1),
            eos1 = mean(eos1))
sos3 <- d %>% 
  group_by(site, year) %>% 
  filter(T1_mean >= 3 & T3_mean >= 3) %>% 
  summarise(sos3 = min(doy),
            eos3 = max(doy)) %>% 
  group_by(site) %>% 
  summarise(sos3 = mean(sos3),
            eos3 = mean(eos3))

# Freeze-thaw cycles

T1_ft <- d %>%
  select(date, year, site, T1_min, T1_max) %>% 
  arrange(site, date) %>% 
  pivot_longer(cols = c(T1_min,T1_max), values_to = "T1") %>% 
  arrange(site, date, name) %>% 
  mutate(diff = T1-lag(T1)) %>% 
  mutate(across(T1, ~ifelse(.x >= 0, 1, -1))) %>% 
  group_by(site, year) %>% 
  mutate(T1_ft = ifelse(T1*lag(T1) == -1 & abs(diff) > 0.5, 1, 0)) %>% 
  summarise(T1_ft = sum(T1_ft, na.rm = T)) %>% 
  group_by(site) %>% 
  summarise(T1_ft = mean(T1_ft))

T3_ft <- d %>%
  select(date, year, site, T3_min, T3_max) %>% 
  arrange(site, date) %>% 
  pivot_longer(cols = c(T3_min,T3_max), values_to = "T3") %>% 
  arrange(site, date, name) %>% 
  mutate(diff = T3-lag(T3)) %>% 
  mutate(across(T3, ~ifelse(.x >= 0, 1, -1))) %>% 
  group_by(site, year) %>% 
  mutate(T3_ft = ifelse(T3*lag(T3) == -1 & abs(diff) >= 1, 1, 0)) %>% 
  summarise(T3_ft = sum(T3_ft, na.rm = T)) %>% 
  group_by(site) %>% 
  summarise(T3_ft = mean(T3_ft))

#####################################################################
# Moisture

moist_summer <- d %>% 
  filter(T1_mean >= 1) %>% 
  filter(month %in% c(6:8)) %>% 
  group_by(site, year) %>% 
  summarise(moist_mean_summer = mean(moist),
            moist_max_summer = max(moist),
            moist_min_summer = min(moist),
            moist_sd_summer = sd(moist)) %>% 
  group_by(site) %>% 
  summarise(moist_mean_summer = mean(moist_mean_summer),
            moist_max_summer = mean(moist_max_summer),
            moist_min_summer = mean(moist_min_summer),
            moist_sd_summer = mean(moist_sd_summer)) %>% 
  mutate(moist_cv_summer = moist_sd_summer/moist_mean_summer)

moist_july <- d %>% 
  filter(T1_mean >= 1) %>% 
  filter(month == 7) %>% 
  group_by(site, year) %>% 
  summarise(moist_mean_july = mean(moist),
            moist_max_july = max(moist),
            moist_min_july = min(moist),
            moist_sd_july = sd(moist)) %>% 
  group_by(site) %>% 
  summarise(moist_mean_july = mean(moist_mean_july),
            moist_max_july = mean(moist_max_july),
            moist_min_july = mean(moist_min_july),
            moist_sd_july = mean(moist_sd_july)) %>% 
  mutate(moist_cv_july = moist_sd_july/moist_mean_july)


# Combine datasets

all <- full_join(am, mm) %>% 
  full_join(.,mins) %>% 
  full_join(.,maxs) %>% 
  full_join(.,tdds) %>% 
  full_join(.,gdd3s) %>% 
  full_join(.,gdd5s) %>% 
  full_join(.,fdds) %>% 
  full_join(.,tdds_days) %>% 
  full_join(.,gdd3s_days) %>% 
  full_join(.,gdd5s_days) %>% 
  full_join(.,fdds_days) %>% 
  full_join(.,sos1) %>% 
  full_join(.,sos3) %>% 
  full_join(.,T1_ft) %>% 
  full_join(.,T3_ft) %>% 
  full_join(.,frost) %>% 
  full_join(.,moist_summer) %>% 
  full_join(.,moist_july)

summary(all)

write_csv(all, "data/microclimate_variables.csv")
