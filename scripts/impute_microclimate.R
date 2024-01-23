library(tidyverse)
library(itertools, lib.loc = "/projappl/project_2003061/Rpackages")
library(missForest, lib.loc = "/projappl/project_2003061/Rpackages")
library(lubridate)
library(zoo)
library(foreach)
library(doMC)
library(parallel)

setwd("/projappl/project_2007415/repos/Model_Kilpisjarvi_microclimates/")

d <- bind_rows(read_csv("../MICROCLIMATES/output/KIL/tomst_data_cleaned.csv") %>% 
                 rename(probl = error_tomst) %>% 
                 select(-tomst_id),
               read_csv("data/tomst_data/Oulu_tomst_data.csv"),
               read_csv("data/tomst_data/Roslin_tomst_data.csv")) %>% 
  mutate(datetime = with_tz(datetime, tzone = "Etc/GMT-2")) %>% 
  filter(minute(datetime) %in% c(0,30,60)) %>% 
  mutate(date = as_date(datetime)) %>% 
  relocate(date, .after = datetime) %>% 
  filter(date >= "2019-09-01") %>% 
  filter(date <= "2023-08-31") %>% 
  filter(!startsWith(site, "L"),
         !startsWith(site, "AE"),
         !startsWith(site, "SE"),
         !startsWith(site, "SAL"),
         !startsWith(site, "X"))

# Calibration function for the moisture count values for unknown soils from Kopecky et al. 2020
cal_funNA <- function(x) {((-1.34e-8) * (x^2) + (2.50e-4) * x + (-1.58e-1))*100 }

# Calibrate the moisture values
d %>% mutate(moist = round(cal_funNA(moist),1)) %>% 
  mutate(moist = ifelse(T1 < 1 | probl == 1, NA, moist)) -> d

d <- d %>% 
  relocate(moist, .after = T3)

# Fill missing observations (max length of three consecutive NAs) with running average
# Linear interpolation
d %>% group_by(site) %>% 
  mutate(across(T1:moist, ~na.approx(.x, datetime, maxgap = 3, na.rm = F))) -> d

d <- d %>% select(-T2)

########################################################################
# AGGREGATE TO DAILY VALUES

notNA_prop <- function(x){ round(sum(is.finite(x)/length(x))*100,1) }
gc()

system.time({
  daily <- mclapply(unique(d$site), function(x){
    print(x)
    suppressWarnings({
      d %>% 
        filter(site == x) %>% 
        group_by(site, date) %>%
        summarise(across(T1:moist, ~notNA_prop(.x), 
                         .names = "{.col}_prop"),
                  across(T1:T3, list(mean = mean, min = min, max = max), 
                         na.rm = T, .names = "{.col}_{.fn}"),
                  moist_mean = mean(moist, na.rm = T),
                  probl = max(probl, na.rm = T),
                  .groups = "drop") -> daily
    })
    
    return(daily)
  }, mc.cores = future::availableCores()) %>% bind_rows()
})
gc()

############################################################
# T1 imputation
registerDoMC(future::availableCores())

# T1 mean
T1mean <- daily %>% 
  select(site, date, probl, T1_prop, T1_mean) %>% 
  filter(T1_prop >= 90) %>% 
  filter(probl %in% c(0,3,5:8,10)) %>% 
  pivot_wider(id_cols = date, names_from = site, values_from = T1_mean) %>% 
  arrange(date)

system.time({
  imp_T1mean <- missForest(T1mean %>% mutate(id = 1:nrow(.)) %>% mutate(doy = yday(date)) %>% select(-c(date)) %>% as.data.frame(),
                     maxiter = 10, ntree = 100, variablewise = F, parallelize = 'variables')
})
imp_T1mean$OOBerror

# T1 max
T1max <- daily %>% 
  select(site, date, probl, T1_prop, T1_max) %>% 
  filter(T1_prop >= 90) %>% 
  filter(probl %in% c(0,3,5:8,10)) %>% 
  pivot_wider(id_cols = date, names_from = site, values_from = T1_max) %>% 
  arrange(date)

system.time({
  imp_T1max <- missForest(T1max %>% mutate(id = 1:nrow(.)) %>% mutate(doy = yday(date)) %>% select(-c(date)) %>% as.data.frame(),
                           maxiter = 10, ntree = 100, variablewise = F, parallelize = 'variables')
})
imp_T1max$OOBerror

# T1 min
T1min <- daily %>% 
  select(site, date, probl, T1_prop, T1_min) %>% 
  filter(T1_prop >= 90) %>% 
  filter(probl %in% c(0,3,5:8,10)) %>% 
  pivot_wider(id_cols = date, names_from = site, values_from = T1_min) %>% 
  arrange(date)

system.time({
  imp_T1min <- missForest(T1min %>% mutate(id = 1:nrow(.)) %>% mutate(doy = yday(date)) %>% select(-c(date)) %>% as.data.frame(),
                          maxiter = 10, ntree = 100, variablewise = F, parallelize = 'variables')
})
imp_T1min$OOBerror

####################################################################
# T3 imputation

# T3 mean
T3mean <- daily %>% 
  select(site, date, probl, T3_prop, T3_mean) %>% 
  filter(T3_prop >= 90) %>% 
  filter(probl %in% c(0,6,9,11)) %>% 
  pivot_wider(id_cols = date, names_from = site, values_from = T3_mean) %>% 
  arrange(date)

system.time({
  imp_T3mean <- missForest(T3mean %>% mutate(id = 1:nrow(.)) %>% mutate(doy = yday(date)) %>% select(-c(date)) %>% as.data.frame(),
                           maxiter = 10, ntree = 100, variablewise = F, parallelize = 'variables')
})
imp_T3mean$OOBerror

# T3 max
T3max <- daily %>% 
  select(site, date, probl, T3_prop, T3_max) %>% 
  filter(T3_prop >= 90) %>% 
  filter(probl %in% c(0,6,9,11)) %>% 
  pivot_wider(id_cols = date, names_from = site, values_from = T3_max) %>% 
  arrange(date)

system.time({
  imp_T3max <- missForest(T3max %>% mutate(id = 1:nrow(.)) %>% mutate(doy = yday(date)) %>% select(-c(date)) %>% as.data.frame(),
                          maxiter = 10, ntree = 100, variablewise = F, parallelize = 'variables')
})
imp_T3max$OOBerror

# T3 min
T3min <- daily %>% 
  select(site, date, probl, T3_prop, T3_min) %>% 
  filter(T3_prop >= 90) %>% 
  filter(probl %in% c(0,6,9,11)) %>% 
  pivot_wider(id_cols = date, names_from = site, values_from = T3_min) %>% 
  arrange(date)

system.time({
  imp_T3min <- missForest(T3min %>% mutate(id = 1:nrow(.)) %>% mutate(doy = yday(date)) %>% select(-c(date)) %>% as.data.frame(),
                          maxiter = 10, ntree = 100, variablewise = F, parallelize = 'variables')
})
imp_T3min$OOBerror

# Moisture
moist <- daily %>% 
  select(site, date, probl, moist_prop, moist_mean) %>% 
  filter(moist_prop >= 90) %>% 
  filter(probl %in% c(0,3:5,7,11)) %>% 
  filter(month(date) %in% 6:8) %>% 
  pivot_wider(id_cols = date, names_from = site, values_from = moist_mean) %>% 
  arrange(date)

system.time({
  imp_moist <- missForest(moist %>% mutate(id = 1:nrow(.)) %>% mutate(doy = yday(date)) %>% select(-c(date)) %>% as.data.frame(),
                          maxiter = 10, ntree = 100, variablewise = F, parallelize = 'variables')
})
imp_moist$OOBerror


dall <- full_join(full_join(T1mean %>% mutate(id = 1:nrow(.)) %>% mutate(doy = yday(date)) %>% select(id, date, doy),
                            imp_T1mean$ximp) %>% select(-id, -doy) %>%
                    pivot_longer(cols = c(-date), names_to = "site", values_to = "T1_mean"),
                  full_join(T1min %>% mutate(id = 1:nrow(.)) %>% mutate(doy = yday(date)) %>% select(id, date, doy),
                            imp_T1min$ximp) %>% select(-id, -doy) %>%
                    pivot_longer(cols = c(-date), names_to = "site", values_to = "T1_min")) %>% 
  full_join(.,
            full_join(T1max %>% mutate(id = 1:nrow(.)) %>% mutate(doy = yday(date)) %>% select(id, date, doy),
                      imp_T1max$ximp) %>% select(-id, -doy) %>%
              pivot_longer(cols = c(-date), names_to = "site", values_to = "T1_max")) %>% 
  full_join(.,
            full_join(T3mean %>% mutate(id = 1:nrow(.)) %>% mutate(doy = yday(date)) %>% select(id, date, doy),
                      imp_T3mean$ximp) %>% select(-id, -doy) %>%
              pivot_longer(cols = c(-date), names_to = "site", values_to = "T3_mean")) %>% 
  full_join(.,
            full_join(T3min %>% mutate(id = 1:nrow(.)) %>% mutate(doy = yday(date)) %>% select(id, date, doy),
                      imp_T3min$ximp) %>% select(-id, -doy) %>%
              pivot_longer(cols = c(-date), names_to = "site", values_to = "T3_min")) %>% 
  full_join(.,
            full_join(T3max %>% mutate(id = 1:nrow(.)) %>% mutate(doy = yday(date)) %>% select(id, date, doy),
                      imp_T3max$ximp) %>% select(-id, -doy) %>%
              pivot_longer(cols = c(-date), names_to = "site", values_to = "T3_max")) %>% 
  full_join(.,
            full_join(moist %>% mutate(id = 1:nrow(.)) %>% mutate(doy = yday(date)) %>% select(id, date, doy),
                      imp_moist$ximp) %>% select(-id, -doy) %>%
              pivot_longer(cols = c(-date), names_to = "site", values_to = "moist")) %>% 
  mutate(moist = ifelse(T1_mean < 1, NA, moist))

write_csv(dall, "data/imputed_tomst_data.csv")
