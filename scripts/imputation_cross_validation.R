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

T1month <- daily %>% 
  select(site, date, probl, T1_prop, T1_mean) %>% 
  filter(T1_prop >= 90) %>% 
  filter(probl %in% c(0,3,5:8,10)) %>% 
  mutate(month = month(date),
         year = year(date)) %>% 
  group_by(site, year, month) %>% 
  summarise(T1_mean = mean(T1_mean, na.rm = FALSE)) %>% 
  drop_na() %>% 
  ungroup()

impsT1mean <- lapply(1:10, function(i){
  print(i)
  T1s <- T1month %>% sample_n(size = 20)
  
  T1mean <- daily %>% 
    select(site, date, probl, T1_prop, T1_mean) %>% 
    filter(T1_prop >= 90) %>% 
    filter(probl %in% c(0,3,5:8,10)) %>% 
    mutate(month = month(date),
           year = year(date)) %>% 
    left_join(.,
              T1s %>% select(-T1_mean) %>% mutate(exc = TRUE))
  
  T1mean_test <- T1mean %>% 
    filter(exc) %>% select(-exc) %>% 
    select(site, date, month, year, T1_mean)
  
  T1mean <- T1mean %>% 
    filter(is.na(exc)) %>% select(-exc) %>% 
    pivot_wider(id_cols = date, names_from = site, values_from = T1_mean) %>% 
    arrange(date)
  
  system.time({
    imp_T1mean <- missForest(T1mean %>% mutate(id = 1:nrow(.)) %>% mutate(doy = yday(date)) %>% select(-c(date)) %>% as.data.frame(),
                             maxiter = 10, ntree = 100, variablewise = F, parallelize = 'variables')
  })
  
  imps <- imp_T1mean$ximp %>% 
    select(-id, -doy) %>% 
    mutate(date = T1mean$date) %>% 
    pivot_longer(cols = -date, names_to = "site", values_to = "T1_mean_imp") %>% 
    left_join(., T1mean_test) %>% 
    drop_na()
  
  return(imps)
}) %>% bind_rows()


# T1 max

T1month <- daily %>% 
  select(site, date, probl, T1_prop, T1_max) %>% 
  filter(T1_prop >= 90) %>% 
  filter(probl %in% c(0,3,5:8,10)) %>% 
  mutate(month = month(date),
         year = year(date)) %>% 
  group_by(site, year, month) %>% 
  summarise(T1_max = mean(T1_max, na.rm = FALSE)) %>% 
  drop_na() %>% 
  ungroup()

impsT1max <- lapply(1:10, function(i){
  print(i)
  T1s <- T1month %>% sample_n(size = 20)
  
  T1max <- daily %>% 
    select(site, date, probl, T1_prop, T1_max) %>% 
    filter(T1_prop >= 90) %>% 
    filter(probl %in% c(0,3,5:8,10)) %>% 
    mutate(month = month(date),
           year = year(date)) %>% 
    left_join(.,
              T1s %>% select(-T1_max) %>% mutate(exc = TRUE))
  
  T1max_test <- T1max %>% 
    filter(exc) %>% select(-exc) %>% 
    select(site, date, month, year, T1_max)
  
  T1max <- T1max %>% 
    filter(is.na(exc)) %>% select(-exc) %>% 
    pivot_wider(id_cols = date, names_from = site, values_from = T1_max) %>% 
    arrange(date)
  
  
  system.time({
    imp_T1max <- missForest(T1max %>% mutate(id = 1:nrow(.)) %>% mutate(doy = yday(date)) %>% select(-c(date)) %>% as.data.frame(),
                            maxiter = 10, ntree = 100, variablewise = F, parallelize = 'variables')
  })
  
  imps <- imp_T1max$ximp %>% 
    select(-id, -doy) %>% 
    mutate(date = T1max$date) %>% 
    pivot_longer(cols = -date, names_to = "site", values_to = "T1_max_imp") %>% 
    left_join(., T1max_test) %>% 
    drop_na()
  
  return(imps)
}) %>% bind_rows()


# T1 min

T1month <- daily %>% 
  select(site, date, probl, T1_prop, T1_min) %>% 
  filter(T1_prop >= 90) %>% 
  filter(probl %in% c(0,3,5:8,10)) %>% 
  mutate(month = month(date),
         year = year(date)) %>% 
  group_by(site, year, month) %>% 
  summarise(T1_min = mean(T1_min, na.rm = FALSE)) %>% 
  drop_na() %>% 
  ungroup()

impsT1min <- lapply(1:10, function(i){
  print(i)
  T1s <- T1month %>% sample_n(size = 20)
  
  T1min <- daily %>% 
    select(site, date, probl, T1_prop, T1_min) %>% 
    filter(T1_prop >= 90) %>% 
    filter(probl %in% c(0,3,5:8,10)) %>% 
    mutate(month = month(date),
           year = year(date)) %>% 
    left_join(.,
              T1s %>% select(-T1_min) %>% mutate(exc = TRUE))
  
  T1min_test <- T1min %>% 
    filter(exc) %>% select(-exc) %>% 
    select(site, date, month, year, T1_min)
  
  T1min <- T1min %>% 
    filter(is.na(exc)) %>% select(-exc) %>% 
    pivot_wider(id_cols = date, names_from = site, values_from = T1_min) %>% 
    arrange(date)
  
  
  system.time({
    imp_T1min <- missForest(T1min %>% mutate(id = 1:nrow(.)) %>% mutate(doy = yday(date)) %>% select(-c(date)) %>% as.data.frame(),
                            maxiter = 10, ntree = 100, variablewise = F, parallelize = 'variables')
  })
  
  imps <- imp_T1min$ximp %>% 
    select(-id, -doy) %>% 
    mutate(date = T1min$date) %>% 
    pivot_longer(cols = -date, names_to = "site", values_to = "T1_min_imp") %>% 
    left_join(., T1min_test) %>% 
    drop_na()
  
  return(imps)
}) %>% bind_rows()

####################################################################
# T3 imputation

# T3 mean

T3month <- daily %>% 
  select(site, date, probl, T3_prop, T3_mean) %>% 
  filter(T3_prop >= 90) %>% 
  filter(probl %in% c(0,6,9,11)) %>% 
  mutate(month = month(date),
         year = year(date)) %>% 
  group_by(site, year, month) %>% 
  summarise(T3_mean = mean(T3_mean, na.rm = FALSE)) %>% 
  drop_na() %>% 
  ungroup()

impsT3mean <- lapply(1:10, function(i){
  print(i)
  T3s <- T3month %>% sample_n(size = 20)
  
  T3mean <- daily %>% 
    select(site, date, probl, T3_prop, T3_mean) %>% 
    filter(T3_prop >= 90) %>% 
    filter(probl %in% c(0,6,9,11)) %>% 
    mutate(month = month(date),
           year = year(date)) %>% 
    left_join(.,
              T3s %>% select(-T3_mean) %>% mutate(exc = TRUE))
  
  T3mean_test <- T3mean %>% 
    filter(exc) %>% select(-exc) %>% 
    select(site, date, month, year, T3_mean)
  
  T3mean <- T3mean %>% 
    filter(is.na(exc)) %>% select(-exc) %>% 
    pivot_wider(id_cols = date, names_from = site, values_from = T3_mean) %>% 
    arrange(date)
  
  system.time({
    imp_T3mean <- missForest(T3mean %>% mutate(id = 1:nrow(.)) %>% mutate(doy = yday(date)) %>% select(-c(date)) %>% as.data.frame(),
                             maxiter = 10, ntree = 100, variablewise = F, parallelize = 'variables')
  })
  
  imps <- imp_T3mean$ximp %>% 
    select(-id, -doy) %>% 
    mutate(date = T3mean$date) %>% 
    pivot_longer(cols = -date, names_to = "site", values_to = "T3_mean_imp") %>% 
    left_join(., T3mean_test) %>% 
    drop_na()
  
  return(imps)
}) %>% bind_rows()


# T3 max

T3month <- daily %>% 
  select(site, date, probl, T3_prop, T3_max) %>% 
  filter(T3_prop >= 90) %>% 
  filter(probl %in% c(0,6,9,11)) %>% 
  mutate(month = month(date),
         year = year(date)) %>% 
  group_by(site, year, month) %>% 
  summarise(T3_max = mean(T3_max, na.rm = FALSE)) %>% 
  drop_na() %>% 
  ungroup()

impsT3max <- lapply(1:10, function(i){
  print(i)
  T3s <- T3month %>% sample_n(size = 20)
  
  T3max <- daily %>% 
    select(site, date, probl, T3_prop, T3_max) %>% 
    filter(T3_prop >= 90) %>% 
    filter(probl %in% c(0,6,9,11)) %>% 
    mutate(month = month(date),
           year = year(date)) %>% 
    left_join(.,
              T3s %>% select(-T3_max) %>% mutate(exc = TRUE))
  
  T3max_test <- T3max %>% 
    filter(exc) %>% select(-exc) %>% 
    select(site, date, month, year, T3_max)
  
  T3max <- T3max %>% 
    filter(is.na(exc)) %>% select(-exc) %>% 
    pivot_wider(id_cols = date, names_from = site, values_from = T3_max) %>% 
    arrange(date)
  
  
  system.time({
    imp_T3max <- missForest(T3max %>% mutate(id = 1:nrow(.)) %>% mutate(doy = yday(date)) %>% select(-c(date)) %>% as.data.frame(),
                            maxiter = 10, ntree = 100, variablewise = F, parallelize = 'variables')
  })
  
  imps <- imp_T3max$ximp %>% 
    select(-id, -doy) %>% 
    mutate(date = T3max$date) %>% 
    pivot_longer(cols = -date, names_to = "site", values_to = "T3_max_imp") %>% 
    left_join(., T3max_test) %>% 
    drop_na()
  
  return(imps)
}) %>% bind_rows()


# T3 min

T3month <- daily %>% 
  select(site, date, probl, T3_prop, T3_min) %>% 
  filter(T3_prop >= 90) %>% 
  filter(probl %in% c(0,6,9,11)) %>% 
  mutate(month = month(date),
         year = year(date)) %>% 
  group_by(site, year, month) %>% 
  summarise(T3_min = mean(T3_min, na.rm = FALSE)) %>% 
  drop_na() %>% 
  ungroup()

impsT3min <- lapply(1:10, function(i){
  print(i)
  T3s <- T3month %>% sample_n(size = 20)
  
  T3min <- daily %>% 
    select(site, date, probl, T3_prop, T3_min) %>% 
    filter(T3_prop >= 90) %>% 
    filter(probl %in% c(0,6,9,11)) %>% 
    mutate(month = month(date),
           year = year(date)) %>% 
    left_join(.,
              T3s %>% select(-T3_min) %>% mutate(exc = TRUE))
  
  T3min_test <- T3min %>% 
    filter(exc) %>% select(-exc) %>% 
    select(site, date, month, year, T3_min)
  
  T3min <- T3min %>% 
    filter(is.na(exc)) %>% select(-exc) %>% 
    pivot_wider(id_cols = date, names_from = site, values_from = T3_min) %>% 
    arrange(date)
  
  
  system.time({
    imp_T3min <- missForest(T3min %>% mutate(id = 1:nrow(.)) %>% mutate(doy = yday(date)) %>% select(-c(date)) %>% as.data.frame(),
                            maxiter = 10, ntree = 100, variablewise = F, parallelize = 'variables')
  })
  
  imps <- imp_T3min$ximp %>% 
    select(-id, -doy) %>% 
    mutate(date = T3min$date) %>% 
    pivot_longer(cols = -date, names_to = "site", values_to = "T3_min_imp") %>% 
    left_join(., T3min_test) %>% 
    drop_na()
  
  return(imps)
}) %>% bind_rows()

####################################################################
# Moisture

moist <- daily %>% 
  select(site, date, probl, moist_prop, moist_mean) %>% 
  filter(moist_prop >= 90) %>% 
  filter(probl %in% c(0,3:5,7,11)) %>% 
  filter(month(date) %in% 6:8) %>% 
  pivot_wider(id_cols = date, names_from = site, values_from = moist_mean) %>% 
  arrange(date)

moistmonth <- daily %>% 
  select(site, date, probl, moist_prop, moist_mean) %>% 
  filter(moist_prop >= 90) %>% 
  filter(probl %in% c(0,3:5,7,11)) %>% 
  mutate(month = month(date),
         year = year(date)) %>% 
  group_by(site, year, month) %>% 
  summarise(moist_mean = mean(moist_mean, na.rm = FALSE)) %>% 
  drop_na() %>% 
  ungroup()

impsmoistmean <- lapply(1:10, function(i){
  print(i)
  moists <- moistmonth %>% sample_n(size = 20)
  
  moistmean <- daily %>% 
    select(site, date, probl, moist_prop, moist_mean) %>% 
    filter(moist_prop >= 90) %>% 
    filter(probl %in% c(0,6,9,11)) %>% 
    mutate(month = month(date),
           year = year(date)) %>% 
    left_join(.,
              moists %>% select(-moist_mean) %>% mutate(exc = TRUE))
  
  moistmean_test <- moistmean %>% 
    filter(exc) %>% select(-exc) %>% 
    select(site, date, month, year, moist_mean)
  
  moistmean <- moistmean %>% 
    filter(is.na(exc)) %>% select(-exc) %>% 
    pivot_wider(id_cols = date, names_from = site, values_from = moist_mean) %>% 
    arrange(date)
  
  system.time({
    imp_moistmean <- missForest(moistmean %>% mutate(id = 1:nrow(.)) %>% mutate(doy = yday(date)) %>% select(-c(date)) %>% as.data.frame(),
                                maxiter = 10, ntree = 100, variablewise = F, parallelize = 'variables')
  })
  
  imps <- imp_moistmean$ximp %>% 
    select(-id, -doy) %>% 
    mutate(date = moistmean$date) %>% 
    pivot_longer(cols = -date, names_to = "site", values_to = "moist_mean_imp") %>% 
    left_join(., moistmean_test) %>% 
    drop_na()
  
  return(imps)
}) %>% bind_rows()

###########################################################
# COMBINE

imps <- bind_rows(impsT1mean %>% rename(imp = T1_mean_imp, ori = T1_mean) %>% mutate(var = "T1_mean"),
          impsT1max %>% rename(imp = T1_max_imp, ori = T1_max) %>% mutate(var = "T1_max"),
          impsT1min %>% rename(imp = T1_min_imp, ori = T1_min) %>% mutate(var = "T1_min"),
          impsT3mean %>% rename(imp = T3_mean_imp, ori = T3_mean) %>% mutate(var = "T3_mean"),
          impsT3max %>% rename(imp = T3_max_imp, ori = T3_max) %>% mutate(var = "T3_max"),
          impsT3min %>% rename(imp = T3_min_imp, ori = T3_min) %>% mutate(var = "T3_min"),
          impsmoistmean %>% rename(imp = moist_mean_imp, ori = moist_mean) %>% mutate(var = "moist")) %>% 
  mutate(dif = ori - imp)

imps %>% 
  group_by(var) %>% 
  summarise(abs_err = mean(abs(dif)),
            median_abs_err = median(abs(dif)),
            bias = mean(dif),
            rmse = sqrt(mean(dif^2))) %>% 
  mutate(across(-var, ~round(.x, 3))) %>% view

imps %>% 
  group_by(site, month, year, var) %>% 
  summarise(imp = mean(imp),
            ori = mean(ori)) %>% 
  mutate(dif = ori - imp) %>% 
  group_by(var) %>% 
  summarise(mean_abs_err = mean(abs(dif)),
            median_abs_err = median(abs(dif)),
            bias = mean(dif),
            rmse = sqrt(mean(dif^2)))

imps %>% 
  filter(var != "moist") %>% 
  group_by(site, month, year, var) %>% 
  summarise(imp = mean(imp),
            ori = mean(ori)) %>% 
  mutate(dif = ori - imp) %>% 
  group_by(month) %>% 
  summarise(mean_abs_err = mean(abs(dif)),
            median_abs_err = median(abs(dif)),
            bias = mean(dif),
            rmse = sqrt(mean(dif^2)))

write_csv(imps, "output/imputaton_cv_stats.csv")
