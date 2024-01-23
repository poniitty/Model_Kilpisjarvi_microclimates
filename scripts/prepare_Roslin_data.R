
library(tidyverse)
library(lubridate)
library(data.table)
library(cowplot)
library(zoo)

# Set date limits to remove implausible dates
mind <- as.Date("2020-06-01", tz = "Etc/GMT-2")
maxd <- as.Date("2022-09-28", tz = "Etc/GMT-2")

raw_data_dir <- "/scratch/project_2007415/microclim/Roslin/Kilpisjarvi_raw"
# raw_data_dir <- "C:/Users/OMISTAJA/OneDrive - University of Helsinki/Kesä2022/pekka2022"


# List logger data files to read
f <- list.files(raw_data_dir, pattern = "data_", full.names = T, recursive = T)

fi <- data.frame(file = f)

fi$file2 <- gsub("_..csv", "", fi$file)

fi$site <- unlist(lapply(lapply(fi$file, function(x) tail(strsplit(x, "/")[[1]],1)), function(x) strsplit(x, "_")[[1]][1]))

fi <- fi[order(fi$site),]

fi$tomst_id <- unlist(lapply(fi$file, function(x) as.numeric(strsplit(gsub("data_","",rev(strsplit(x, "/")[[1]])[1]), "_")[[1]][2])))

readdata <- function(i){
  nn <- sum(grepl(i, fi$file2))
  
  if(nn > 1){
    
    fi2 <- fi %>% filter(grepl(i, fi$file2))
    
    df2 <- data.frame()
    for(ii in fi2$file2){
      print(ii)
      d <- fread(ii)
      
      d %>% select(V2,V3,V4,V5,V6,V7) -> d
      
      d %>% filter(!duplicated(.$V2, fromLast = T)) -> d
      
      df2 <- bind_rows(df2, d)
    }
    
    df2 %>% filter(!duplicated(.$V2, fromLast = T)) -> df2
    
    df2$site <- fi[which(fi$file2 == ii),"id"]
    d$tomst_id <- fi[which(fi$file2 == i),"tomst_id"]
    
    df2 %>% mutate(across(V4:V6, ~as.numeric(gsub(",",".\\",.)))) -> df2
    
    df2 %>% mutate(V2 = ymd_hm(V2, tz = "UTC")) %>% 
      mutate(V2 = with_tz(V2, tzone = "Etc/GMT-2")) -> df2
    
    
    return(df2)
    
  } else {
    
    print(i)
    d <- fread(fi$file[fi$file2 == i])
    
    d %>% select(V2,V3,V4,V5,V6,V7) -> d
    
    d %>% filter(!duplicated(.$V2, fromLast = T)) -> d
    
    d %>% mutate(across(V4:V6, ~as.numeric(gsub(",",".\\",.)))) -> d
    
    d$site <- fi[which(fi$file2 == i),"site"]
    d$tomst_id <- fi[which(fi$file2 == i),"tomst_id"]
    
    d %>% mutate(V2 = ymd_hm(V2, tz = "UTC")) %>% 
      mutate(V2 = with_tz(V2, tzone = "Etc/GMT-2")) -> d
    
    return(d)
    
  }
  
}

mylist <- lapply(fi$file2, readdata)
df <- rbindlist( mylist )

# Rename columns
df %>% rename(datetime = V2,
              zone = V3,
              T1 = V4,
              T2 = V5,
              T3 = V6,
              moist = V7) -> df

df %>% arrange(site, datetime) -> df

# Exceptions

df %>% group_by(site, tomst_id) %>% 
  summarise(maxdt = max(datetime)) -> maxdt
# maxdt <- full_join(maxdt, fi %>% select(site, tomst_id) %>% filter(!duplicated(.)))

maxdt %>% arrange(maxdt)
maxdt %>% arrange(desc(maxdt))

# Remove implausible dates
df %>% filter(datetime > mind,
              datetime < maxd) -> df

sites <- unique(df$site)

# Remove the measurements of the time of the 2022 visiting
# as reading the logger may influence the measurements
maxdt <- maxdt %>% group_by(site) %>% arrange(desc(maxdt)) %>% slice_head(n = 1) %>% ungroup()
full_join(df, maxdt %>% select(-tomst_id) %>% mutate(visit = 1) %>% rename(datetime = maxdt)) %>% 
  group_by(site, tomst_id) %>%
  mutate(visit = rollapply(visit, width=3, FUN=sum, na.rm = T, 
                           fill = NA, partial = T, align = "center")) %>% 
  mutate(across(T1:moist, ~ifelse(visit == 0, ., NA))) %>%
  select(-visit) -> df


# Calculate different daily values for diagnostics
df %>% mutate(date = as_date(datetime)) %>% 
  group_by(site,date,tomst_id) %>% 
  mutate(roll_diff_T1 = T1 - lead(T1),
         roll_diff_T2 = T2 - lead(T2),
         roll_diff_T3 = T3 - lead(T3)) %>%
  summarise(soil_sd = sd(T1),
            air_sd = sd(T3),
            soil_mean = mean(T1),
            air_mean = mean(T3),
            soil_min = min(T1),
            soil_max = max(T1),
            surf_max = max(T2),
            air_max = max(T3),
            mean_diff = mean(abs(T3-T1)),
            moist = min(moist),
            corr = cor(T1,T3, use = "pairwise.complete.obs"),
            max_diff_T1 = max(roll_diff_T1, na.rm = T),
            max_diff_T2 = max(roll_diff_T2, na.rm = T),
            max_diff_T3 = max(roll_diff_T3, na.rm = T)) %>% 
  mutate(sa_diff = air_mean-soil_mean,
         sa_max_diff = air_max-soil_max,
         ss_max_diff = surf_max-soil_max) %>% 
  mutate(sd_ratio = soil_sd/air_sd) %>% 
  mutate(sd_ratio = ifelse(soil_sd < 1, 0, sd_ratio)) %>% 
  mutate(corr = ifelse(is.na(corr), 0, corr)) %>% 
  as.data.frame() -> df2

# create column for error codes
df2 %>% mutate(probl = 0) -> df2

############################################################################
# PLOTTINGS
############################################################################

# Plot temperatures
pdf("visuals/Temperature_graphs_Roslin.pdf", 12, 5)
for(i in sites){
  # i <- "SAA895"
  print(i)
  
  temp <- df %>% filter(site == i)
  
  if(length(na.omit(unique(temp$tomst_id))) > 1){
    
    for(ii in na.omit(unique(temp$tomst_id))){
      temp %>% 
        filter(tomst_id == ii) %>% 
        #group_by(date) %>% 
        #summarise_at(vars(i, "soil"), funs(mean, min, max), na.rm = T) %>% 
        #lapply(function(x) replace(x, is.infinite(x),NA)) %>% as_tibble() %>% 
        ggplot(aes_string(x="datetime")) +
        geom_line(aes_string(y = "T3"), col = "cornflowerblue") +
        geom_line(aes_string(y = "T2"), col = "brown1") +
        geom_line(aes_string(y = "T1"), col = "darkgoldenrod") +
        theme_minimal() +
        ylab("Temperature") + xlab("Date")+
        scale_y_continuous(limits = c(-20, 35))+
        ggtitle(paste0(i,"_",ii)) -> GG
      print(GG)
    }
    
  } else {
    temp %>% 
      #group_by(date) %>% 
      #summarise_at(vars(i, "soil"), funs(mean, min, max), na.rm = T) %>% 
      #lapply(function(x) replace(x, is.infinite(x),NA)) %>% as_tibble() %>% 
      ggplot(aes_string(x="datetime")) +
      geom_line(aes_string(y = "T3"), col = "cornflowerblue") +
      geom_line(aes_string(y = "T2"), col = "brown1") +
      geom_line(aes_string(y = "T1"), col = "darkgoldenrod") +
      theme_minimal() +
      ylab("Temperature") + xlab("Date")+
      scale_y_continuous(limits = c(-20, 35))+
      ggtitle(i) -> GG
    print(GG)
  }
}
dev.off()

# pdf("visuals/Temperature_diagnose_graphs.pdf", 18, 11)
# for(i in sites){
#   print(i)
#   df %>% filter(site == i) %>%
#     ggplot(aes_string(x="datetime")) +
#     geom_line(aes_string(y = "T3"), col = "cornflowerblue") +
#     geom_line(aes_string(y = "T2"), col = "brown1") +
#     geom_line(aes_string(y = "T1"), col = "darkgoldenrod") +
#     theme_minimal() +
#     ylab("Temperature") + xlab("Date")+
#     scale_y_continuous(limits = c(-20, 35))+
#     ggtitle(i) -> GG1
#   
#   df2 %>% filter(site == i) %>%
#     ggplot(aes_string(x="date")) +
#     geom_line(aes_string(y = "air_sd"), col = "cornflowerblue") +
#     geom_line(aes_string(y = "soil_sd"), col = "darkgoldenrod") +
#     theme_minimal() +
#     ylab("Temperature SD") + xlab("Date") -> GG2
#   
#   df2 %>% filter(site == i) %>%
#     ggplot(aes_string(x="date")) +
#     geom_line(aes_string(y = "sd_ratio"), col = "black") +
#     theme_minimal() +
#     ylab("SD ratio") + xlab("Date") -> GG3
#   
#   df2 %>% filter(site == i) %>%
#     ggplot(aes_string(x="date")) +
#     geom_line(aes_string(y = "air_mean"), col = "cornflowerblue") +
#     geom_line(aes_string(y = "soil_mean"), col = "darkgoldenrod") +
#     geom_line(aes_string(y = "air_max"), col = "cornflowerblue") +
#     geom_line(aes_string(y = "surf_max"), col = "brown1") +
#     geom_line(aes_string(y = "soil_max"), col = "darkgoldenrod") +
#     theme_minimal() +
#     ylab("Temperature mean + max") + xlab("Date")+
#     scale_y_continuous(limits = c(-25, 40))-> GG4
#   
#   df2 %>% filter(site == i) %>%
#     ggplot(aes_string(x="date")) +
#     geom_line(aes_string(y = "corr")) +
#     theme_minimal() +
#     ylab("Air - soil correlation") + xlab("Date")+
#     scale_y_continuous(limits = c(-1, 1)) -> GG5
#   
#   df2 %>% filter(site == i) %>%
#     ggplot(aes_string(x="date")) +
#     geom_line(aes_string(y = "max_diff_T3"), col = "cornflowerblue") +
#     geom_line(aes_string(y = "max_diff_T2"), col = "brown1") +
#     geom_line(aes_string(y = "max_diff_T1"), col = "darkgoldenrod") +
#     theme_minimal() +
#     ylab("rolling diffs") + xlab("Date") -> GG6
#   
#   df2 %>% filter(site == i) %>%
#     ggplot(aes_string(x="date")) +
#     geom_line(aes_string(y = "sa_diff"), col = "cornflowerblue") +
#     geom_line(aes_string(y = "sa_max_diff"), col = "brown1") +
#     geom_line(aes_string(y = "ss_max_diff"), col = "darkgoldenrod") +
#     theme_minimal() +
#     ylab("Differences between") + xlab("Date") -> GG7
#   
#   df2 %>% filter(site == i) %>%
#     ggplot(aes_string(x="date")) +
#     geom_line(aes_string(y = "moist"), col = "black") +
#     theme_minimal() +
#     ylab("moisture") + xlab("Date") -> GG8
#   
#   print(plot_grid(plotlist = list(GG1,GG2,GG3,GG4,GG5,GG6,GG7,GG8), nrow = 4))
# }
# dev.off()

###################################################################################################################

# Months to plot
times <- seq(floor_date(as_date(min(df2$date)), "month"),
             ceiling_date(as_date(max(df2$date)), "month") + months(1) - days(1),
             by = "month")

# Plot each site month by month
for(siteid in sites){
  # siteid <- "SAA1195"
  print(siteid)
  pdf(paste0("visuals/monthly_", gsub("/",".",siteid), ".pdf"), 10, 6)
  temp <- df %>% filter(site == siteid)
  
  if(length(na.omit(unique(temp$tomst_id))) > 1){
    
    for(ii in na.omit(unique(temp$tomst_id))){
      
      for (tt in 1:(length(times) - 1)) {
        temp %>% filter(tomst_id == ii) %>%
          filter(datetime >= ymd(times[tt]),
                 datetime < ymd(times[tt + 1])) -> dft
        
        if(nrow(dft %>% filter(complete.cases(.)) > 0)){
          dft %>%
            ggplot(aes_string(x = "datetime")) +
            geom_line(aes_string(y = "T3"), col = "cornflowerblue") +
            geom_line(aes_string(y = "T2"), col = "brown1") +
            geom_line(aes_string(y = "T1"), col = "darkgoldenrod") +
            theme_minimal() +
            ylab("Temperature") + xlab("Date") +
            ggtitle(paste("Site: ", siteid, "; Tomst: ", ii, "; Time: ", times[tt])) +
            scale_x_datetime(date_minor_breaks = "1 day") -> GG1
          
          dft %>%
            ggplot(aes_string(x = "datetime")) +
            geom_line(aes_string(y = "moist"), col = "blue") +
            theme_minimal() +
            ylab("Moisture") + xlab("Date") +
            ggtitle(paste("Site: ", siteid, "; Time: ", times[tt])) +
            scale_x_datetime(date_minor_breaks = "1 day") -> GG2
          
          print(plot_grid(plotlist = list(GG1, GG2), nrow = 2))
        }
      }
    }
  } else {
    for (tt in 1:(length(times) - 1)) {
      
      temp %>% 
        filter(datetime >= ymd(times[tt]),
               datetime < ymd(times[tt + 1])) -> dft
      
      if(nrow(dft %>% filter(complete.cases(.)) > 0)){
        dft %>%
          ggplot(aes_string(x = "datetime")) +
          geom_line(aes_string(y = "T3"), col = "cornflowerblue") +
          geom_line(aes_string(y = "T2"), col = "brown1") +
          geom_line(aes_string(y = "T1"), col = "darkgoldenrod") +
          theme_minimal() +
          ylab("Temperature") + xlab("Date") +
          ggtitle(paste("Site: ", siteid, "; Time: ", times[tt])) +
          scale_x_datetime(date_minor_breaks = "1 day") -> GG1
        
        dft %>%
          ggplot(aes_string(x = "datetime")) +
          geom_line(aes_string(y = "moist"), col = "blue") +
          theme_minimal() +
          ylab("Moisture") + xlab("Date") +
          ggtitle(paste("Site: ", siteid, "; Time: ", times[tt])) +
          scale_x_datetime(date_minor_breaks = "1 day") -> GG2
        
        print(plot_grid(plotlist = list(GG1, GG2), nrow = 2))
      }
    }
  }
  dev.off()
}  

#################################################################################
# Screening each site for possible errors

# SITE = K01
siteid <- "K01"

office <- c(as_date(as_date(min(df$datetime)):as_date("2020-06-17")))
probls <- c()
hattu <- c()

df2 %>% mutate(probl = ifelse(site == siteid &
                                date %in% office,
                              2, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% probls,
                        1, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% hattu,
                        3, probl)) -> df2

# SITE = K02
siteid <- "K02"

office <- c(as_date(as_date(min(df$datetime)):as_date("2020-06-12")))
probls <- c()
hattu <- c()

df2 %>% mutate(probl = ifelse(site == siteid &
                                date %in% office,
                              2, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% probls,
                        1, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% hattu,
                        3, probl)) -> df2

# SITE = K03
siteid <- "K03"

office <- c(as_date(as_date(min(df$datetime)):as_date("2020-06-12")))
probls <- c()
hattu <- c()

df2 %>% mutate(probl = ifelse(site == siteid &
                                date %in% office,
                              2, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% probls,
                        1, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% hattu,
                        3, probl)) -> df2

# SITE = K04
siteid <- "K04"

office <- c(as_date(as_date(min(df$datetime)):as_date("2020-06-18")))
probls <- c()
hattu <- c()

df2 %>% mutate(probl = ifelse(site == siteid &
                                date %in% office,
                              2, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% probls,
                        1, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% hattu,
                        3, probl)) -> df2

# SITE = K05
siteid <- "K05"

office <- c(as_date(as_date(min(df$datetime)):as_date("2020-06-13")))
probls <- c()
hattu <- c()

df2 %>% mutate(probl = ifelse(site == siteid &
                                date %in% office,
                              2, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% probls,
                        1, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% hattu,
                        3, probl)) -> df2

# SITE = K06
siteid <- "K06"

office <- c(as_date(as_date(min(df$datetime)):as_date("2020-06-19")))
probls <- c()
hattu <- c()

df2 %>% mutate(probl = ifelse(site == siteid &
                                date %in% office,
                              2, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% probls,
                        1, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% hattu,
                        3, probl)) -> df2

# SITE = K07
siteid <- "K07"

office <- c(as_date(as_date(min(df$datetime)):as_date("2020-06-12")))
probls <- c()
hattu <- c()

df2 %>% mutate(probl = ifelse(site == siteid &
                                date %in% office,
                              2, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% probls,
                        1, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% hattu,
                        3, probl)) -> df2

# SITE = K08
siteid <- "K08"

office <- c(as_date(as_date(min(df$datetime)):as_date("2020-06-19")))
probls <- c()
hattu <- c()

df2 %>% mutate(probl = ifelse(site == siteid &
                                date %in% office,
                              2, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% probls,
                        1, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% hattu,
                        3, probl)) -> df2

# SITE = K09
siteid <- "K09"

office <- c(as_date(as_date(min(df$datetime)):as_date("2020-06-16")))
probls <- c()
hattu <- c()

df2 %>% mutate(probl = ifelse(site == siteid &
                                date %in% office,
                              2, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% probls,
                        1, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% hattu,
                        3, probl)) -> df2

# SITE = K10
siteid <- "K10"

office <- c(as_date(as_date(min(df$datetime)):as_date("2020-06-12")))
probls <- c()
hattu <- c()

df2 %>% mutate(probl = ifelse(site == siteid &
                                date %in% office,
                              2, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% probls,
                        1, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% hattu,
                        3, probl)) -> df2

# SITE = K11
siteid <- "K11"

office <- c(as_date(as_date(min(df$datetime)):as_date("2020-06-12")))
probls <- c()
hattu <- c()

df2 %>% mutate(probl = ifelse(site == siteid &
                                date %in% office,
                              2, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% probls,
                        1, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% hattu,
                        3, probl)) -> df2

# SITE = K12
siteid <- "K12"

office <- c(as_date(as_date(min(df$datetime)):as_date("2020-06-15")))
probls <- c()
hattu <- c()

df2 %>% mutate(probl = ifelse(site == siteid &
                                date %in% office,
                              2, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% probls,
                        1, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% hattu,
                        3, probl)) -> df2

# SITE = K13
siteid <- "K13"

office <- c(as_date(as_date(min(df$datetime)):as_date("2020-06-18")))
probls <- c()
hattu <- c()

df2 %>% mutate(probl = ifelse(site == siteid &
                                date %in% office,
                              2, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% probls,
                        1, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% hattu,
                        3, probl)) -> df2

# SITE = K14
siteid <- "K14"

office <- c(as_date(as_date(min(df$datetime)):as_date("2020-06-12")))
probls <- c()
hattu <- c()

df2 %>% mutate(probl = ifelse(site == siteid &
                                date %in% office,
                              2, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% probls,
                        1, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% hattu,
                        3, probl)) -> df2

# SITE = K15
siteid <- "K15"

office <- c(as_date(as_date(min(df$datetime)):as_date("2020-06-12")))
probls <- c()
hattu <- c()

df2 %>% mutate(probl = ifelse(site == siteid &
                                date %in% office,
                              2, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% probls,
                        1, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% hattu,
                        3, probl)) -> df2

# SITE = K16
siteid <- "K16"

office <- c(as_date(as_date(min(df$datetime)):as_date("2020-06-17")))
probls <- c()
hattu <- c()

df2 %>% mutate(probl = ifelse(site == siteid &
                                date %in% office,
                              2, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% probls,
                        1, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% hattu,
                        3, probl)) -> df2

# SITE = K17
siteid <- "K17"

office <- c(as_date(as_date(min(df$datetime)):as_date("2020-06-19")))
probls <- c()
hattu <- c()

df2 %>% mutate(probl = ifelse(site == siteid &
                                date %in% office,
                              2, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% probls,
                        1, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% hattu,
                        3, probl)) -> df2

# SITE = K18
siteid <- "K18"

office <- c(as_date(as_date(min(df$datetime)):as_date("2020-06-12")))
probls <- c()
hattu <- c()

df2 %>% mutate(probl = ifelse(site == siteid &
                                date %in% office,
                              2, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% probls,
                        1, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% hattu,
                        3, probl)) -> df2

# SITE = K19
siteid <- "K19"

office <- c(as_date(as_date(min(df$datetime)):as_date("2020-06-18")))
probls <- c()
hattu <- c()

df2 %>% mutate(probl = ifelse(site == siteid &
                                date %in% office,
                              2, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% probls,
                        1, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% hattu,
                        3, probl)) -> df2

# SITE = K20
siteid <- "K20"

office <- c(as_date(as_date(min(df$datetime)):as_date("2020-06-12")))
probls <- c(as_date(as_date("2020-09-08"):as_date("2021-06-09")))
hattu <- c()

df2 %>% mutate(probl = ifelse(site == siteid &
                                date %in% office,
                              2, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% probls,
                        1, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% hattu,
                        3, probl)) -> df2

# SITE = K21
siteid <- "K21"

office <- c(as_date(as_date(min(df$datetime)):as_date("2020-06-12")))
probls <- c()
hattu <- c()

df2 %>% mutate(probl = ifelse(site == siteid &
                                date %in% office,
                              2, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% probls,
                        1, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% hattu,
                        3, probl)) -> df2

# SITE = K22
siteid <- "K22"

office <- c(as_date(as_date(min(df$datetime)):as_date("2020-06-16")))
probls <- c()
hattu <- c()

df2 %>% mutate(probl = ifelse(site == siteid &
                                date %in% office,
                              2, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% probls,
                        1, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% hattu,
                        3, probl)) -> df2

# SITE = K23
siteid <- "K23"

office <- c(as_date(as_date(min(df$datetime)):as_date("2020-06-12")))
probls <- c()
hattu <- c()

df2 %>% mutate(probl = ifelse(site == siteid &
                                date %in% office,
                              2, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% probls,
                        1, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% hattu,
                        3, probl)) -> df2

# SITE = K24
siteid <- "K24"

office <- c(as_date(as_date(min(df$datetime)):as_date("2020-06-15")))
probls <- c()
hattu <- c()

df2 %>% mutate(probl = ifelse(site == siteid &
                                date %in% office,
                              2, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% probls,
                        1, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% hattu,
                        3, probl)) -> df2

# SITE = K25
siteid <- "K25"

office <- c(as_date(as_date(min(df$datetime)):as_date("2020-06-12")))
probls <- c()
hattu <- c()

df2 %>% mutate(probl = ifelse(site == siteid &
                                date %in% office,
                              2, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% probls,
                        1, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% hattu,
                        3, probl)) -> df2

# SITE = K26
siteid <- "K26"

office <- c(as_date(as_date(min(df$datetime)):as_date("2020-06-16")))
probls <- c(as_date(as_date("2021-08-07"):as_date("2021-08-12")))
hattu <- c()

df2 %>% mutate(probl = ifelse(site == siteid &
                                date %in% office,
                              2, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% probls,
                        1, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% hattu,
                        3, probl)) -> df2

# SITE = K27
siteid <- "K27"

office <- c(as_date(as_date(min(df$datetime)):as_date("2020-06-12")))
probls <- c(as_date(as_date("2020-06-12"):as_date("2020-08-24")))
hattu <- c()
probls2 <- c(as_date(as_date("2020-08-25"):as_date(max(df$datetime))))

df2 %>% mutate(probl = ifelse(site == siteid &
                                date %in% office,
                              2, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% probls,
                        1, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% hattu,
                        3, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% probls2,
                        11, probl)) -> df2

# SITE = K28
siteid <- "K28"

office <- c(as_date(as_date(min(df$datetime)):as_date("2020-06-15")))
probls <- c()
hattu <- c()

df2 %>% mutate(probl = ifelse(site == siteid &
                                date %in% office,
                              2, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% probls,
                        1, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% hattu,
                        3, probl)) -> df2

# SITE = K29
siteid <- "K29"

office <- c(as_date(as_date(min(df$datetime)):as_date("2020-06-12")))
probls <- c()
hattu <- c()

df2 %>% mutate(probl = ifelse(site == siteid &
                                date %in% office,
                              2, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% probls,
                        1, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% hattu,
                        3, probl)) -> df2

# SITE = K30
siteid <- "K30"

office <- c(as_date(as_date(min(df$datetime)):as_date("2020-06-16")))
probls <- c()
hattu <- c()

df2 %>% mutate(probl = ifelse(site == siteid &
                                date %in% office,
                              2, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% probls,
                        1, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% hattu,
                        3, probl)) -> df2

# SITE = K31
siteid <- "K31"

office <- c(as_date(as_date(min(df$datetime)):as_date("2020-06-17")))
probls <- c()
hattu <- c()

df2 %>% mutate(probl = ifelse(site == siteid &
                                date %in% office,
                              2, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% probls,
                        1, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% hattu,
                        3, probl)) -> df2

# SITE = K32
siteid <- "K32"

office <- c(as_date(as_date(min(df$datetime)):as_date("2020-06-12")))
probls <- c()
hattu <- c()

df2 %>% mutate(probl = ifelse(site == siteid &
                                date %in% office,
                              2, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% probls,
                        1, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% hattu,
                        3, probl)) -> df2

# SITE = K33
siteid <- "K33"

office <- c(as_date(as_date(min(df$datetime)):as_date("2020-06-12")))
probls <- c()
hattu <- c()

df2 %>% mutate(probl = ifelse(site == siteid &
                                date %in% office,
                              2, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% probls,
                        1, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% hattu,
                        3, probl)) -> df2

# SITE = K34
siteid <- "K34"

office <- c(as_date(as_date(min(df$datetime)):as_date("2020-06-19")))
probls <- c()
hattu <- c()

df2 %>% mutate(probl = ifelse(site == siteid &
                                date %in% office,
                              2, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% probls,
                        1, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% hattu,
                        3, probl)) -> df2

# SITE = K35
siteid <- "K35"

office <- c(as_date(as_date(min(df$datetime)):as_date("2020-06-12")))
probls <- c()
hattu <- c()

df2 %>% mutate(probl = ifelse(site == siteid &
                                date %in% office,
                              2, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% probls,
                        1, probl)) %>% 
  mutate(probl = ifelse(site == siteid &
                          date %in% hattu,
                        3, probl)) -> df2

########################################################################
# FILL MISSING TIMESTAMPS WITH NA

df2 <- df2 %>% 
  mutate(site_id = paste0(site, "_", tomst_id))
df <- df %>% 
  mutate(site_id = paste0(site, "_", tomst_id))

sites2 <- unique(df$site_id)
sites2 <- sites2[!grepl("_NA", sites2)]

df3 <- data.frame()
for(i in sites2){
  #i <- "AIL152_94194008"
  
  df %>% filter(site_id == i) -> temp
  
  temp %>% mutate(timediff = as.numeric(datetime - lag(datetime))) -> temp
  temp[1,"timediff"] <- 15
  holes <- table(temp$timediff)
  
  if(max(temp$timediff, na.rm = T) > 15){
    
    print(i)
    
    missingt <- c()
    for(ii in which(temp %>% pull(timediff) > 15)){
      
      temp %>% slice((ii-1):(ii+1)) %>% pull(timediff) -> diffs
      
      if(diffs[1] %% 15 == 0L){
        seq(temp %>% slice(ii-1) %>% pull(datetime),
            temp %>% slice(ii) %>% pull(datetime), by = "15 mins") -> seqs
        
        missingt <- c(missingt, 
                      as.character(seqs[which(!seqs %in% (temp %>% slice((ii-1):(ii+1)) %>% pull(datetime)))]))
        
      } else {
        seq(temp %>% slice(ii-1) %>% pull(datetime),
            temp %>% slice(ii) %>% pull(datetime), by = "10 mins") -> seqs
        
        missingt <- c(missingt, 
                      as.character(seqs[which(!seqs %in% (temp %>% slice((ii-1):(ii+1)) %>% pull(datetime)))]))
        
      }
    }
    
    missingdf <- data.frame(datetime = ymd_hms(missingt),
                            site_id = i)
    
    print(NROW(missingdf))
    
    temp %>% full_join(., missingdf, by = c("datetime", "site_id")) %>% 
      arrange(datetime) %>% 
      select(-timediff) -> temp
    
    df3 <- bind_rows(df3, temp)
    
  } else {
    
    temp %>% select(-timediff) -> temp
    
    df3 <- bind_rows(df3, temp)
  }
}

###################################################################################
# Delete erroneous data
#

df3 %>% mutate(date = as_date(datetime)) %>%
  left_join(., df2 %>% select(site, date, probl)) %>% 
  filter(probl != 2) -> df3

df3 %>% mutate(h = hour(datetime)) %>% 
  group_by(site, date) %>% 
  summarise(nh = length(unique(h))) %>% 
  filter(nh != 24) -> missh

df3 %>% left_join(., missh) %>% 
  mutate(T1 = replace(T1, !is.na(nh), NA),
         T2 = replace(T2, !is.na(nh), NA),
         T3 = replace(T3, !is.na(nh), NA),
         moist = replace(moist, !is.na(nh), NA)) %>%
  select(-nh,-zone) -> df3

##############################################################
# DETECT ANOMALIES CROSS-RELATING THE SITES

my_sd = function(x) {
  if(length(x) %% 2 == 0L) { return(sd(x, na.rm = T)) }
  if(length(x) %% 2 == 1L) { if(length(x) == 3){
    return(sd(x[-ceiling(0.5*length(x))]))
  } else {
    mid <- ceiling(0.5*length(x))
    return(sd(x[-c(mid-1, mid, mid+1)], na.rm = T))
  }  }
}
my_mean = function(x) {
  if(length(x) %% 2 == 0L) { return(mean(x)) }
  if(length(x) %% 2 == 1L) {
    if(length(x) == 3){
      return(mean(x[-ceiling(0.5*length(x))]))
    } else {
      if(length(x) == 5){
        mid <- ceiling(0.5*length(x))
        return(mean(x[-c(mid-1, mid, mid+1)]))
      } else {
        mid <- ceiling(0.5*length(x))
        return(mean(x[-c(mid-2, mid-1, mid, mid+1, mid+2)]))
      }
    }
  }
}

df3 %>% mutate(my = paste0(year(date),"_",month(date))) -> df3

dfall <- data.frame()
pdf("visuals/Temperature_graphs_spikes.pdf", 10, 12)
for(i in sites){
  #i <- 12
  
  print(i)
  df3 %>% filter(site == i) %>%
    filter(probl != 1) %>%
    mutate(timediff1 = as.numeric(datetime - lag(datetime)),
           timediff2 = as.numeric(lead(datetime) - datetime)) %>%
    filter(timediff1 %in% c(10,15)|timediff2 %in% c(10,15)) %>%
    mutate(timediff1 = as.numeric(datetime - lag(datetime)),
           timediff2 = as.numeric(lead(datetime) - datetime)) -> temp
  
  if(temp %>% pull(timediff1) %>% min(., na.rm = T) < 10){
    temp %>% filter(timediff1 + timediff2 >= 20 | is.na(timediff1 + timediff2)) -> temp
  }
  
  dftemp <- data.frame()
  for(ii in unique(temp$my)){
    # ii <- "2019_9"
    print(ii)
    
    ####################################
    # T1
    
    temp %>% filter(my == ii) %>%
      filter(probl %in% c(0,3,5:8)) %>% 
      select(datetime, T1, site) %>%
      filter(complete.cases(.)) %>%
      rename(T1f = T1,
             site2 = site) -> temp1
    
    if(length(unique(as_date(temp1$datetime))) > 3 & NROW(temp1) > 100){
      
      # Cross correlate site to others and extract the site with highest correlation
      temp1 %>%
        left_join(., df3, by = "datetime") %>%
        filter(site != i) %>%
        mutate(T1 = ifelse(probl %in% c(0,3,5:8), T1, NA)) %>% 
        arrange(site, datetime) %>% 
        group_by(site) %>% 
        summarise(cor = cor(T1f, T1)) %>% 
        arrange(desc(cor)) %>% 
        slice(1) %>% pull(site) -> site_to_compare
      
      temp1 %>%
        left_join(., df3 %>% filter(site == site_to_compare), by = "datetime") %>%
        select(datetime, T1f, T1) %>% 
        mutate(me = T1f-T1) %>%
        mutate(lag_T1 = me - lag(me)) %>%
        mutate(lead_T1 = me - lead(me)) %>%
        mutate(lag_T1f = T1f - lag(T1f)) %>%
        mutate(lead_T1f = T1f - lead(T1f)) %>%
        rowwise() %>% mutate(diff = sum(c(abs(lag_T1),abs(lead_T1)), na.rm = T),
                             change = sum(c(abs(lag_T1f),abs(lead_T1f)), na.rm = T)) %>% 
        ungroup() %>% 
        mutate(roll_diff = rollapply(diff, width=48*3, FUN=my_mean, fill = NA, partial = T, align = "center")+0.001,
               change_diff = rollapply(change, width=48*3, FUN=my_mean, fill = NA, partial = T, align = "center")+0.001) %>% 
        mutate(rel_diff = diff/roll_diff,
               rel_change = change/change_diff) %>%
        mutate(rel_lag_T1f = lag_T1f/change_diff,
               rel_lead_T1f = lead_T1f/change_diff) %>%
        mutate(error = ifelse((rel_diff > 10 & rel_change > 10 & rel_lag_T1f > 10 & rel_lead_T1f > 5), 1, 0)) %>% 
        mutate(error = ifelse((rel_diff > 10 & rel_change > 10 & rel_lag_T1f < (-10) & rel_lead_T1f < (-5)), 1, error)) %>% 
        mutate(error = ifelse((rel_diff > 20 & rel_change > 20 & rel_lag_T1f > 10), 1, error)) %>%
        mutate(error = ifelse((rel_diff > 20 & rel_change > 20 & rel_lead_T1f < (-10)), 1, error)) %>%
        mutate(error = ifelse((lag_T1f > 3 & lead_T1f > 3), 1, error)) %>% 
        mutate(error = ifelse((lag_T1f < (-3) & lead_T1f < (-3)), 1, error)) %>% 
        mutate(error = ifelse((diff < 1), 0, error)) %>%
        mutate(error = ifelse((diff > 20 | change > 20), 1, error)) %>%
        mutate(error = error + lag(error,1) + lag(error,2) + lag(error,3) + lag(error,4) + lead(error)) %>%
        mutate(error = ifelse(error > 1, 1, error)) %>%
        mutate(T1c = ifelse(!is.na(error) & error == 1, NA, T1f)) %>% 
        mutate(me = ifelse(!is.na(error) & error == 1, NA, me)) %>% 
        mutate(fill = rollapply(me, width=11, FUN=mean, na.rm = T, fill = NA, partial = T, align = "center") + T1) %>% 
        mutate(T1c = ifelse(is.na(T1c) & !is.na(error) & error == 1, fill, T1c)) -> temp1
      
      temp1 %>% 
        mutate(error = as.numeric(ifelse(error == 1, 0, NA))) %>% 
        ggplot(aes_string(x="datetime")) +
        geom_line(aes(y = T1f), col = "cornflowerblue") +
        geom_point(aes(y = error), col = "black") +
        theme_minimal() +
        ylab("T1") + xlab("Date")+
        scale_y_continuous(limits = c(ifelse(min(temp1$T1f, na.rm = T) > 0, (-0.1), min(temp1$T1f, na.rm = T)),
                                      ifelse(max(temp1$T1f, na.rm = T) < 0, 0.1, max(temp1$T1f, na.rm = T)))) +
        ggtitle(paste(i, ii)) -> GG1
      
      temp1 %>% 
        mutate(error = as.numeric(ifelse(error == 1, 0, NA))) %>% 
        ggplot(aes_string(x="datetime")) +
        geom_line(aes(y = T1c), col = "brown1") +
        geom_point(aes(y = error), col = "black") +
        theme_minimal() +
        ylab("T1 corrected") + xlab("Date")+
        scale_y_continuous(limits = c(ifelse(min(temp1$T1f, na.rm = T) > 0, (-0.1), min(temp1$T1f, na.rm = T)),
                                      ifelse(max(temp1$T1f, na.rm = T) < 0, 0.1, max(temp1$T1f, na.rm = T)))) +
        ggtitle(paste(i, ii)) -> GG2
      
    } else {
      temp1 %>% mutate(T1c = T1f) -> temp1
      
      temp1 %>% ggplot(aes_string(x="datetime")) +
        theme_minimal() +
        ggtitle(paste(i, ii, ", no correction")) -> GG1 -> GG2
    }
    
    ####################################
    # T2
    
    temp %>% filter(my == ii) %>%
      filter(probl %in% c(0,3,4:6,9)) %>% 
      select(datetime, T2, site) %>%
      filter(complete.cases(.)) %>%
      rename(T2f = T2,
             site2 = site) -> temp2
    
    if(length(unique(as_date(temp2$datetime))) > 3 & NROW(temp2) > 100){
      
      # Cross correlate site to others and extract the site with highest correlation
      temp2 %>%
        left_join(., df3, by = "datetime") %>%
        filter(site != i) %>%
        mutate(T2 = ifelse(probl %in% c(0,3,4:6,9), T2, NA)) %>% 
        arrange(site, datetime) %>% 
        group_by(site) %>% 
        summarise(cor = cor(T2f, T2)) %>% 
        arrange(desc(cor)) %>% 
        slice(1) %>% pull(site) -> site_to_compare
      
      temp2 %>%
        left_join(., df3 %>% filter(site == site_to_compare), by = "datetime") %>%
        select(datetime, T2f, T2) %>% 
        mutate(me = T2f-T2) %>%
        mutate(lag_T2 = me - lag(me)) %>%
        mutate(lead_T2 = me - lead(me)) %>%
        mutate(lag_T2f = T2f - lag(T2f)) %>%
        mutate(lead_T2f = T2f - lead(T2f)) %>%
        rowwise() %>% mutate(diff = sum(c(abs(lag_T2),abs(lead_T2)), na.rm = T),
                             change = sum(c(abs(lag_T2f),abs(lead_T2f)), na.rm = T)) %>% 
        ungroup() %>% 
        mutate(roll_diff = rollapply(diff, width=48*3, FUN=my_mean, fill = NA, partial = T, align = "center")+0.001,
               change_diff = rollapply(change, width=48*3, FUN=my_mean, fill = NA, partial = T, align = "center")+0.001) %>% 
        mutate(rel_diff = diff/roll_diff,
               rel_change = change/change_diff) %>%
        mutate(rel_lag_T2f = lag_T2f/change_diff,
               rel_lead_T2f = lead_T2f/change_diff) %>%
        mutate(error = ifelse((rel_diff > 10 & rel_change > 10 & rel_lag_T2f > 10 & rel_lead_T2f > 5), 1, 0)) %>% 
        mutate(error = ifelse((rel_diff > 10 & rel_change > 10 & rel_lag_T2f < (-10) & rel_lead_T2f < (-5)), 1, error)) %>% 
        mutate(error = ifelse((rel_diff > 15 & rel_change > 15 & rel_lag_T2f > 7 & rel_lead_T2f > 5), 1, 0)) %>% 
        mutate(error = ifelse((rel_diff > 15 & rel_change > 15 & rel_lag_T2f < (-7) & rel_lead_T2f < (-5)), 1, error)) %>% 
        mutate(error = ifelse((rel_diff > 20 & rel_change > 20 & rel_lag_T2f > 10), 1, error)) %>%
        mutate(error = ifelse((rel_diff > 20 & rel_change > 20 & rel_lead_T2f < (-10)), 1, error)) %>%
        mutate(error = ifelse((lag_T2f > 5 & lead_T2f > 5), 1, error)) %>% 
        mutate(error = ifelse((lag_T2f < (-5) & lead_T2f < (-5)), 1, error)) %>% 
        mutate(error = ifelse((diff < 1), 0, error)) %>%
        mutate(error = ifelse((diff > 20 | change > 20), 1, error)) %>%
        mutate(error = error + lag(error,1) + lag(error,2) + lag(error,3) + lag(error,4) + lead(error)) %>%
        mutate(error = ifelse(error > 1, 1, error)) %>%
        mutate(T2c = ifelse(!is.na(error) & error == 1, NA, T2f)) %>% 
        mutate(me = ifelse(!is.na(error) & error == 1, NA, me)) %>% 
        mutate(fill = rollapply(me, width=11, FUN=mean, na.rm = T, fill = NA, partial = T, align = "center") + T2) %>% 
        mutate(T2c = ifelse(is.na(T2c) & !is.na(error) & error == 1, fill, T2c)) -> temp2
      # edit(temp2 %>% mutate(datetime = as.character(datetime)))
      temp2 %>% 
        mutate(error = as.numeric(ifelse(error == 1, 0, NA))) %>% 
        ggplot(aes_string(x="datetime")) +
        geom_line(aes(y = T2f), col = "cornflowerblue") +
        geom_point(aes(y = error), col = "black") +
        theme_minimal() +
        ylab("T2") + xlab("Date")+
        scale_y_continuous(limits = c(ifelse(min(temp2$T2f, na.rm = T) > 0, (-0.1), min(temp2$T2f, na.rm = T)),
                                      ifelse(max(temp2$T2f, na.rm = T) < 0, 0.1, max(temp2$T2f, na.rm = T)))) +
        ggtitle(paste(i, ii)) -> GG3
      
      temp2 %>% 
        mutate(error = as.numeric(ifelse(error == 1, 0, NA))) %>% 
        ggplot(aes_string(x="datetime")) +
        geom_line(aes(y = T2c), col = "brown1") +
        geom_point(aes(y = error), col = "black") +
        theme_minimal() +
        ylab("T2 corrected") + xlab("Date")+
        scale_y_continuous(limits = c(ifelse(min(temp2$T2f, na.rm = T) > 0, (-0.1), min(temp2$T2f, na.rm = T)),
                                      ifelse(max(temp2$T2f, na.rm = T) < 0, 0.1, max(temp2$T2f, na.rm = T)))) +
        ggtitle(paste(i, ii)) -> GG4
      
    } else {
      temp2 %>% mutate(T2c = T2f) -> temp2
      
      temp2 %>% ggplot(aes_string(x="datetime")) +
        theme_minimal() +
        ggtitle(paste(i, ii, ", no correction")) -> GG3 -> GG4
    }
    
    ####################################
    # T3
    
    temp %>% filter(my == ii) %>%
      filter(probl %in% c(0,6,9,11)) %>% 
      select(datetime, T3, site) %>%
      filter(complete.cases(.)) %>%
      rename(T3f = T3,
             site2 = site) -> temp3
    
    if(length(unique(as_date(temp3$datetime))) > 3 & NROW(temp3) > 100){
      
      # Cross correlate site to others and extract the site with highest correlation
      temp3 %>%
        left_join(., df3, by = "datetime") %>%
        filter(site != i) %>%
        mutate(T3 = ifelse(probl %in% c(0,3,6,9,11), T3, NA)) %>% 
        arrange(site, datetime) %>% 
        group_by(site) %>% 
        summarise(cor = cor(T3f, T3)) %>% 
        arrange(desc(cor)) %>% 
        slice(1) %>% pull(site) -> site_to_compare
      
      temp3 %>%
        left_join(., df3 %>% filter(site == site_to_compare), by = "datetime") %>%
        select(datetime, T3f, T3) %>% 
        mutate(me = T3f-T3) %>%
        mutate(lag_T3 = me - lag(me)) %>%
        mutate(lead_T3 = me - lead(me)) %>%
        mutate(lag_T3f = T3f - lag(T3f)) %>%
        mutate(lead_T3f = T3f - lead(T3f)) %>%
        rowwise() %>% mutate(diff = sum(c(abs(lag_T3),abs(lead_T3)), na.rm = T),
                             change = sum(c(abs(lag_T3f),abs(lead_T3f)), na.rm = T)) %>% 
        ungroup() %>% 
        mutate(roll_diff = rollapply(diff, width=48*3, FUN=my_mean, fill = NA, partial = T, align = "center")+0.001,
               change_diff = rollapply(change, width=48*3, FUN=my_mean, fill = NA, partial = T, align = "center")+0.001) %>% 
        mutate(rel_diff = diff/roll_diff,
               rel_change = change/change_diff) %>%
        mutate(rel_lag_T3f = lag_T3f/change_diff,
               rel_lead_T3f = lead_T3f/change_diff) %>%
        mutate(error = ifelse((rel_diff > 10 & rel_change > 10 & rel_lag_T3f > 10 & rel_lead_T3f > 5), 1, 0)) %>% 
        mutate(error = ifelse((rel_diff > 15 & rel_change > 15 & rel_lag_T3f < (-10) & rel_lead_T3f < (-5)), 1, error)) %>% 
        mutate(error = ifelse((rel_diff > 20 & rel_change > 20 & rel_lag_T3f > 10), 1, error)) %>%
        mutate(error = ifelse((rel_diff > 20 & rel_change > 20 & rel_lead_T3f < (-10)), 1, error)) %>%
        mutate(error = ifelse((lag_T3f > 6 & lead_T3f > 6), 1, error)) %>%
        mutate(error = ifelse((lag_T3f < (-6) & lead_T3f < (-6)), 1, error)) %>%
        mutate(error = ifelse((diff < 1), 0, error)) %>%
        mutate(error = ifelse((diff > 20 | change > 20), 1, error)) %>%
        mutate(error = error + lag(error,1) + lag(error,2) + lag(error,3) + lag(error,4) + lead(error)) %>%
        mutate(error = ifelse(error > 1, 1, error)) %>%
        mutate(T3c = ifelse(!is.na(error) & error == 1, NA, T3f)) %>% 
        mutate(me = ifelse(!is.na(error) & error == 1, NA, me)) %>% 
        mutate(fill = rollapply(me, width=11, FUN=mean, na.rm = T, fill = NA, partial = T, align = "center") + T3) %>% 
        mutate(T3c = ifelse(is.na(T3c) & !is.na(error) & error == 1, fill, T3c)) -> temp3
      
      temp3 %>% 
        mutate(error = as.numeric(ifelse(error == 1, 0, NA))) %>% 
        ggplot(aes_string(x="datetime")) +
        geom_line(aes(y = T3f), col = "cornflowerblue") +
        geom_point(aes(y = error), col = "black") +
        theme_minimal() +
        ylab("T3") + xlab("Date")+
        scale_y_continuous(limits = c(ifelse(min(temp3$T3f, na.rm = T) > 0, (-0.1), min(temp3$T3f, na.rm = T)),
                                      ifelse(max(temp3$T3f, na.rm = T) < 0, 0.1, max(temp3$T3f, na.rm = T)))) +
        ggtitle(paste(i, ii)) -> GG5
      
      temp3 %>% 
        mutate(error = as.numeric(ifelse(error == 1, 0, NA))) %>% 
        ggplot(aes_string(x="datetime")) +
        geom_line(aes(y = T3c), col = "brown1") +
        geom_point(aes(y = error), col = "black") +
        theme_minimal() +
        ylab("T3 corrected") + xlab("Date")+
        scale_y_continuous(limits = c(ifelse(min(temp3$T3f, na.rm = T) > 0, (-0.1), min(temp3$T3f, na.rm = T)),
                                      ifelse(max(temp3$T3f, na.rm = T) < 0, 0.1, max(temp3$T3f, na.rm = T)))) +
        ggtitle(paste(i, ii)) -> GG6
      
    } else {
      temp3 %>% mutate(T3c = T3f) -> temp3
      
      temp3 %>% ggplot(aes_string(x="datetime")) +
        theme_minimal() +
        ggtitle(paste(i, ii, ", no correction")) -> GG5 -> GG6
    }
    
    #############################################
    
    print(plot_grid(plotlist = list(GG1,GG2,GG3,GG4,GG5,GG6), nrow = 3))
    
    
    temp4 <- full_join(temp %>% filter(my == ii) %>% select(site, datetime, moist, date, probl),
                       temp1 %>% select(datetime, T1c) %>%
                         rename(T1 = T1c), by = "datetime")
    temp4 <- full_join(temp4,
                       temp2 %>% select(datetime, T2c) %>%
                         rename(T2 = T2c), by = "datetime")
    temp4 <- full_join(temp4,
                       temp3 %>% select(datetime, T3c) %>%
                         rename(T3 = T3c), by = "datetime")
    
    dftemp <- bind_rows(dftemp,temp4)
    
  }
  dfall <- bind_rows(dfall, dftemp)
}
dev.off()

####################################################################
# MASK IMPOSSIBLE VALUES

dfall %>% mutate(T1 = ifelse(T1 < (-50) | T1 > 50, NA, T1),
                 T2 = ifelse(T2 < (-50) | T2 > 50, NA, T2),
                 T3 = ifelse(T3 < (-50) | T3 > 50, NA, T3),
                 moist = ifelse(moist < 600 | moist >= 4096, NA, moist)) -> dfall

dfall <- dfall %>% 
  mutate(moist = ifelse(moist < 1000 & site == "K25", NA, moist))

###############################################################################
# PLOT CORRECTED

pdf("visuals/Temperature_graphs_corrected_Roslin.pdf", 12, 10)
for(i in sites){
  #i <- "L12
  print(i)
  dfall %>% filter(site == i) %>% 
    mutate(T1 = as.numeric(ifelse(probl %in% c(1,4,9,11), NA, T1))) %>% 
    mutate(T2 = as.numeric(ifelse(probl %in% c(1,7,8,11), NA, T2))) %>% 
    mutate(T3 = as.numeric(ifelse(probl %in% c(1,4,5,7,8), NA, T3))) %>% 
    #group_by(date) %>% 
    #summarise_at(vars(i, "soil"), funs(mean, min, max), na.rm = T) %>% 
    #lapply(function(x) replace(x, is.infinite(x),NA)) %>% as_tibble() %>% 
    ggplot(aes_string(x="datetime")) +
    geom_line(aes_string(y = "T3"), col = "cornflowerblue") +
    geom_line(aes_string(y = "T2"), col = "brown1") +
    geom_line(aes_string(y = "T1"), col = "darkgoldenrod") +
    theme_minimal() +
    ylab("Temperature") + xlab("Date")+
    scale_y_continuous(limits = c(-20, 35))+
    ggtitle(i) -> GG1
  
  dfall %>% filter(site == i) %>% 
    mutate(moist = as.numeric(ifelse(probl %in% c(1,6,8,9), NA, moist))) %>% 
    mutate(moist = as.numeric(ifelse(T1 <= 1, NA, moist))) %>% 
    #group_by(date) %>% 
    #summarise_at(vars(i, "soil"), funs(mean, min, max), na.rm = T) %>% 
    #lapply(function(x) replace(x, is.infinite(x),NA)) %>% as_tibble() %>% 
    ggplot(aes_string(x="datetime")) +
    geom_line(aes_string(y = "moist"), col = "black") +
    theme_minimal() +
    ylab("Soil moisture count") + xlab("Date")+
    scale_y_continuous(limits = c(500, 4000))+
    ggtitle(i) -> GG2
  
  print(plot_grid(plotlist = list(GG1, GG2), nrow = 2))
  
}
dev.off()


round2 <- function(x) round(x,2)

fwrite(dfall %>% select(-date) %>% 
         mutate(across(T1:T3, round2)), "data/tomst_data/Roslin_tomst_data.csv")

####################################################################################