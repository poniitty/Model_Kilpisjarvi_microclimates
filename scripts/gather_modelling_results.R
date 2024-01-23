library(tidyverse)
library(terra)

f <- list.files("output/model_summaries", "cvmetrics", full.names = T)
cv <- lapply(f, read_csv) %>% 
  bind_rows()
f <- list.files("output/model_summaries", "varimp", full.names = T)
vi <- lapply(f, read_csv) %>% 
  bind_rows()

cv %>% select(-c(mtry:min.node.size)) %>% arrange(Rsquared) %>% print(n = 20)
cv %>% select(-c(mtry:min.node.size)) %>% arrange(desc(Rsquared)) %>% print(n = 20)

mean(cv$Rsquared)
median(cv$Rsquared)

cv %>% 
  select(RMSE:MAE, resp_var) %>% 
  arrange(desc(Rsquared)) %>% 
  relocate(resp_var, Rsquared) %>% 
  mutate(Rsquared = round(Rsquared, 3),
         RMSE = round(RMSE, 2),
         MAE = round(MAE, 2)) %>% 
  view
  
vi1 <- vi %>% 
  group_by(predictor) %>% 
  count %>% 
  arrange(desc(n))

vi2 <- vi %>% 
  pivot_wider(id_cols = predictor, names_from = resp_var, values_from = varimp) %>% 
  mutate(across(everything(), ~ifelse(is.na(.x), 0, .x))) %>% 
  pivot_longer(cols = where(is.numeric), names_to = "resp_var", values_to = "varimp") %>% 
  group_by(resp_var) %>% 
  mutate(varimp = scales::rescale(varimp, c(0, 1))) %>% 
  arrange(resp_var, desc(varimp)) %>% 
  group_by(predictor) %>% 
  summarise(varimp = mean(varimp)) %>% 
  arrange(desc(varimp)) %>% 
  mutate(varimp = round(varimp, 3))

full_join(vi1, vi2) %>% 
  arrange(desc(n)) %>% 
  view

########################################
# Correlation matrix
# devtools::install_github("kassambara/ggcorrplot", lib = "/projappl/project_2007415/Rpackages/")
library(ggcorrplot, lib.loc = "/projappl/project_2007415/Rpackages/")

f <- list.files("output/predictions", pattern = ".tif$", full.names = T)

# for(i in f){
#   print(i)
#   r <- rast(i)
#   if(names(r) != gsub(".tif","",tail(str_split(i,"/")[[1]],1))){
#     r <- r*1
#     names(r) <- gsub(".tif","",tail(str_split(i,"/")[[1]],1))
#     writeRaster(r, i, overwrite = TRUE)
#     gc()
#   }
# }

s <- lapply(f, rast)
s <- rast(s)
names(s) <- gsub(".tif","",basename(f))

d <- spatSample(s, 10000, "random", na.rm = T)

corr <- round(cor(d), 2)

gg <- ggcorrplot(corr, hc.order = TRUE, show.diag = F,
           outline.col = "white")
ggsave("visuals/cormat.pdf", plot = gg, device = "pdf", height = 15, width = 15)

# Find best non-correlated set

f <- list.files("output/model_summaries", "cvmetrics", full.names = T)
cv <- lapply(f, read_csv) %>% 
  bind_rows()

d2 <- d
corr <- abs(cor(d2))
diag(corr) <- NA
success <- max(corr, na.rm = T) > 0.7
while(success){
  
  coln <- which(corr == max(corr, na.rm = T), arr.ind = TRUE) %>% rownames()
  
  nn <- cv %>% filter(resp_var %in% coln) %>% 
    arrange(Rsquared) %>% slice(1) %>% pull(resp_var)
  
  d2 <- d2 %>% select(-all_of(nn))
  
  corr <- abs(cor(d2))
  diag(corr) <- NA
  success <- max(corr, na.rm = T) > 0.7
}



#######################################################
# Are correlations preserved

dd <- read_csv("data/microclimate_variables.csv") %>% 
  select(-site) %>% 
  select(sort(names(.)))

d <- d %>% 
  select(names(dd))

corr1 <- cor(dd, use = "pairwise.complete.obs")
corr2 <- cor(d)

corrdiff <- corr2 - corr1
median(abs(corrdiff))
mean(abs(corrdiff))
gg <- ggcorrplot(corrdiff, show.diag = F,
                 outline.col = "white")
ggsave("visuals/cordiffmat.pdf", plot = gg, device = "pdf", height = 15, width = 15)


d2 <- dd
corr <- abs(cor(d2, use = "pairwise.complete.obs"))
diag(corr) <- NA
success <- max(corr, na.rm = T) > 0.5
while(success){
  
  coln <- which(corr == max(corr, na.rm = T), arr.ind = TRUE) %>% rownames()
  
  nn <- cv %>% filter(resp_var %in% coln) %>% 
    arrange(Rsquared) %>% slice(1) %>% pull(resp_var)
  
  d2 <- d2 %>% select(-all_of(nn))
  
  corr <- abs(cor(d2, use = "pairwise.complete.obs"))
  diag(corr) <- NA
  success <- max(corr, na.rm = T) > 0.5
}

