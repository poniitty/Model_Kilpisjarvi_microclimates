library(tidyverse)
library(caret)
library(doParallel)
library(terra)

cl <- makePSOCKcluster(future::availableCores())
registerDoParallel(cl)

d <- full_join(read_csv("data/microclimate_variables.csv"),
               read_csv("data/predictors/env_data.csv")) %>% 
  full_join(., read_csv("data/predictors/env_buf20m_data.csv")) %>% 
  full_join(., read_csv("data/predictors/planet_pca_data.csv")) %>% 
  mutate(dem = round(dem/10))

d <- d %>% select(-nearZeroVar(d))

correlationMatrix <- cor(d[,63:ncol(d)], use = "pairwise.complete.obs")
dim(correlationMatrix)
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.9, names = TRUE, exact = TRUE)
print(highlyCorrelated)

d <- d %>% select(-all_of(highlyCorrelated))
names(d)

resp_vars <- names(d)[2:64]
all_predictors <- names(d)[65:ncol(d)]

# Tuning parameters for Random Forests
tunegrid_ranger <- expand.grid(mtry=c(2:8), splitrule = "extratrees", min.node.size = c(1:5))


###################################################################
# Raster stack
s <- rast("data/predictors/scd.tif")
names(s) <- "scd"

pca <- rast("/scratch/project_2007415/Planet/KilpisjarviAll/extra/planet_PCA.tif")
names(pca) <- gsub("PC","PCA",names(pca))

s <- c(s, pca)

for(i in all_predictors[!all_predictors %in% names(s)]){
  print(i)
  if(grepl("20",i)){
    r <- rast(paste0("data/predictors/",gsub("20","",i),".tif"))
    r <- round(project(r, s))
    fm <- focalMat(r, 20, "circle")
    fm[fm > 0] <- 1
    r <- round(focal(r, fm, "mean", na.policy = "omit", na.rm = T))
  } else {
    r <- rast(paste0("data/predictors/",i,".tif"))
    r <- round(project(r, s[[1]]))
    if(max(minmax(r)) > 200000){
      r[r > 200000] <- NA
    }
  }
  names(r) <- i
  s <- c(s, r)
}
gc()
s[["dem"]] <- round(s[["dem"]]/10)
writeRaster(s, "data/predictors/predstack.tif", overwrite = T, datatype = "INT2S")
rm(s)

s <- rast("data/predictors/predstack.tif")

# Cross-validation options
cv_method <- "repeatedcv"
cv_number <- 10
cv_repeats <- 1

ctrl <- trainControl(
  method = cv_method,
  number = cv_number,
  repeats = cv_repeats,
  savePredictions = F
)

for(resp_var in resp_vars){
  # resp_var <- "moist_mean_july"
  if(!file.exists(paste0("output/predictions/", resp_var, ".tif"))){
    print(resp_var)
    modd <- d %>% select(all_of(c(resp_var, all_predictors))) %>% 
      drop_na()
    # modd %>% filter(!complete.cases(.))
    # define the control using a random forest selection function
    control <- rfeControl(functions=rfFuncs, method="cv", number=10, repeats = 4)
    results <- rfe(x = modd %>% select(-all_of(resp_var)), 
                   y = modd %>% pull(all_of(resp_var)),
                   sizes=c(3:length(all_predictors)), rfeControl=control)
    
    plot(results, type=c("g", "o"))
    
    modm <- train(x = modd %>% select(all_of(predictors(results))), 
                  y = modd %>% pull(all_of(resp_var)),
                  method = "ranger",
                  importance = "permutation",
                  trControl = ctrl,
                  tuneGrid=tunegrid_ranger,
                  verbose=FALSE,
                  metric = "RMSE",
                  na.action = "na.omit")
    
    varImp(modm, scale = F)$importance %>% 
      rownames_to_column("predictor") %>% 
      rename(varimp = Overall) %>% 
      mutate(resp_var = resp_var) %>% 
      write_csv(paste0("output/model_summaries/varimp_",resp_var,".csv"))
    
    modm$results %>% 
      arrange(RMSE) %>% 
      slice_head(n = 1) %>% 
      mutate(resp_var = resp_var) %>% 
      write_csv(paste0("output/model_summaries/cvmetrics_",resp_var,".csv"))
    
    # Spatial predictions
    system.time({
      pr1 <- predict(s, modm$finalModel, na.rm = T, type = "response")
    })
    
    gc()
    plot(pr1$prediction, main = resp_var)
    writeRaster(round(pr1$prediction,2), 
                paste0("output/predictions/", resp_var, ".tif"),
                overwrite = T)
    rm(pr1)
    tmpFiles(current = TRUE, remove = TRUE)
    unlink(list.files(tempdir(), recursive = T, full.names = T), force = T)
    gc()
  }
}
