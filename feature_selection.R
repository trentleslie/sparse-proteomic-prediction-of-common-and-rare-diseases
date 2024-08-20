# Feature Selection Script

library(caret)
library(glmnet)
library(ROSE)
library(parallel)

perform_feature_selection <- function(u.train, dz, protein_list) {
  jj <- lapply(1:200, function(x) sample(nrow(u.train), round(nrow(u.train)*0.50), replace=F))
  
  las.morb <- mclapply(1:length(jj), function(i) {
    print(i)
    tmp <- train(u.train[jj[[i]], c(protein_list, colnames(rand.vars))],
                 u.train[jj[[i]], dz],
                 metric = "Kappa",
                 method = "glmnet",
                 family = "binomial",
                 tuneGrid = as.data.frame(expand_grid(
                   alpha = 1, lambda = 10^-seq(4, .25, -.5))),
                 trControl = trainControl(method = "repeatedcv",
                                          number = 3, repeats = 5,
                                          sampling = "rose",
                                          allowParallel = T)
    )
    cf <- as.matrix(coef(tmp$finalModel, s = tmp$finalModel$lambdaOpt))
    gc()
    return(cf)
  }, mc.cores = 12, mc.allow.recursive = T)
  
  return(las.morb)
}

generate_feature_selection_ranking <- function(las.morb) {
  p.select <- matrix(nrow = nrow(las.morb[[1]]), ncol = length(las.morb))
  
  for(i in 1:length(las.morb)) {
    tryCatch({
      tmp <- las.morb[[i]][,1]
      if(is.null(tmp)) {
        return(tmp <- rep(0, nrow(p.select)))
      } else {
        NULL
      }
    }, error = function(e) {
      NULL
    })
    p.select[,i] <- tmp
  }
  row.names(p.select) <- names(las.morb[[1]][,1])
  
  p.select <- p.select %>% 
    as_tibble(rownames = "mrc_olink.id") %>% 
    filter(mrc_olink.id != "(Intercept)") %>% 
    mutate(select = abs(rowSums(across(where(is.numeric))))) %>% 
    arrange(desc(select)) 
  
  p.select <- as.data.frame(p.select)
  p.select$select.perc <- (p.select$select) / max(p.select$select)
  
  return(p.select)
}