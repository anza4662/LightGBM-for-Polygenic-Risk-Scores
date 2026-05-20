#!/bin/Rscript
library(data.table)
library(dplyr)
library(pROC)
library(ggplot2)
library(stringr)
library(lightgbm)
library(knitr)
library(forcats)

args = commandArgs(trailingOnly=TRUE)
dataPath <- args[1]
outDir <- args[2]

# Change paramters here
lrs <- c(0.2, 0.1, 0.001)
depths <- c(2, 4, 6, 8, 10)
featureFrac <- 0.8
maxRounds <- 10000
nThreads <- 36
earlyStopRounds <- 100
evalMetric <- "auc"

# Read in data
fullData <- fread(dataPath)

# Split data into train test and validation
train <- fullData %>% sample_frac(0.6)

temp <- anti_join(fullData, train, by = c("FID","IID"))

val <- temp %>% sample_frac(0.5)
test <- anti_join(temp, val, by = c("FID","IID"))

# Translate from plink's 1 and 2 binary enconding to 0 and 1 for cases and controls
train$PHENOTYPE <- train$PHENOTYPE - 1
val$PHENOTYPE <- val$PHENOTYPE - 1

# Phenotypes
trainY <- train$PHENOTYPE
valY <- val$PHENOTYPE

# Remove unwanted cols
trainX <- train %>% select(-c(PAT,MAT,SEX,PHENOTYPE,FID,IID))
testX <- test %>% select(-c(PAT,MAT,SEX,PHENOTYPE,FID,IID))
valX <- val %>% select(-c(PAT,MAT,SEX,PHENOTYPE,FID,IID))

# Fix colnames

# This is neccesary if the id's for the variants are names like "CHROM:POS", if you have other names for the variants
# you may need to change the code below since the LigthGBM package does not freely allow all symbols for the variable
# names.

originalColnames <- colnames(trainX)

varTrain <- str_remove(colnames(trainX), "_.*$") %>% str_remove(":")
varTest <- str_remove(colnames(testX), "_.*$") %>% str_remove(":")
varVal <- str_remove(colnames(valX), "_.*$") %>% str_remove(":")

colnames(testX) <- varTest
colnames(trainX) <- varTrain
colnames(valX) <- varVal

colnamesNewAndOrig <- data.table(origName = originalColnames, newName = varTrain)

# The lightgbm package want the data in a special format so we convert it to this format using lgb.Dataset

dtrain <-as.matrix(trainX)
dval <- as.matrix(valX)
dtest <- as.matrix(testX)

dtrain2 <- lgb.Dataset(data = dtrain,
                       label = trainY
                       )
dval2 <- lgb.Dataset(data = dval,
                     label = valY
                     )

watchlist <- list(val = dval2)

results <- c()

npar <- length(lrs) * length(depths)

cat(paste0("\n Running grid search over ", npar, " params."))
cat("\n Grid = ")
cnt <- 0

for (lr in lrs) {
  for (d in depths) {
    cat(paste0(as.integer(round(cnt/npar, 2)*100), "% "))
    cnt <- cnt + 1

    lgbFit <- lgb.train(
      data = dtrain2,
      nrounds = maxRounds,
      early_stopping_rounds = earlyStopRounds,
      valids = watchlist,
      params = list(
        task = "train",
        objective = "binary",
        learning_rate = lr,
        num_leaves = 2^d - 1,
        num_threads = nThreads,
        max_depth = d,
        metric = evalMetric,
        feature_fraction = featureFrac
      ),
      verbose = -1
    )

    tmp <- data.frame(depth = d, lr = lr,
                      M = lgbFit$best_iter,
                      auc = round(lgbFit$best_score, 5)
    )
    results <- rbind(results, tmp)

  }
}
cat("DONE\n")

write.table(results,
            paste0(outDir,"/validationResults.tab"),
            row.names = F,
            quote = F,
            col.names = T
            )

bestModel <- results[results$auc == max(results$auc),]

cat("\n \nBest hyperparameters:")
knitr::kable(bestModel, format = "simple")

cat("\n Fitting final model..")

lgbFinal <- lightgbm(
  data = dtrain,
  label = trainY,
  nrounds = bestModel$M,
  params = list(
    num_threads = nThreads,
    max_depth = bestModel$depth,
    objective = "binary",
    learning_rate = bestModel$lr
  ),
  verbose = -1
)

predictionsLGBM <- predict(lgbFinal, dtest, type = "response")

cat("\n \nAUC test data")
r1 <- roc(testY, predictionsLGB, ci = T)
r1

imp <- lgb.importance(lgbFinal, percentage = TRUE)

colnames(imp)[1] <- "newName"
imp <- merge(imp, colnamesNewAndOrig)

write.table(imp, "/nfs/GENETEC/PRSes/Scripts/MLmethods/lgbm/results/featureImportance.tab",
           quote = F, col.names = T, row.names = F)

top <- imp %>% top_n(n = 20, wt = Gain) %>% mutate(origName = fct_reorder(origName, Gain))

p1 <- ggplot(data = top, aes(x = origName, y = Gain)) + geom_col(fill = "lightblue") +
  theme_minimal() + coord_flip() +
  theme(axis.title.y = element_text(margin = margin(r = 10)),
        legend.position = "none") +
labs(
  title = "Feature importance (top 20 SNPs) for LightGBM",
  x = "SNP",
  y = "Gain"
)

ggsave(p1, paste0(outDir,"/featureImportance.png"))

cat("\nLightGBM tool finished.")