setwd("/Volumes/YZ")
library(tidyr)
update_merge <- read.delim("cohortdata.txt")

## 1. Prepare the covariate matrix
whenscreen_year_month <- separate( update_merge["first_clinic_date"],
                                   first_clinic_date,
                                   into = c("year", "month", "day", "hour", "min", "second"),
                                   sep = "[^[:alnum:]]+" )
whenscreen_year_month <- as.data.frame(sapply(whenscreen_year_month ,as.numeric))
update_merge_trimmed <- update_merge[whenscreen_year_month$year>2007,]
whenscreen_year_month <- subset(whenscreen_year_month, year>2007)

# Prepare exposure variable
## 3. Prepare exposure variable
pollution <- update_merge_trimmed[,4:243]
pollution_index <- data.frame("index" = colnames(pollution))
pollution_index <- separate(pollution_index, index, into = c("M", "Pollutant", "year"), sep = "[_.]")
pollution_index$index <- paste(pollution_index$M, pollution_index$Pollutant, pollution_index$year, sep = "_")

SID <- update_merge_trimmed$SID
create_exposurevar <- function(pollutant){ # define the exposure window
        colindex <- pollution_index$index
        Z <- integer()
        for (i in 1:nrow(whenscreen_year_month)){
                temp_year <- whenscreen_year_month$year[i]
                temp_index <- paste("M01", pollutant, temp_year, sep = "_")
                start <- match(temp_index, colindex)
                end <- start+11
                Z[i] <- rowMeans(pollution[i,start:end])
        }
        Z <- data.frame(Z, SID)
        return(Z)
}


load("betas_and_samples.rda")
library(glmnet)
Lasso_Pollutant <- function(pollutant){
        pollutiondata <- create_exposurevar(pollutant)
        # IQR normalisation
        iqr_normalise <- function(x) {
                (x - median(x, na.rm = TRUE)) / IQR(x, na.rm = TRUE)
        }
        pollutiondata$Z <- iqr_normalise(pollutiondata$Z)
        ## 4. merge to sampledata
        sampledata <- data.frame(y = pollutiondata$Z)
        rownames(sampledata) <- pollutiondata$SID
        
        ## 1. data preparation
        infile <- paste("QC_easyEWAS_lmer_ANAV", pollutant, "noanx_nobmi_train.txt", sep = "_")
        ewas <- read.delim(infile)
        ewas <- ewas[order(ewas$FDR, decreasing = F),]
        cutoff <- round(nrow(ewas)*0.01)
        cpgs <- na.omit(ewas$probe[1:cutoff])
        methylationdata <- t(betas[cpgs,])
        methylationdata <- methylationdata[intersect(rownames(sampledata), rownames(methylationdata)),]
        model_mtx <- cbind(sampledata, methylationdata)
        colnames(model_mtx)[1] <- "y"
        
        ## 2. train-test split
        test_idx  <- sample(nrow(model_mtx), 0.2 * nrow(model_mtx))
        train_idx <- setdiff(seq_len(nrow(model_mtx)), test_idx)
        
        mets_train <- model_mtx[train_idx,]
        mets_train <- makeX(mets_train, na.impute = TRUE) ## impute missing values
        
        mets_test <- model_mtx[test_idx,]
        mets_test <- makeX(mets_test, na.impute = TRUE)
        ## 3. model training
        library(glmnet)
        penalty.factor <- c(
                rep(1, length(cpgs))  # CpGs — LASSO penalized
        )
     
        model <- glmnet(x = as.matrix(mets_train[,-1]), y = mets_train[,1],
                        alpha = 1, penalty.factor = penalty.factor)
        
        lambda_use <- model$lambda[length(model$lambda)]  # smallest lambda in path
        coefs <- coef(model, s = lambda_use)
        
        # ---- generate predicted values ----
        train_pred <- as.vector(predict(model, newx = as.matrix(mets_train[,-1]), s = lambda_use))
        test_pred  <- as.vector(predict(model, newx = as.matrix(mets_test[,-1]),  s = lambda_use))
        
        # ---- train performance ----
        test_cor  <- cor(mets_test[,1], test_pred)
        test_r2   <- test_cor^2
        test_rmse <- sqrt(mean((mets_test[,1] - test_pred)^2))
        
        # ---- train performance ----
        train_cor  <- cor(mets_train[,1], train_pred)
        train_r2   <- train_cor^2
        train_rmse <- sqrt(mean((mets_train[,1] - train_pred)^2))
        
        res <- list("coefs" = coefs,
                    "test_cor" = test_cor,
                    "test_r2" = test_r2,
                    "test_rmse" = test_rmse,
                    "train_cor" = train_cor,
                    "train_r2"= train_r2,
                    "train_rmse" = train_rmse)
        return(res)
}


library(ggplot2)
coef_df <- rep(1, 7753)
test_cor <- integer()
test_r2 <- integer()
test_rmse <- integer()

train_r2 <- integer()
train_cor <- integer()
train_rmse <- integer()
for (i in 1:100){
        print(i)
        res <- Lasso_Pollutant("PM25")
        coef_df <- cbind(coef_df, res$coefs)
        test_cor <- c(test_cor, res$test_cor)
        test_r2 <- c(test_r2, res$test_r2)
        test_rmse <- c(test_rmse, res$test_rmse)
        train_cor <- c(train_cor, res$train_cor)
        train_r2 <- c(train_r2, res$train_r2)
        train_rmse <- c(train_rmse, res$train_rmse)
}

coef_df <- coef_df[,-1]
selection_freq <- rowMeans(coef_df != 0)
selection_freq <- sort(selection_freq, decreasing = TRUE)

freq_df <- data.frame(
        variable  = names(selection_freq),
        frequency = selection_freq
) 

write.table(freq_df, "Updated_Iteration_LASSO_PM25.txt", sep = "\t", quote = F, row.names = F)

library(ggplot2)
ggplot(freq_df, aes(x = reorder(variable, frequency), y = frequency)) +
        geom_bar(stat = "identity", fill = "steelblue", width = 0.7) +
        geom_hline(yintercept = 0.5, colour = "red", linetype = "dashed", linewidth = 0.8) +
        annotate("text", x = 1, y = 0.52, label = "50% threshold", colour = "red", size = 3, hjust = 0) +
        coord_flip() +
        labs(x = "CpG Sites", y = "Selection Frequency",
             title = "LASSO Feature Selection Stability across 100 Bootstrap Models") +
        theme_classic() +
        theme(axis.text.y = element_text(size = 7))


library(ggplot2)

# Sample data
df <- data.frame(test_cor,
                 test_r2,
                 test_rmse,
                 train_r2,
                 train_cor,
                 train_rmse)

ggplot(df, aes(x = "", y = test_rmse)) +
        geom_boxplot(outlier.shape = NA) +  # avoid double-plotting outliers
        geom_jitter(width = 0.1, alpha = 0.6, size = 1.5) +
        labs(x = NULL, y = "Value", title = "Test RMSE") +
        theme_minimal()
