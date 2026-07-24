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
covariates <- c("age_when_screened", "gender", "smoking3","edu_code", "hhincome_code", "drinker_category", "phy_code")
covs <- update_merge_trimmed[,c(covariates, "SID")]
covs$gender <- ifelse(covs$gender=="M", 1, 0)
covs$smoking3 <- as.numeric(factor(covs$smoking3, levels = c("never smoker", "former smoker", "current smoker")))
covs$edu_code <- as.numeric(factor(covs$edu_code))
covs$hhincome_code <- as.numeric(factor(covs$hhincome_code))
covs$drinker_category <- as.numeric(factor(covs$drinker_category))

## 2. Prepare the cell prop matrix
cell <- c("CD8.naive", "CD8pCD28nCD45RAn","NK", "PlasmaBlast", "Mono", "CD4T", "SID")
cellp <- update_merge_trimmed[,cell]

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

pollutiondata <- create_exposurevar("NO2")
# IQR normalisation
iqr_normalise <- function(x) {
        (x - median(x, na.rm = TRUE)) / IQR(x, na.rm = TRUE)
}
pollutiondata$Z <- iqr_normalise(pollutiondata$Z)
## 4. merge to sampledata
sampledata <- merge(pollutiondata, covs, by = "SID")
sampledata <- merge(sampledata,cellp, by = "SID")
rownames(sampledata) <- sampledata$SID
sampledata <- sampledata[,-1]
load("betas_and_samples.rda")
library(glmnet)
Lasso_Pollutant <- function(pollutant){
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
        set.seed(1234) # make sure train-test random split the same everytime as EWAS
        test_idx  <- sample(nrow(model_mtx), 0.2 * nrow(model_mtx))
        train_idx <- setdiff(seq_len(nrow(model_mtx)), test_idx)
        
        mets_train <- model_mtx[train_idx,]
        mets_train <- makeX(mets_train, na.impute = TRUE) ## impute missing values
        
        mets_test <- model_mtx[test_idx,]
        mets_test <- makeX(mets_test, na.impute = TRUE)
        ## 3. model training
        library(glmnet)
        penalty.factor <- c(
                rep(0, 15),  # covariates — never shrunk out // factors are treated as individual binary features
                rep(1, length(cpgs))  # CpGs — LASSO penalized
        )
        
        Nbootstrap <- 500
        ensemble_model_list <- list()
        for (i in 1:Nbootstrap){
                print(paste("we are now running bootstrap", i, sep = " "))
                select_id <- sample.int(nrow(mets_train) ,replace = T)
                mets_train_sample <- mets_train[select_id,]
                glmnet_Modelfit <- cv.glmnet(x = as.matrix(mets_train_sample[,-1]),  y= mets_train_sample[,1], 
                                             alpha = 1,
                                             nfolds = 10,
                                             penalty.factor = penalty.factor, ## whether to shrunk out features (covariates no need for this)
                                             trace.it = F)
                
                model <- glmnet ( x = as.matrix(mets_train_sample[,-1]),  y= mets_train_sample[,1], 
                                  lambda = glmnet_Modelfit$lambda.1se, alpha = 1,
                                  penalty.factor = penalty.factor,
                                  nfolds = 10)
                
                test_value <- unname(predict(model,
                                             as.matrix(mets_test[,-1])))[,1]
                
                ensemble_model_list[[i]] <- model
        }
        
        # feature selection stability
        # extract coefficients from all 500 models
        coef_matrix <- do.call(cbind, lapply(ensemble_model_list, function(x) {
                as.vector(coef(x))
        }))
        rownames(coef_matrix) <- rownames(coef(ensemble_model_list[[1]]))

        selection_freq <- rowMeans(coef_matrix != 0)
        selection_freq <- sort(selection_freq, decreasing = TRUE)
        mean_coef <- rowMeans(coef_matrix)[names(selection_freq)]
        write.table(cbind(selection_freq, mean_coef), paste("LASSO_selection_stability_", pollutant, ".txt", sep = ""), row.names = T, quote = F, sep = "\t")
        
        # library(ggplot2)
        # covariate_names <- colnames(mets_train)[2:16]
        # freq_df <- data.frame(
        #         variable  = names(selection_freq),
        #         frequency = selection_freq
        # ) |> subset(!variable %in% c("(Intercept)", covariate_names))
        # # remove covariates and cell proportions
        # 
        # p <- ggplot(freq_df, aes(x = reorder(variable, frequency), y = frequency)) +
        #         geom_bar(stat = "identity", fill = "steelblue", width = 0.7) +
        #         geom_hline(yintercept = 0.5, colour = "red", linetype = "dashed", linewidth = 0.8) +
        #         annotate("text", x = 1, y = 0.52, label = "50% threshold", colour = "red", size = 3, hjust = 0) +
        #         coord_flip() +
        #         labs(x = "CpG Sites", y = "Selection Frequency",
        #              title = "LASSO Feature Selection Stability across 500 Bootstrap Models") +
        #         theme_classic() +
        #         theme(axis.text.y = element_text(size = 7))
        # 
        # plot(p)
        
        model_predictions <-
                do.call(cbind.data.frame, lapply(ensemble_model_list, function(x) {
                        unname(predict(x, as.matrix(mets_test[,-1])))[, 1]
                }))
        methyExposure <- data.frame('exposure' = mets_test[,1])
        row.names(methyExposure) <- row.names(mets_test)
        methyExposure["Lasso_bagged_ensemble"] <- rowMeans(model_predictions)
        return(methyExposure)
}


library(ggplot2)
results_df <- Lasso_Pollutant("NO2")

# correlation and RMSE for annotation
r    <- cor(results_df$exposure, results_df$Lasso_bagged_ensemble)
rmse <- sqrt(mean((results_df$exposure - results_df$Lasso_bagged_ensemble)^2))

ggplot(results_df, aes(x = exposure, y = Lasso_bagged_ensemble)) +
        geom_point(alpha = 0.6, size = 1.5) +
        geom_abline(slope = 1, intercept = 0, colour = "red", linetype = "dashed") +  # perfect prediction line
        geom_smooth(method = "lm", colour = "blue", se = TRUE) +                      # actual fit line
        annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -1,
                 label = sprintf("r = %.3f\nRMSE = %.3f", r, rmse)) +
        labs(x = "Actual Exposure", y = "Predicted Exposure",
             title = "DNAm Predicted vs Actual Air Pollution Exposure") +
        theme_classic()

# funnel plot for selection frequency and mean coefficients
library(ggplot2)
df <- read.delim("LASSO_selection_stability_NO2.txt")
df <- df[-c(1:16),]
df$CpG <- rownames(df)
df$highlight <- df$selection_freq > 0.8

library(ggplot2)
library(ggrepel)

ggplot(df, aes(x = selection_freq, y = mean_coef)) +
        geom_point(aes(color = highlight), size = 1.8) +
        scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "#1F4E79")) +
        
        # framed, non-overlapping labels for highlighted CpGs
        geom_label_repel(
                data = subset(df, highlight),
                aes(label = CpG),
                color = "#1F4E79",
                fill = "white",
                label.size = 0.5,   # square frame thickness
                label.padding = unit(0.15, "lines"),
                size = 3,
                max.overlaps = Inf,
                box.padding = 0.5,
                point.padding = 0.3,
                segment.color = "grey50",
                segment.size = 0.3,
                min.segment.length = 0,
                force = 1,
                max.time = 2,
                seed = 42
        ) +
        
        theme_bw() +
        labs(
                x = "Bootstrap model selection proportion",
                y = "Mean standardized coefficient",
                color = "Selected > 0.8",
                title = "Bootstrap Aggreate LASSO model"
        ) +
        theme(
                legend.position = "none",
                panel.grid.minor = element_blank()
        )


## plot performance for training model
Lasso_Pollutant_test_performance <- function(pollutant){
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
        set.seed(1234) # make sure train-test random split the same everytime as EWAS
        test_idx  <- sample(nrow(model_mtx), 0.2 * nrow(model_mtx))
        train_idx <- setdiff(seq_len(nrow(model_mtx)), test_idx)
        
        mets_train <- model_mtx[train_idx,]
        mets_train <- makeX(mets_train, na.impute = TRUE) ## impute missing values
        
        mets_test <- model_mtx[test_idx,]
        mets_test <- makeX(mets_test, na.impute = TRUE)
        ## 3. model training
        library(glmnet)
        penalty.factor <- c(
                rep(0, 15),  # covariates — never shrunk out // factors are treated as individual binary features
                rep(1, length(cpgs))  # CpGs — LASSO penalized
        )
        
        # Modify your loop to also store test predictions and performance metrics
        Nbootstrap <- 500
        ensemble_model_list <- list()
        perf_list <- vector("list", Nbootstrap)
        
        y_test <- mets_test[,1]
        
        for (i in 1:Nbootstrap){
                print(paste("we are now running bootstrap", i, sep = " "))
                select_id <- sample.int(nrow(mets_train), replace = TRUE)
                mets_train_sample <- mets_train[select_id,]
                
                glmnet_Modelfit <- cv.glmnet(x = as.matrix(mets_train_sample[,-1]), y = mets_train_sample[,1],
                                             alpha = 1, nfolds = 10,
                                             penalty.factor = penalty.factor,
                                             trace.it = FALSE)
                
                model <- glmnet(x = as.matrix(mets_train_sample[,-1]), y = mets_train_sample[,1],
                                lambda = glmnet_Modelfit$lambda.1se, alpha = 1,
                                penalty.factor = penalty.factor)
                
                train_pred <- unname(predict(model, as.matrix(mets_train_sample[,-1])))[,1]
                y_train_sample <- mets_train_sample[,1]
                resid <- y_train_sample - train_pred
                
                # --- performance metrics for this bootstrap model ---
                corr <- cor(y_train_sample, train_pred)
                r2   <- 1 - sum(resid^2) / sum((y_train_sample - mean(y_train_sample))^2)
                rmse  <- sqrt(mean(resid^2))
                perf_list[[i]] <- data.frame(bootstrap = i, R2 = r2, Corr = corr, RMSE = rmse)
                
                ensemble_model_list[[i]] <- model
        }
        
        perf_df <- do.call(rbind, perf_list)
        return(perf_df)
}

library(ggplot2)
library(tidyr)
library(glmnet)
perf_df <- Lasso_Pollutant_test_performance("NO2")

perf_long <- pivot_longer(perf_df, cols = c(R2, Corr, RMSE),
                          names_to = "metric", values_to = "value")

ggplot(perf_long, aes(x = metric, y = value)) +
        geom_boxplot(fill = "steelblue", alpha = 0.6, outlier.alpha = 0.3) +
        geom_jitter(width = 0.15, alpha = 0.2, size = 0.8) +
        facet_wrap(~metric, scales = "free") +
        theme_bw() +
        labs(title = paste0("Bootstrap model performance (n = ", 500, ")"),
             x = NULL, y = "Value")

all_test_preds <- sapply(ensemble_model_list, function(m) {
        unname(predict(m, as.matrix(mets_test[,-1])))[,1]
})  # matrix: nrow(mets_test) x Nbootstrap

ensemble_pred <- rowMeans(all_test_preds)
ensemble_resid <- y_test - ensemble_pred
ensemble_r2 <- 1 - sum(ensemble_resid^2) / sum((y_test - mean(y_test))^2)
ensemble_rmse <- sqrt(mean(ensemble_resid^2))
