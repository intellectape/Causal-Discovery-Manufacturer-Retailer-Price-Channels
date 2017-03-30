# This project is created by:
# Name: Aditya Bhardwaj
# Unity ID: abhardw2
# Project Title: Project 7 - Causal Discovery between Manufacturer-Retailer Price Channels

# Load the libraries 

library(vars)
library(urca)
library(pcalg)
library(stats)
library(lmtest)
library(forecast)
library(tseries)

# Read the input data 

csv_file_path <- "BI/Projects/Project7/Input Data/data.csv"
tetrd_file_path <- "BI/Projects/Project7/Input Data/tetrd.csv"

data_causal <- read.csv(csv_file_path, sep = ",", header = TRUE)
ts_causal <- ts(data_causal, frequency = 52)


# Build a VAR model 

qty <- data_causal[, 1]
rprice <- data_causal[, 2]
mprice <- data_causal[, 3]

# Plot of the data

plot.ts(qty)
plot.ts(rprice)
plot.ts(mprice)

VARselect(data_causal, lag.max = 10, type="const")
var_model <- VAR(data_causal, p=1, type="const")

# Extract the residuals from the VAR model 

qty_res <- var_model$varresult$Move$residuals
rprice_res <- var_model$varresult$RPRICE$residuals
mprice_res <- var_model$varresult$MPRICE$residuals

# Check for stationarity using the Augmented Dickey-Fuller test 

qty_lag <- ndiffs(qty_res)
mprice_lag <- ndiffs(mprice_res)
rprice_lag <- ndiffs(rprice_res)

# Augmented Dickey-Fuller test
# used urca library as it gives more precise results of the test than the adf.test
qty_adf <- ur.df(qty_res, lags = qty_lag, type="none")
mprice_adf <- ur.df(mprice_res, lags = mprice_lag, type="none")
rprice_adf <- ur.df(rprice_res, lags = rprice_lag, type="none")

summary(qty_adf)
summary(mprice_adf)
summary(rprice_adf)


# Check whether the variables follow a Gaussian distribution  

ks.test(qty_res,"pnorm", mean = mean(qty_res), sd = sd(qty_res))
ks.test(mprice_res,"pnorm", mean = mean(mprice_res), sd = sd(mprice_res))
ks.test(rprice_res,"pnorm", mean = mean(rprice_res), sd = sd(rprice_res))

# Write the residuals to a csv file to build causal graphs using Tetrad software
residual_data <- cbind(qty_res, rprice_res, mprice_res)
write.csv(residual_data, file = tetrd_file_path, row.names = FALSE)

# PC Algorithm
suffStat =  list(C=cor(residual_data), n = nrow(residual_data))
pc_algo_res <- pc(suffStat, indepTest = gaussCItest, alpha = 0.1, labels = colnames(residual_data), skel.method = "original")
plot(pc_algo_res, main="PC Algorithm Output")

# LinGam Algorithm
lingam_algo_res <- lingam(residual_data, verbose = TRUE)
show(lingam_algo_res)
plotling <- as(lingam_algo_res, "amat")
