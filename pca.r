library(RNCEP)

library(tidyverse)  # data manipulation and visualization
library(gridExtra)  # plot arrangement
library(ggplot2)


## wx.ag variable
load(file='wx_agg.Rdata')

dim(wx.ag)

tabular <- read.csv("ncep-data.csv")

head(tabular)
variance <- apply(tabular, 2, var)

dropped_df <-subset(tabular, select = -c(datetime, X))
head(dropped_df)
names(dropped_df)[3] <- "geopotential_ht"
scaled_df <- apply(dropped_df, 2, scale)
head(scaled_df)
tail(scaled_df)

## ggplot(as.data.frame(scaled_df)) + geom_point(mapping = aes(x=latitude, y=longitude))

## principle component analysis
# Calculate eigenvalues & eigenvectors
scaled.cov <- cov(scaled_df)
scaled.eigen <- eigen(scaled.cov)
str(scaled.eigen)

# Extract the loadings
(phi <- scaled.eigen$vectors[,1:2])
phi

row.names(phi) <- names(dropped_df)
colnames(phi) <- c("PC1", "PC2")
phi

PC1 <- as.matrix(scaled_df) %*% phi[,1]
PC2 <- as.matrix(scaled_df) %*% phi[,2]

# Create data frame with Principal Components scores
PC <- data.frame(Day = row.names(tabular), PC1, PC2)
head(PC)

## Plot Principal Components for each State
ggplot(PC, aes(PC1, PC2)) + 
    modelr::geom_ref_line(h = 0) +
    modelr::geom_ref_line(v = 0) +
    geom_point() +
    ## geom_text(aes(label = Day), size = 3) +
    xlab("First Principal Component") + 
    ylab("Second Principal Component") + 
    ggtitle("First Two Principal Components of USArrests Data")

