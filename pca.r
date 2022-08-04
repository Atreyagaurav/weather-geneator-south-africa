library(RNCEP)

## library(tidyverse)  # data manipulation and visualization
library(gridExtra)  # plot arrangement
library(ggplot2)
library(plotly)

## wx.ag variable
## load(file='wx_agg.Rdata')

## dim(wx.ag)

tabular <- read.csv("ncep-data.csv")

head(tabular)

tabular$latlon <- sprintf("%+.2f,%+.2f", tabular$longitude, tabular$latitude)


dropped_df <-subset(tabular, select = c(time_days, latlon, geopotential_ht))
head(dropped_df)

## 5717.384
mean(dropped_df$geopotential_ht)
MEAN_GHT = 5717.384

table_wide <- spread(apply(dropped_df, 2, scale),
                     key=latlon,
                     value=geopotential_ht)
rownames(table_wide) <- str(table_wide$time_days)
table_wide <- subset(table_wide, select = -c(time_days))
table_wide <- table_wide - MEAN_GHT
head(table_wide)

variance <- apply(table_wide, 2, var)

## scaled_df <- apply(table_wide, 2, scale)
## head(scaled_df)
## tail(scaled_df)

## ggplot(as.data.frame(scaled_df)) + geom_point(mapping = aes(x=latitude, y=longitude))

## principle component analysis
# Calculate eigenvalues & eigenvectors
table_wide.cov <- cov(table_wide)
table_wide.eigen <- eigen(table_wide.cov)
str(table_wide.eigen)

# Extract the loadings
(phi <- table_wide.eigen$vectors[,1:3])
phi

row.names(phi) <- names(dropped_df)
colnames(phi) <- c("PC1", "PC2", "PC3")
phi

ALL_PC <- as.matrix(table_wide) %*% phi
PC1 <- ALL_PC[,1]
PC2 <- ALL_PC[,2]
PC3 <- ALL_PC[,3]

# Create data frame with Principal Components scores
PC <- data.frame(Day = row.names(table_wide), PC1, PC2, PC3)
head(PC)

## plots

options(browser = 'firefox')
plot_ly(x=PC$PC1, y=PC$PC2, z=PC$PC3, type="scatter3d", mode="markers", color=PC1)

## ## Plot Principal Components for each State
## ggplot(PC, aes(PC1, PC2)) + 
##     ## modelr::geom_ref_line(h = 0) +
##     ## modelr::geom_ref_line(v = 0) +
##     geom_point() +
##     ## geom_text(aes(label = Day), size = 3) +
##     xlab("First Principal Component") + 
##     ylab("Second Principal Component") + 
##     ggtitle("First Two Principal Components of NCEP data")


## get back original data (approx)
table_appx = t(phi %*% t(ALL_PC))
table_appx[1:10, 1]
table_wide[1:10, 1]
