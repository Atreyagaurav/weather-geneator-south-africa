library(raster)
library(tidyverse)
## wx.ag variable
## load(file='wx_agg.Rdata')

## dim(wx.ag)

tabular <- read.csv("ncep-data.csv")

head(tabular)

tabular$lonlat <- sprintf("%+.2f,%+.2f", tabular$longitude, tabular$latitude)


dropped_df <-subset(tabular, select = c(time_days, lonlat, geopotential_ht))
head(dropped_df)


table_wide <- spread(dropped_df,
                     key=lonlat,
                     value=geopotential_ht)
rownames(table_wide) <- str(table_wide$time_days)
table_wide <- subset(table_wide, select = -c(time_days))
table_wide <- table_wide
head(table_wide, 1)

pca <- prcomp(table_wide, rank = 50, scale. = T)
summary(pca)

table_appx <- t(t(pca$x %*% t(pca$rotation)) * pca$scale + pca$center)

table_appx[1:4,1:5]
table_wide[1:4,1:5]

table_appx_long = pivot_longer(as.data.frame(table_appx), cols = -c(), names_to = "lonlat", values_to = "geopotential_ht")
head(table_appx_long)
## org_raster <- rasterFromXYZ()

km <- kmeans(pca$x, centers = 8, nstart = 5)

kmeans_centers = km$centers
kmeans_revert <- t(t(kmeans_centers %*% t(pca$rotation)) * pca$scale + pca$center)


plot_raster <- function(df) {
    df <- as.data.frame(df)
    df$lat <- 0
    df$lon <- 0

    for (i in 1:length(df$lonlat)){
        ll = df$lonlat[i]
        ll_split <- str_split_fixed(ll, ",", n=2)
        lon = as.numeric(ll_split[1])
        lat = as.numeric(ll_split[2])
        df$lat[i] <- lat
        df$lon[i] <- lon
    }

    xyz <- subset(df, select = c("lon", "lat", "geopotential_ht"))

    dfr <- rasterFromXYZ(xyz = xyz)
    plot(dfr)
}

N = 4
means_N <- data.frame(lonlat = names(kmeans_revert[N,]), geopotential_ht=as.numeric(kmeans_revert[N,]))
plot_raster(means_N)

day=20
means_N <- data.frame(lonlat = names(table_appx[day,]), geopotential_ht=as.numeric(table_appx[day,]))
plot_raster(means_N)

day=20
means_N <- data.frame(lonlat = names(table_wide[day,]), geopotential_ht=as.numeric(table_wide[day,]))
plot_raster(means_N)

