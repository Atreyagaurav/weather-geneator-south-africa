library(raster)
library(tidyverse)
library(stringr)
## wx.ag variable
## load(file='wx_agg.Rdata')

## dim(wx.ag)

tabular <- read.csv("csvs/daily-500hpa-anomaly.csv")

head(tabular)
tabular$lon <- tabular$lon-360*(tabular$lon > 180)
tabular$lonlat <- sprintf("%+.2f,%+.2f", tabular$lon, tabular$lat)


dropped_df <-subset(tabular, select = c(time, lonlat, delta_hgt))
head(dropped_df)


table_wide <- spread(dropped_df,
                     key=lonlat,
                     value=delta_hgt)
rownames(table_wide) <- str(table_wide$time)
table_wide <- subset(table_wide, select = -c(time))
table_wide <- table_wide
head(table_wide, 1)

## save this for scott's code
saveRDS(table_wide, 'data/hgt-SA-anomaly.rds')

pca <- prcomp(table_wide, rank = 10, scale. = TRUE)

summary(pca)

table_appx <- t(t(pca$x %*% t(pca$rotation)) * pca$scale + pca$center)

table_appx[1:4,1:5]
table_wide[1:4,1:5]

table_appx_long <- pivot_longer(as.data.frame(table_appx), cols = -c(), names_to = "lonlat", values_to = "delta_hht")

head(table_appx_long)
## org_raster <- rasterFromXYZ()

km <- kmeans(pca$x, centers = 6, nstart = 10)

kmeans_centers <- km$centers
kmeans_revert <- t(t(kmeans_centers %*% t(pca$rotation)) * pca$scale + pca$center)

regimes <- data.frame(date = sort(unique(tabular$time)), clusters = km$cluster)

regimes$date <- str_trunc(str_replace_all(regimes$date, "_", "-"), 10, ellipsis = "")
head(regimes)
write.csv(regimes, "csvs/simple-kmeans-6.csv")


get_raster <- function(df) {
    df <- as.data.frame(df)
    df$lat <- 0
    df$lon <- 0

    for (i in 1:length(df$lonlat)){
        ll <- df$lonlat[i]
        ll_split <- str_split_fixed(ll, ",", n=2)
        lon <- as.numeric(ll_split[1])
        lat <- as.numeric(ll_split[2])
        df$lat[i] <- lat
        df$lon[i] <- lon
    }

    xyz <- subset(df, select = c("lon", "lat", "geopotential_ht"))

    dfr <- rasterFromXYZ(xyz = xyz)
    return(dfr)
}


## for (N in 1:50){
##     pca_X <- data.frame(lonlat = names(pca$x[N,]), geopotential_ht=as.numeric(pca$x[N,]))
##     dfr1 <- get_raster(pca_X)
##     writeRaster(dfr1, sprintf("./rasters/pca/component-%02d.tif", N))
##     ## plot(dfr1)
## }


for (N in 1:13){
    means_N <- data.frame(lonlat = names(kmeans_revert[N,]), geopotential_ht=as.numeric(kmeans_revert[N,]))
    dfr1 <- get_raster(means_N)
    writeRaster(dfr1, sprintf("rasters/kmeans.old/cluster-%d.tif", N))
    ## plot(dfr1)
}

for (day in 1:20) {
    original <- data.frame(lonlat = names(table_wide[day,]), geopotential_ht=as.numeric(table_wide[day,]))
    dfr1 <- get_raster(original)
    writeRaster(dfr1, sprintf("./rasters/original/day-%02d.tif", day))
}

for (day in 1:20) {
    recreated <- data.frame(lonlat = names(table_appx[day,]), geopotential_ht=as.numeric(table_appx[day,]))
    dfr1 <- get_raster(recreated)
    writeRaster(dfr1, sprintf("./rasters/recreated/day-%02d.tif", day))
}
