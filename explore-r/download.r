library(RNCEP)
library(fields)
library(readr)

wx.extent <- NCEP.gather(variable = "hgt", level=500,
                         months.minmax =c(1,12), years.minmax = c(1979,2021),
                         lat.southnorth = c(-50, -15), lon.westeast = c(330, 60),
                         reanalysis2 = TRUE, return.units = TRUE)

## The data array may be saved directly as an R object ##
save(wx.extent, file='wx_extent.Rdata')
## And then later recalled ##
load(file='wx_extent.Rdata')


## Now calculate the average temperature perday ##
wx.ag <- NCEP.aggregate(wx.data=wx.extent, YEARS=TRUE, MONTHS=TRUE,
                        DAYS=TRUE, HOURS=FALSE, fxn='mean')


## The data array may be saved directly as an R object ##
save(wx.ag, file='wx_agg.Rdata')
## And then later recalled ##
load(file='wx_agg.Rdata')



wx.ag.min = floor(min(wx.ag)/100)* 100
wx.ag.max = ceiling(max(wx.ag)/100)* 100
zlim = c(wx.ag.min, wx.ag.max)

for (i in 1:50) {
png(sprintf("/tmp/plt-%03d.png", i))
## NCEP.vis.area(wx.ag, layer = i, show.pts = FALSE, draw.contours = FALSE,
##               cols = terrain.colors(64, zlim=c(4e3,6e3)),
##               map.args = list(regions="south africa"),
##               image.plot.args = list(legend.args=list(zlim=c(5e3,6e3))))
image.plot(wx.ag[,,i], zlim=zlim)
dev.off()}


NCEP.vis.area(wx.ag, layer = 2, show.pts = FALSE, draw.contours = FALSE,
              cols = terrain.colors(64),
              image.plot.args = list(zlim=zlim),
              map.args = list(regions="south africa"))

image.plot(wx.ag[,,200], zlim=zlim)

image.plot.args = list(wx.ag[,,200], zlim=zlim)
do.call("image.plot", image.plot.args)


## Convert annual aggregate to tabular format.
df <- NCEP.array2df(wx.ag)


## function to parse the date
parse_date_fmt <- function (datetime_str){
    date <- readr::parse_date(sprintf("%.10s", datetime_str), "%Y_%m_%d")
    return (date)
}

start <- parse_date_fmt("1979_01_01_XX")
days_after_1979 <- function(datetime_str){
    days <- as.numeric(parse_date_fmt(datetime_str) - start)
}

df$time_days <- days_after_1979(df$datetime)
names(df)[5] <- "geopotential_ht"
write.csv(df, "ncep-data.csv")
