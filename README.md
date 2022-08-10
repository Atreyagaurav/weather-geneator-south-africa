# Principal component analysis on NCEP data

## Objective
We're trying to use PCA on NCEp data for the weather generator.

## Download 
Data is downloaded using `NCEP` library, you can download it using:

```R
install.packages("RNCEP")
```

The file `download.r` has the code to download, calculate yearly aggregate and save tabular data. It's made to be run interactively.

## PCA analysis

PCA analysis is done to reduce the dimensionality of the data. It took me a while to understand the dimensions of this data as initially I thought it was like point data with dimensions lat, lon and time. Hence no reason to reduce the dimensionality. 

Now I'm come to the conclusion that, the data isn't the point but rather a state, which includes all the gridded data at single time frame. Which means at a single time we have `lat × lon` number of points, and it's a matrix data, hence we have  `lat × lon` dimensions for each data.

## Kmeans clustering
After PCA analysis was done, then those ordinates were used for k means clustering to get 8 points. After doing the reverse transformation from them, we got the 8 clusters in original raster format. cluster rasters are in rasters directory.

## Markov Chain
After we had clusters, and categorization of each points in the time series. Then each pattern was simply counted and then converted to observed probabilities.

The results are sent in slack. I'll save the image files and put it here later.

Results are in the file: [file:./csvs/tables.ods]
