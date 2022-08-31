#!/usr/bin/env python3
import netCDF4
import numpy as np
import os
# import sys

# hgt_nc_file = sys.argv[1]
hgt_nc_file = "ncep-reanalysis-2/4xdaily-500hpa.nc"
nc = netCDF4.Dataset(hgt_nc_file, mode='r')

lat = nc.variables['lat'][:]
lon = nc.variables['lon'][:]
time_var = nc.variables['time']
dtime = netCDF4.num2date(time_var[:], time_var.units)

# originally 5 for 500hpa, but already filtered so there is only one level
hgt = nc.variables['hgt'][:, 0, :, :]

hgt_daily = np.zeros((hgt.shape[0]//4, *hgt.shape[1:]))

for i in range(hgt.shape[0]//4):
    hgt_daily[i] = sum(hgt[4*i:4*(i+1)])/4

del(hgt)

days = sorted({d.strftime("%Y-%m-%d") for d in dtime})

hgt_ltm_days = dict()
for i, d in enumerate(days):
    day = d.split("-", maxsplit=1)[1]
    if day in hgt_ltm_days:
        hgt_ltm_days[day].append(i)
    else:
        hgt_ltm_days[day] = [i]

hgt_anomaly = np.zeros(hgt_daily.shape)

for day, ind in hgt_ltm_days.items():
    hgt_ltm = sum(hgt_daily[ind, :, :])/len(ind)
    for i in ind:
        hgt_anomaly[i] = hgt_daily[i] - hgt_ltm


fname, _ = os.path.splitext(hgt_nc_file)
outfile = f'{fname}-anomaly.csv'

with open(outfile, "w") as writer:
    writer.write('time,lat,lon,delta_hgt\n')
    for lat_i, lat_v in enumerate(lat):
        for lon_i, lon_v in enumerate(lon):
            for time_i, time_v in enumerate(days):
                line = f'{time_v},{lat_v},{lon_v}'
                writer.write(f'{line},{hgt_anomaly[time_i,lat_i,lon_i]}\n')
