import numpy as np
import netCDF4 as nc4
import glob
import gzip
import os
import shutil
import tempfile
import pandas as pd
import csv

#===================================================================
# Loading and preparing data
#===================================================================

# Empty lists for storing variables
time_list = []; day_list = []; granule_list = []; year_list = []; surface_list = []
lat_co_list = []; lon_co_list = []; co_pix_i_list = []; co_pix_f_list = []
co_length_list = []; core_cog_list = []; core_rcog_list = []
core_pix_i_list = []; core_pix_f_list = []; num_core_list = []; core_length_list = [] 

years = ['2009', '2010'] # Edit this to add your own years

file_path = '/../doppler/data8/pilewskie/convective_object_database/Version_2/' # Edit this text to add your own file path

# Loading in CO database
for m, year in enumerate(years):

    files = sorted(glob.glob(file_path + year + '/' + '*.nc'))

    if len(files) == 0:
        files = sorted(glob.glob(file_path + year + '/' + '*.nc.gz'))

    for f in files:
        
        print(f)

        # Loading in data
        if f.endswith(".gz"):
            infile = gzip.open(f, 'rb')
            tmp = tempfile.NamedTemporaryFile(dir='/dev/shm/')
            shutil.copyfileobj(infile, tmp)
            data = nc4.Dataset(tmp.name, 'r', format = 'NETCDF4')
            #os.unlink(tmp.name)
            infile.close()
            tmp.close()
        else:
            data = nc4.Dataset(f, 'r', format = 'NETCDF4')
   
        #============
        # Grabbing data stored in NetCDFs
        #============

        # General information
        granule = data.variables['CloudSat Granule'][:]
        time = data.variables['Time of Day'][:]; 
        day = data.variables['Julian Day'][:]
        surface_flag = data.variables['Surface flag'][:] # Determining if it is over land or ocean
        lat = data.variables['Latitudes of CO'][:]
        lon = data.variables['Longitudes of CO'][:]
        freeze = data.variables['Mean Freezing Level'][:]

        # Convective object variables
        co_length = data.variables['CO Length'][:]
        co_pix_i = data.variables['Initial CO Profile Index'][:]
        co_pix_f = data.variables['Final CO Profile Index'][:]

        # Convective core variables
        core_pix_i = data.variables['Initial Core Index'][:] 
        core_pix_f = data.variables['Final Core Index'][:]
        num_core = data.variables['No of Convective Cores'][:]
        core_length = data.variables['Length of Convective Cores'][:] # Multiply by 1.1 km to get in units of km 
        CoG = data.variables['Z-weighted Mean Core CoG'][:]

        # Make data workable by reshaping and getting rid of masks
        for i in range(len(co_length)): # indexing over convective objects in each CloudSat granule

            if (CoG[i] is np.ma.masked) == True:
                continue
            elif (co_length[i] is np.ma.masked) == True:
                continue
            elif (num_core[i] is np.ma.masked) == True:
                continue
            else:
    
                # Initial latitude and longitude of the CO
                lat_co_list.append(lat[i,0,:co_length[i]])
                lon_co_list.append(lon[i,0,:co_length[i]])              

                # Temporal information
                time_list.append(time[i].astype(int))
                day_list.append(day.astype(int))
                granule_list.append(granule[0].astype(int))
                year_list.append(int(year))               

                # Convective object variables
                co_pix_i_list.append(co_pix_i[i])
                co_pix_f_list.append(co_pix_f[i]) 
                surface_list.append(surface_flag[i].astype(int))
                co_length_list.append(co_length[i]*1.1)

                # Convective core variables
                core_pix_i_list.append([x for x in core_pix_i[i] if np.ma.is_masked(x) == False])  
                core_pix_f_list.append([x for x in core_pix_f[i] if np.ma.is_masked(x) == False]) 
                core_cog_list.append([x for x in CoG[i] if np.ma.is_masked(x) == False]) 
                num_core_list.append(num_core[i].astype(int))
                core_length_list.append([x*1.1 for x in core_length[i] if np.ma.is_masked(x) == False])           
            
        data.close()
            
        del granule, time, day, surface_flag, lat, lon, freeze
        del co_length, co_pix_i, co_pix_f, core_pix_i, core_pix_f, num_core, core_length, CoG

#============
# Sorting values by longitude to ease with plotting
#============

# Index sorted by longitude to use for reordering the variables
idx_sorted = [i[0] for i in sorted(enumerate(lon_co), key=lambda x:x[1])]

# General information   
lon_sorted = [lon_co[i] for i in idx_sorted]
lat_sorted = [lat_co[i] for i in idx_sorted]
time_sorted = [time_list[i] for i in idx_sorted]
day_sorted = [day_list[i] for i in idx_sorted]
granule_sorted = [granule_list[i] for i in idx_sorted]
year_sorted = [year_list[i] for i in idx_sorted]
surface_sorted = [surface_list[i] for i in idx_sorted]

# Convective object variables
co_pix_i_sorted = [co_pix_i_list[i] for i in idx_sorted]
co_pix_f_sorted = [co_pix_f_list[i] for i in idx_sorted]
co_length_sorted = [co_length_list[i] for i in idx_sorted]

# Convective core variables
core_pix_i_sorted = [core_pix_i_list[i] for i in idx_sorted]
core_pix_f_sorted = [core_pix_f_list[i] for i in idx_sorted]
core_cog_sorted = [core_cog_list[i] for i in idx_sorted]
num_core_sorted = [num_core_list[i] for i in idx_sorted]
core_size_sorted = [core_size_list[i] for i in idx_sorted]

# Store mean coordinates and unmasked/reshaped data as dictionaries to call for plotting
keys = ['Mean longitude', 'Mean latitude', 'Julian day', 'Day or Night', 'Granule', 'Year', 'Surface flag', \
        'Initial CO pixel', 'Final CO pixel', 'CO length', 'Initial core index', 'Final core index', \
        'Mean Core CoG', 'Number of cores', 'Core sizes']

data = [lon_sorted, lat_sorted, day_sorted, time_sorted, granule_sorted, year_sorted, surface_sorted, \
        co_pix_i_sorted, co_pix_f_sorted, co_length_sorted, core_pix_i_sorted, core_pix_f_sorted, \
        core_cog_sorted, num_core_sorted, core_size_sorted]

co_dict = {key: value for key,value in zip(keys, data)}

df = pd.DataFrame(co_dict)
df.to_csv('/../doppler/data8/pilewskie/convective_object_database/postprocessed_data/example_COs.csv', index=False) # Enter your own file path & file name here


