import numpy as np
import scipy.io as sio
import calendar
from getconstants import getconstants
from timeit import default_timer as timer
import os
import netCDF4 as nc
import xarray as xr
from affine import Affine
from rasterio import features
import geopandas as gpd
import matplotlib.pyplot as plt
import sys
import time
import pandas as pd
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
from matplotlib.colors import LinearSegmentedColormap
import math



start1 = timer()

def conv(array):
  arrayf=np.zeros(np.shape(array),dtype=float)
  arrayf[:,:]=np.nan
  arrayf[:,0:180]=array[:,180:360]
  arrayf[:,180:360]=array[:,0:180]
  return arrayf


def transform_from_latlon(lat,lon):
  lat=np.asarray(lat)
  lon=np.asarray(lon)
  trans=Affine.translation(lon[0],lat[0])
  scale= Affine.scale(lon[1]-lon[0],lat[1]-lat[0])
  return trans*scale

def rasterize(shapes, coords, latitude='latitude', longitude='longitude', fill=np.nan, **kwargs):
  transform = transform_from_latlon(coords[latitude], coords[longitude])
  out_shape= (len(coords[latitude]),len(coords[longitude]))
  raster = features.rasterize(shapes, out_shape=out_shape, fill=fill, transform=transform, dtype=float, **kwargs)
  spatial_coords={latitude:coords[latitude], longitude: coords[longitude]}
  return xr.DataArray(raster, coords=spatial_coords, dims=(latitude, longitude))



def acumulator2(array,lat,lon,cats,Region):
  flatarray=array.flatten()
  lat1D=np.tile(lat,len(lon))
  lon1D=np.tile(lon,len(lat))
  #print(len(flatarray))
  #print(len(lat1D))
  #print(len(lon1D))
  cats=np.array(cats)
  df = pd.DataFrame({'latitude':lat1D,'longitude':lon1D, 'Var':flatarray})
  df = df.sort_values('Var')
  varsort=np.array(df['Var'])
  accum=np.cumsum(varsort)
  print(np.nanmin(accum))
  print(np.nanmax(accum))
  accum=accum*100/np.nanmax(accum)
  print(np.nanmax(accum))
  id_cats=np.zeros(len(accum))
  for i in range(0,len(accum)):
    if accum[i]<=0:
      id_cats[i]=0
    for j in range(1,len(cats)):
      if math.trunc(accum[i])<=cats[j] and accum[i]>=cats[j-1]:
        id_cats[i]=j

  df['Cat']=id_cats  
  df['Accumulator']=accum
  df = df.sort_index()
  categories=np.array(df['Cat'])
  accumulator=np.array(df['Accumulator'])
  variable=np.array(df['Var'])
  categories=categories.reshape(np.shape(array))
  accumulator=accumulator.reshape(np.shape(array))
  variable=variable.reshape(np.shape(array)) 
  print(categories)
  #print(categories)
  #categories[Region == 0]=100
  #categories[Region[0:180,:] == 0]=0        
  return variable, accumulator, categories


def save_nc(new_filename,lats,lons,rh):
  from netCDF4 import Dataset
  # write netCDF file
  #------------------#
  # open a netCDF file to write
  ncout1 = Dataset(new_filename, 'w', format='NETCDF4')
  # define axis size
  ncout1.createDimension('lat', len(lats))
  ncout1.createDimension('lon', len(lons))
  # create latitude axis
  lat = ncout1.createVariable('lat', 'f4', ('lat'))
  lat.standard_name = 'latitude'
  lat.long_name = 'latitude'
  lat.units = 'degrees_north'
  lat.axis = 'Y'#
  # create longitude axis
  lon = ncout1.createVariable('lon', 'f4', ('lon'))
  lon.standard_name = 'longitude'
  lon.long_name = 'longitude'
  lon.units = 'degrees_east'
  lon.axis = 'X'#


  # create variable array
  vout1 = ncout1.createVariable('recycling ratio', 'f4', ('lat', 'lon'))
  vout1.long_name = 'Total accumulated precipitation'
  vout1.units = '--'#

  lon[:] = lons[:]
  lat[:] = lats[:]
  vout1[:] = rh[:]
  ncout1.close()


#list_basins=[354,355,356,357,358,359,360,361,362,363,364,366,367,368,369,370,371,372,373,374,375,377,378,379,380,382,383,384,385,386,387,388,391,392,393,395,396,397,398,399,400,403,404]

list_basins=range(1,406,1)

##CARGA DE ARCHIVOS
#Cargar shapefile
#watershed = '/media/user/Data/CCI-Scripts/Cauca-Magdalena.shp'
watershed = '../GRDC_405_basins_from_mouth_Pr.shp'
wshed=gpd.read_file(watershed)
Rr=np.zeros(len(wshed['BASIN_ID']))
data={'ID':wshed['BASIN_ID'],'RecyclinR':Rr}
df_data=pd.DataFrame(data=data)
#polygon = wshed['geometry']



for basin in list_basins:
  polygon = wshed[wshed.BASIN_ID==(basin)].geometry
  print(polygon)



  # obtain the constants
  invariant_data = './Invariant.nc'

  latnrs = np.arange(0,180)
  lonnrs = np.arange(0,360)
  # load the latitude and longitude from the invariants file
  latitude1 = nc.Dataset(invariant_data, mode = 'r').variables['latitude'][latnrs] # [degrees north]
  longitude1 = nc.Dataset(invariant_data, mode = 'r').variables['longitude'][lonnrs] # [degrees east]


  # BEGIN OF INPUT 2 (FILL THIS IN)
  ##cargar shapefile


  ##Rasterizar
  #obtener raster de referencia
  inv=nc.Dataset(invariant_data)
  ar_lsm=(inv['lsm'][0,:,:])
  latitude=np.asarray(inv['latitude'])
  longitude=np.asarray(inv['longitude'])

  xr_lsm=xr.DataArray(ar_lsm, coords=[latitude, longitude], dims=['y','x'])
  #convertir shape
  pshed_mask= rasterize(polygon, xr_lsm.coords, longitude='x', latitude= 'y')




  ###Archivos
  ###Archivos

  datadir='../Utrack1.0/'
  arr=sorted(os.listdir(datadir))
  print(arr)



  #####Flujos superficiales
  file2='../ERA5/P_ERA5_MA08-17.nc'
  file3='../ERA5/E_ERA5_MA08-17.nc'
  ds2=xr.open_dataset(file2)
  ds3=xr.open_dataset(file3)

  ###Variables
  pcp=ds2['tp']
  evp=abs(ds3['e'])

  ####Aplicar máscara de la cuenca a datos de lluvia
  mask = np.empty(np.shape(pcp[0]))
  #mask[:]=np.nan
  mask[pshed_mask == 1] = 1
  
  ls_mask = np.empty(np.shape(pcp[0]))
  ls_mask[:]=1
  ls_mask[xr_lsm == 0] = 0
  pcp=pcp*mask
  #evp=evp*ls_mask





  #####Extraer coordenadas para la mascara de la cuenca
  id_la=list()
  id_lo=list()
  la=list()
  lo=list()


  for i in range(0,len(latitude)):
    for j in range(0,len(longitude)):
      if (np.isnan(pcp[0,i,j].values)) != True:
        lat_v=pcp[0,i,j].latitude.values
        lon_v=pcp[0,i,j].longitude.values
        id_la.append(list(latitude).index(lat_v))
        id_lo.append(list(longitude).index(lon_v))
        la.append(lat_v)
        lo.append(lon_v)

  reciclaje_m = np.empty([12,180,360])
  #reciclaje_m[:] = np.nan
  print('Numero de puntos', len(id_la))

  c=0
  for i in arr:
    file=(datadir+i)
    print(file)
    with xr.open_dataset(file) as ds:
      rec_e=ds['moisture_flow']
      #print(rec_e)
      latt=ds['targetlat']
      lont=ds['targetlon']

      #print(np.shape(reciclaje))
      reciclaje = np.empty([len(id_la),180,360])
      #reciclaje[:] = np.nan
      cc=0
      for i,j in zip(id_la,id_lo):
          recic=np.exp(rec_e[:,:,i,j]*(-0.1))
          recic_e=recic/(recic.sum(axis=0,skipna=True).sum(axis=0,skipna=True))
          #recic_e=recic_e*ls_mask[0:180,:]
          #print(i,j)
          reciclaje[cc,:,:]=recic_e
          cc=cc+1
      recic_f=np.nansum(reciclaje,axis=0)###fraction Ec total por mes a la cuenca
    reciclaje_m[c,:,:]=recic_f*evp[c,0:180,:]###flujo Ec a la cuenca para los 12 meses
    c=c+1



  rec_an=np.nansum(reciclaje_m,axis=0)####Flujo Ec total anual para la cuenca
  pan_basin=np.nansum(pcp)####P anual acumulada en la cuenca
  rec_an2=rec_an*ls_mask[0:180,:]
  prec_r=np.nansum(rec_an2)*100/np.nansum(rec_an)
  print('@@@@@@@@@@@@@@@@@@@')
  print('Recycling Ratio = ', prec_r)


  df_data.loc[df_data['ID']==(basin),'RecyclinR']=(prec_r)


  cats=[0,10,20,30,40,50,60,70,80,90,100]

  if len(id_la) != 0:
    var, accum, cat=(acumulator2(np.asarray(rec_an2),latitude1,longitude1,cats, ls_mask))
    print(np.min(cat))
    print(np.max(cat))



    x_shift,y_shift = np.meshgrid(longitude,latitude)
    x = x_shift 
    y = y_shift
    clevs2 = np.linspace(1,5,6)
    tcl=['100','90','80','70','60','50','40','30','20','10','0']
    #cs = axs1.contourf(x,y,cat,clevs2, cmap=cmap)


   

    #cat[cat > 5] = 1
    #cat[cat <=5] = np.nan



    longitude2=list()

    for i in longitude:
      if i >=0 and i<180:
        longitude2.append(i)
      if i>=180 and i<360:
        longitude2.append(i-360)


    longitude2=sorted(longitude2)
    catf=conv(cat)


    print(np.shape(latitude))
    print(np.shape(longitude))
    print(np.shape(cat))
    save_nc('./MIPS/MIP'+str(basin)+'Final.nc',latitude[0:180],longitude2,catf)
    del pcp, mask,pshed_mask,polygon,xr_lsm



df_data.to_csv('Recycling1-99.csv')
print(len(id_lo))

end1 = timer()
print ('The total runtime is',(end1-start1),' seconds.')
