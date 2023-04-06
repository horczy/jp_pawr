
import os
import re
import bz2
import struct
import time
from datetime import datetime
import sys
import numpy as np
from netCDF4 import Dataset
from scipy.interpolate import griddata
from dateutil.relativedelta import relativedelta

class readwasn():
    def __init__(self, filename):
        # print(filename)
        self.filename = filename
        f = bz2.BZ2File(filename, "rb")
        self.ZonName = struct.unpack("12s", f.read(12))
        self.DataName = struct.unpack("38s", f.read(38))  # 数据说明(例如 2008年5月19日雷达三维拼图)38个字节
        self.Flag = struct.unpack("8s", f.read(8))  # 文件标志，"swan"
        self.Version = struct.unpack("8s", f.read(8))  # 数据版本号，"1.0"
        self.year = struct.unpack("H", f.read(2))
        self.month = struct.unpack("H", f.read(2))
        self.day = struct.unpack("H", f.read(2))
        self.hour = struct.unpack("H", f.read(2))
        self.minute = struct.unpack("H", f.read(2))
        self.interval = struct.unpack("H", f.read(2))
        self.XNumGrids = struct.unpack("H", f.read(2))[0]
        self.YNumGrids = struct.unpack("H", f.read(2))[0]
        self.ZNumGrids = struct.unpack("H", f.read(2))[0]
        # print(self.ZNumGrids)
        # self.ZNumGrids=1
        self.RadarCount = struct.unpack("i", f.read(4))  # 拼图雷达数 四个字节
        self.StartLon = struct.unpack("f", f.read(4))[0] # 网格开始经度（左上角） 四个字节
        self.StartLat = struct.unpack("f", f.read(4))[0]  # 网格开始纬度（左上角） 四个字节
        self.CenterLon = struct.unpack("f", f.read(4))[0]  # 网格中心经度 四个字节
        self.CenterLat = struct.unpack("f", f.read(4))[0]  # 网格中心纬度 四个字节
        self.XReso = struct.unpack("f", f.read(4))[0]  # 经度方向分辨率 四个字节
        self.YReso = struct.unpack("f", f.read(4))[0]  # 纬度方向分辨率 四个字节
        self.ZhighGrids = struct.unpack("40f", f.read(40 * 4))
        self.RadarStationName = []
        for i in range(20):
            self.RadarStationName.append(struct.unpack("16s", f.read(16)))
        self.RadarLongitude = struct.unpack("20f", f.read(20 * 4))
        self.RadarLatitude = struct.unpack("20f", f.read(20 * 4))
        self.RadarAltitude = struct.unpack("20f", f.read(20 * 4))
        self.MosaicFlag = struct.unpack("20B", f.read(20))
        f.read(172)
        tempdata = np.frombuffer(f.read(self.ZNumGrids * self.YNumGrids * self.XNumGrids), dtype="B")
        self.data = np.array(tempdata, dtype=int)
        # self.data = (self.data - 66.) / 2.
        self.data.shape = self.ZNumGrids, self.YNumGrids, self.XNumGrids


def read_wrf_latlon(ifile):
 
    ncObj = Dataset(ifile,'r',format='NETCDF3_CLASSIC')

    lat = ncObj.variables['XLAT_M'][:]
    lon = ncObj.variables['XLONG_M'][:]

    ncObj.close()
    return  lat, lon


def read_area_dbz(ifile):
 
    ncObj = Dataset(ifile,'r',format='NETCDF3_CLASSIC')

    lat = ncObj.variables['lat'][:]
    lon = ncObj.variables['lon'][:]
    dbz = ncObj.variables['com_ref'][:][:][:]

    ncObj.close()
    return  dbz, lat, lon 


def write_cross_data1(outfile,nlat,nlon,nlev,grid_lat,grid_lon,ext_var):

    #---define dimensions
    ncfile = Dataset(outfile,'w',format='NETCDF3_CLASSIC')    
    ncfile.createDimension('south_north',nlat)
    ncfile.createDimension('west_east',nlon)
    ncfile.createDimension('bottom_top',nlev)
    
    #---define variables
    lev     = ncfile.createVariable('lev','i4',('bottom_top'))
    lat     = ncfile.createVariable('lat','f4',('south_north'))
    lon     = ncfile.createVariable('lon','f4',('west_east'))
    ele     = ncfile.createVariable('dbz','i4',('bottom_top','south_north','west_east'))

    #---assign values
    lev[:] = np.array([500,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000,7000,8000,9000,10000,12000,14000,15500,17000,19000])
    lon[:] = grid_lon
    lat[:]  = grid_lat
    ele[:,:,:] = ext_var
    
    #---close ncfile
    ncfile.close()


def write_cross_data2(outfile,nlat,nlon,nlev,grid_lat,grid_lon,ext_var):

    #---define dimensions
    ncfile = Dataset(outfile,'w',format='NETCDF3_CLASSIC')    
    ncfile.createDimension('south_north',nlat)
    ncfile.createDimension('west_east',nlon)
    ncfile.createDimension('lev',nlev)
    
    #---define variables
    ncfile.createVariable('lev','i4',('lev'))
    ncfile.createVariable('lat','f4',('south_north','west_east'))
    ncfile.createVariable('lon','f4',('south_north','west_east'))
    ncfile.createVariable('dbz','i4',('lev','south_north','west_east'))

    #---assign values
    lev = np.array([500,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000,7000,8000,9000,10000,12000,14000,15500,17000,19000])
    ncfile.variables['lev'][:] = lev
    ncfile.variables['lon'][:,:] = grid_lon
    ncfile.variables['lat'][:,:] = grid_lat
    ncfile.variables['dbz'][:,:,:] = ext_var
    
    #---close ncfile
    ncfile.close()


def write_cross_data3(outfile,time_s,nlat,nlon,nlev,grid_lat,grid_lon,ext_var,ext_var2):

    #---define dimensions
    ncfile = Dataset(outfile,'w',format='NETCDF4')    
    ncfile.createDimension('south_north',nlat)
    ncfile.createDimension('west_east',nlon)
    ncfile.createDimension('lev',nlev)
    ncfile.createDimension('time',1)
    # ncfile.createDimension('DateStrLen',None)
    
    #---define variables
    ncfile.createVariable('Times','f4',('time'))
    ncfile.createVariable('lev','i4',('lev'))
    ncfile.createVariable('lat','f4',('south_north','west_east'))
    ncfile.createVariable('lon','f4',('south_north','west_east'))
    ncfile.createVariable('dbz','i4',('time','lev','south_north','west_east'))
    ncfile.createVariable('maxdbz','i4',('time','lev','south_north','west_east'))

    #---assign values
    lev = np.array([500,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000,7000,8000,9000,10000,12000,14000,15500,17000,19000])
    ncfile.variables['Times'][0] = unix_time(time_s)
    ncfile.variables['lev'][:] = lev
    ncfile.variables['lon'][:,:] = grid_lon
    ncfile.variables['lat'][:,:] = grid_lat
    ncfile.variables['dbz'][0,:,:,:] = ext_var
    ncfile.variables['maxdbz'][0,0,:,:] = ext_var2

    ncfile.variables['Times'].description = "seconds since 1970-01-01 00:00:00"
    ncfile.variables['Times'].units = 's'
    ncfile.variables['lev'].description = "height"
    ncfile.variables['lev'].units = "m"
    ncfile.variables['lat'].description = "latitude"
    ncfile.variables['lat'].units = "degrees_north"
    ncfile.variables['lon'].description = "longitude"
    ncfile.variables['lon'].units = "degrees_east"
    ncfile.variables['dbz'].description = 'cappi dbz'
    ncfile.variables['dbz'].units = 'dBZ'
    ncfile.variables['dbz'].description = 'Max dbz'
    ncfile.variables['dbz'].units = 'dBZ'
    
    #---close ncfile
    ncfile.close()


def unix_time(dt):
    # 转换成时间数组
    timeArray = time.strptime(dt, "%Y%m%d%H%M%S")
    # 转换成时间戳
    timestamp = int(time.mktime(timeArray))
    return timestamp


if __name__ == "__main__":

    time_date = sys.argv[1]
    file = sys.argv[2]

    bin_path,filename = os.path.split(file)
    # "Z_OTHE_RADAMOSAIC_20200718005400.bin.bz2"
    time_string = filename.split(".")[0].split("_")[3]
    print(time_string)

    dd = readwasn(file)
    ext_var = dd.data
    mask = np.where(ext_var == 0, 1, 0)
    ext_var[mask==1] = -999

    lon0 = dd.StartLon
    lat0 = dd.StartLat

    nlon = dd.XNumGrids
    nlat = dd.YNumGrids
    nlev = dd.ZNumGrids

    dx = dd.XReso
    dy = dd.YReso

    lon1 = lon0 + (nlon - 1) * dx
    lat1 = lat0 - (nlat - 1) * dy

    grid_lon = np.linspace(lon0, lon1, nlon)
    grid_lat = np.linspace(lat0, lat1, nlat)

    lat = grid_lat[:]
    lon = grid_lon[:]
    dbz = ext_var[:,:,:]

    num = len(lat)*len(lon)
    
    points = np.zeros([num,2])
    i=0
    for ilat in lat:
        for ilon in lon:
            points[i,0]=ilat
            points[i,1]=ilon
            i+=1

    lat1,lon1 = read_wrf_latlon("/public/home/premopr/data/SWAN_radar_process/GEJBGD/cappi_to_wrf/geo_em.d04.nc")
    nlat = np.shape(lat1)
    nlon = np.shape(lon1)

    new_cappi = np.zeros((21, nlat[1], nlat[2]))
    maxdbz = np.zeros((nlat[1], nlat[2]))

    for il in range(21):
        print("-------------" + str(il) )
        cappi_ref = dbz[il,:,:]
        values = np.reshape(cappi_ref, (num,))
        interp_data = griddata( points, values, (lat1, lon1), method='nearest')
        interp_data[ np.isnan(interp_data) ] = -999
        mask = np.where(interp_data < 0, 1, 0)
        mask1 = np.where(interp_data == 1, 1, 0)
        interp_data = (interp_data - 66) / 2
        interp_data[ mask == 1 ] = -999
        interp_data[ mask1 == 1 ] = 0
        new_cappi[il,:,:] = interp_data[0,:,:]

    maxdbz = np.amax(new_cappi, axis=0)
    dir_out = "/public/home/premopr/data/SWAN_radar_process/GEJBGD/outdata/cappi_dbz_wrf/data_d04/%s" %(time_date)
    if not os.path.exists(dir_out):
        os.makedirs(dir_out)
    fileout = os.path.join(dir_out, "cappi_d04_dbz_%s.nc") %(time_string)
    print(fileout)
    write_cross_data3(fileout,time_string,nlat[1],nlon[2],21,lat1[0,:,:],lon1[0,:,:],new_cappi,maxdbz)

                



