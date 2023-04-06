import os
import sys
import wrf
import time
import xarray as xr
import numpy as np
from wrf import getvar
import netCDF4 as nc
from netCDF4 import Dataset

os.environ["HDF5_USE_FILE_LOCKING"] = 'FALSE'

def unix_time(dt):
    # 转换成时间数组
    timeArray = time.strptime(dt, "%Y%m%d%H%M%S")
    # 转换成时间戳
    timestamp = int(time.mktime(timeArray))
    return timestamp

def write_nc_data(fileoutput,BADPT,times,nlevel,nlat,nlon,QRNEW,QROLD,QSNEW,QSOLD,QGNEW,QGOLD,UNEW,UOLD,VNEW,VOLD,TNEW,TOLD,QVNEW,QVOLD,REFNEW,REFOLD,QCNEW,QCOLD,QINEW,QIOLD):

    #---define dimensions
    ncfile = Dataset(fileoutput,'w',format='NETCDF4')
    ncfile.createDimension('south_north',nlat)
    ncfile.createDimension('west_east',nlon)
    ncfile.createDimension('south_north_stag',nlat+1)
    ncfile.createDimension('west_east_stag',nlon+1)
    ncfile.createDimension('bottom_top',nlevel)
    ncfile.createDimension('Time')
    ncfile.createDimension('DateStrLen',19)

    #---define variables
    ncfile.createVariable('Times','c',('Time','DateStrLen'))
    #ncfile.createVariable('REF_OLD','f4',('Time','bottom_top','south_north','west_east'))
    #ncfile.createVariable('REF_NEW','f4',('Time','bottom_top','south_north','west_east'))
    ncfile.createVariable('QR_NDG_OLD','f4',('Time','bottom_top','south_north','west_east'))
    ncfile.createVariable('QR_NDG_NEW','f4',('Time','bottom_top','south_north','west_east'))
    ncfile.createVariable('QS_NDG_OLD','f4',('Time','bottom_top','south_north','west_east'))
    ncfile.createVariable('QS_NDG_NEW','f4',('Time','bottom_top','south_north','west_east'))
    ncfile.createVariable('QG_NDG_OLD','f4',('Time','bottom_top','south_north','west_east'))
    ncfile.createVariable('QG_NDG_NEW','f4',('Time','bottom_top','south_north','west_east'))
    ncfile.createVariable('U_NDG_OLD','f4',('Time','bottom_top','south_north','west_east_stag'))
    ncfile.createVariable('U_NDG_NEW','f4',('Time','bottom_top','south_north','west_east_stag'))
    ncfile.createVariable('V_NDG_OLD','f4',('Time','bottom_top','south_north_stag','west_east'))
    ncfile.createVariable('V_NDG_NEW','f4',('Time','bottom_top','south_north_stag','west_east'))
    ncfile.createVariable('Q_NDG_OLD','f4',('Time','bottom_top','south_north','west_east'))
    ncfile.createVariable('Q_NDG_NEW','f4',('Time','bottom_top','south_north','west_east'))
    ncfile.createVariable('T_NDG_OLD','f4',('Time','bottom_top','south_north','west_east'))
    ncfile.createVariable('T_NDG_NEW','f4',('Time','bottom_top','south_north','west_east'))
    ncfile.createVariable('QC_NDG_OLD','f4',('Time','bottom_top','south_north','west_east'))
    ncfile.createVariable('QC_NDG_NEW','f4',('Time','bottom_top','south_north','west_east'))
    ncfile.createVariable('QI_NDG_OLD','f4',('Time','bottom_top','south_north','west_east'))
    ncfile.createVariable('QI_NDG_NEW','f4',('Time','bottom_top','south_north','west_east'))
    #---assign values1
    ncfile.variables['Times'][0,:] = times
    #ncfile.variables['REF_OLD'][0,:,:,:] = REFOLD
    #ncfile.variables['REF_NEW'][0,:,:,:] = REFNEW
    ncfile.variables['QR_NDG_OLD'][0,:,:,:] = QROLD
    ncfile.variables['QR_NDG_NEW'][0,:,:,:] = QRNEW
    ncfile.variables['QS_NDG_OLD'][0,:,:,:] = QSOLD
    ncfile.variables['QS_NDG_NEW'][0,:,:,:] = QSNEW
    ncfile.variables['QG_NDG_OLD'][0,:,:,:] = QGOLD
    ncfile.variables['QG_NDG_NEW'][0,:,:,:] = QGNEW

    ncfile.variables['U_NDG_OLD'][0,:,:,:] = UOLD
    ncfile.variables['U_NDG_NEW'][0,:,:,:] = UNEW
    ncfile.variables['V_NDG_OLD'][0,:,:,:] = VOLD
    ncfile.variables['V_NDG_NEW'][0,:,:,:] = VNEW

    ncfile.variables['Q_NDG_OLD'][0,:,:,:] = QVOLD
    ncfile.variables['Q_NDG_NEW'][0,:,:,:] = QVNEW
    ncfile.variables['T_NDG_OLD'][0,:,:,:] = TOLD
    ncfile.variables['T_NDG_NEW'][0,:,:,:] = TNEW
    ncfile.variables['QC_NDG_OLD'][0,:,:,:] = QCOLD
    ncfile.variables['QC_NDG_NEW'][0,:,:,:] = QCNEW
    ncfile.variables['QI_NDG_OLD'][0,:,:,:] = QIOLD
    ncfile.variables['QI_NDG_NEW'][0,:,:,:] = QINEW
    #---variables info
    ncfile.variables['QR_NDG_NEW'].FieldType = 104
    ncfile.variables['QR_NDG_NEW'].MemoryOrder = 'XYZ'
    ncfile.variables['QR_NDG_NEW'].units = 'kg/kg'
    ncfile.variables['QR_NDG_NEW'].stagger = ' '
    ncfile.variables['QR_NDG_NEW'].coordinates = 'XLONG  XLAT'

    ncfile.variables['QR_NDG_OLD'].FieldType = 104
    ncfile.variables['QR_NDG_OLD'].MemoryOrder = 'XYZ'
    ncfile.variables['QR_NDG_OLD'].units = 'kg/kg'
    ncfile.variables['QR_NDG_OLD'].stagger = ' '
    ncfile.variables['QR_NDG_OLD'].coordinates = 'XLONG  XLAT'

    ncfile.variables['QS_NDG_NEW'].FieldType = 104
    ncfile.variables['QS_NDG_NEW'].MemoryOrder = 'XYZ'
    ncfile.variables['QS_NDG_NEW'].units = 'kg/kg'
    ncfile.variables['QS_NDG_NEW'].stagger = ' '
    ncfile.variables['QS_NDG_NEW'].coordinates = 'XLONG  XLAT'

    ncfile.variables['QS_NDG_OLD'].FieldType = 104
    ncfile.variables['QS_NDG_OLD'].MemoryOrder = 'XYZ'
    ncfile.variables['QS_NDG_OLD'].units = 'kg/kg'
    ncfile.variables['QS_NDG_OLD'].stagger = ' '
    ncfile.variables['QS_NDG_OLD'].coordinates = 'XLONG  XLAT'

    ncfile.variables['QG_NDG_NEW'].FieldType = 104
    ncfile.variables['QG_NDG_NEW'].MemoryOrder = 'XYZ'
    ncfile.variables['QG_NDG_NEW'].units = 'kg/kg'
    ncfile.variables['QG_NDG_NEW'].stagger = ' '
    ncfile.variables['QG_NDG_NEW'].coordinates = 'XLONG  XLAT'

    ncfile.variables['QG_NDG_OLD'].FieldType = 104
    ncfile.variables['QG_NDG_OLD'].MemoryOrder = 'XYZ'
    ncfile.variables['QG_NDG_OLD'].units = 'kg/kg'
    ncfile.variables['QG_NDG_OLD'].stagger = ' '
    ncfile.variables['QG_NDG_OLD'].coordinates = 'XLONG  XLAT'

    ncfile.variables['U_NDG_NEW'].FieldType = 104
    ncfile.variables['U_NDG_NEW'].MemoryOrder = 'XYZ'
    ncfile.variables['U_NDG_NEW'].units = 'm s-1'
    ncfile.variables['U_NDG_NEW'].stagger = ' '
    ncfile.variables['U_NDG_NEW'].coordinates = 'XLONG  XLAT'

    ncfile.variables['U_NDG_OLD'].FieldType = 104
    ncfile.variables['U_NDG_OLD'].MemoryOrder = 'XYZ'
    ncfile.variables['U_NDG_OLD'].units = 'm s-1'
    ncfile.variables['U_NDG_OLD'].stagger = ' '
    ncfile.variables['U_NDG_OLD'].coordinates = 'XLONG  XLAT'

    ncfile.variables['V_NDG_NEW'].FieldType = 104
    ncfile.variables['V_NDG_NEW'].MemoryOrder = 'XYZ'
    ncfile.variables['V_NDG_NEW'].units = 'm s-1'
    ncfile.variables['V_NDG_NEW'].stagger = ' '
    ncfile.variables['V_NDG_NEW'].coordinates = 'XLONG  XLAT'

    ncfile.variables['V_NDG_OLD'].FieldType = 104
    ncfile.variables['V_NDG_OLD'].MemoryOrder = 'XYZ'
    ncfile.variables['V_NDG_OLD'].units = 'm s-1'
    ncfile.variables['V_NDG_OLD'].stagger = ' '
    ncfile.variables['V_NDG_OLD'].coordinates = 'XLONG  XLAT'

    ncfile.variables['Q_NDG_NEW'].FieldType = 104
    ncfile.variables['Q_NDG_NEW'].MemoryOrder = 'XYZ'
    ncfile.variables['Q_NDG_NEW'].units = 'kg/kg'
    ncfile.variables['Q_NDG_NEW'].stagger = ' '
    ncfile.variables['Q_NDG_NEW'].coordinates = 'XLONG  XLAT'

    ncfile.variables['Q_NDG_OLD'].FieldType = 104
    ncfile.variables['Q_NDG_OLD'].MemoryOrder = 'XYZ'
    ncfile.variables['Q_NDG_OLD'].units = 'kg/kg'
    ncfile.variables['Q_NDG_OLD'].stagger = ' '
    ncfile.variables['Q_NDG_OLD'].coordinates = 'XLONG  XLAT'

    ncfile.variables['T_NDG_NEW'].FieldType = 104
    ncfile.variables['T_NDG_NEW'].MemoryOrder = 'XYZ'
    ncfile.variables['T_NDG_NEW'].units = 'K'
    ncfile.variables['T_NDG_NEW'].stagger = ' '
    ncfile.variables['T_NDG_NEW'].coordinates = 'XLONG  XLAT'

    ncfile.variables['T_NDG_OLD'].FieldType = 104
    ncfile.variables['T_NDG_OLD'].MemoryOrder = 'XYZ'
    ncfile.variables['T_NDG_OLD'].units = 'K'
    ncfile.variables['T_NDG_OLD'].stagger = ' '
    ncfile.variables['T_NDG_OLD'].coordinates = 'XLONG  XLAT'

    ncfile.variables['QC_NDG_NEW'].FieldType = 104
    ncfile.variables['QC_NDG_NEW'].MemoryOrder = 'XYZ'
    ncfile.variables['QC_NDG_NEW'].units = 'kg/kg'
    ncfile.variables['QC_NDG_NEW'].stagger = ' '
    ncfile.variables['QC_NDG_NEW'].coordinates = 'XLONG  XLAT'

    ncfile.variables['QI_NDG_OLD'].FieldType = 104
    ncfile.variables['QI_NDG_OLD'].MemoryOrder = 'XYZ'
    ncfile.variables['QI_NDG_OLD'].units = 'kg/kg'
    ncfile.variables['QI_NDG_OLD'].stagger = ' '
    ncfile.variables['QI_NDG_OLD'].coordinates = 'XLONG  XLAT'


    #ncfile.variables['REF_NEW'].FieldType = 104
    #ncfile.variables['REF_NEW'].MemoryOrder = 'XYZ'
    #ncfile.variables['REF_NEW'].units = 'dBZ'
    #ncfile.variables['REF_NEW'].stagger = ' '
    #ncfile.variables['REF_NEW'].coordinates = 'XLONG  XLAT'

    #ncfile.variables['REF_OLD'].FieldType = 104
    #ncfile.variables['REF_OLD'].MemoryOrder = 'XYZ'
    #ncfile.variables['REF_OLD'].units = 'dBZ'
    #ncfile.variables['REF_OLD'].stagger = ' '
    #ncfile.variables['REF_OLD'].coordinates = 'XLONG  XLAT'

    #---close ncfile
    ncfile.close()

if __name__ == "__main__":
    file_old = sys.argv[1]
    file_new = sys.argv[2]
    filename = os.path.basename(file_old)
    filepath = os.path.dirname(file_old)
    filestring = filename.split("_")[1:]
    string_lgt = len(filestring)
    fileoutname = 'wrffdda'
    for nn in filestring:
        fileoutname = fileoutname + '_' + nn
    filepath = '/public/home/huozhaoyang/research/wrffdda/fddafile/WRFOUT_TO_WRFFDDA/wrffdda_noise/'
    fileoutput = filepath+'/'+fileoutname
    print('file_old:' + os.path.basename(file_old))
    print('file_new:' + os.path.basename(file_new))
    print('fileoutput:'+fileoutput)


    filedata_old=Dataset(file_old)
    filedata_new=Dataset(file_new)
    tc = getvar(filedata_old, "tc", meta=False)
    nlevel = tc.shape[0]
    nlat = tc.shape[1]
    nlon = tc.shape[2]

    times = getvar(filedata_old, "Times",meta=False)
    timestr = str(times)
    TIME = timestr[0:4] +'-'+ timestr[5:7] +'-'+ timestr[8:10] +'_'+ timestr[11:13] +':'+ timestr[14:16] +':'+ timestr[17:19]
    print(TIME)

    BADPT=-999.0
    REFNEW = getvar(filedata_new, "REFL_10CM",meta=False)
    REFOLD = getvar(filedata_old, "REFL_10CM",meta=False)

    QRNEW = getvar(filedata_new, "QRAIN",meta=False)
    QROLD = getvar(filedata_old, "QRAIN",meta=False)
    QSNEW = getvar(filedata_new, "QSNOW",meta=False)
    QSOLD = getvar(filedata_old, "QSNOW",meta=False)
    QGNEW = getvar(filedata_new, "QGRAUP",meta=False)
    QGOLD = getvar(filedata_old, "QGRAUP",meta=False)

    UNEW = getvar(filedata_new, "U",meta=False)
    UOLD = getvar(filedata_old, "U",meta=False)
    VNEW = getvar(filedata_new, "V",meta=False)
    VOLD = getvar(filedata_old, "V",meta=False)

    TNEW = getvar(filedata_new, "T",meta=False)
    TOLD = getvar(filedata_old, "T",meta=False)
    QVNEW = getvar(filedata_new, "QVAPOR",meta=False)
    QVOLD = getvar(filedata_old, "QVAPOR",meta=False)


    QCNEW = getvar(filedata_new, "QCLOUD",meta=False)
    QCOLD = getvar(filedata_old, "QCLOUD",meta=False)
    QINEW = getvar(filedata_new, "QICE",meta=False)
    QIOLD = getvar(filedata_old, "QICE",meta=False)


    REFN=np.array([[[BADPT for i in range(nlon)] for j in range(nlat)] for k in range(nlevel)])
    REFO=np.array([[[BADPT for i in range(nlon)] for j in range(nlat)] for k in range(nlevel)])

    QRN=np.array([[[BADPT for i in range(nlon)] for j in range(nlat)] for k in range(nlevel)])
    QSN=np.array([[[BADPT for i in range(nlon)] for j in range(nlat)] for k in range(nlevel)])
    QGN=np.array([[[BADPT for i in range(nlon)] for j in range(nlat)] for k in range(nlevel)])
    QRO=np.array([[[BADPT for i in range(nlon)] for j in range(nlat)] for k in range(nlevel)])
    QSO=np.array([[[BADPT for i in range(nlon)] for j in range(nlat)] for k in range(nlevel)])
    QGO=np.array([[[BADPT for i in range(nlon)] for j in range(nlat)] for k in range(nlevel)])

    UN=np.array([[[BADPT for i in range(nlon)] for j in range(nlat+1)] for k in range(nlevel)])
    VN=np.array([[[BADPT for i in range(nlon+1)] for j in range(nlat)] for k in range(nlevel)])
    UO=np.array([[[BADPT for i in range(nlon)] for j in range(nlat+1)] for k in range(nlevel)])
    VO=np.array([[[BADPT for i in range(nlon+1)] for j in range(nlat)] for k in range(nlevel)])

    QVN=np.array([[[BADPT for i in range(nlon)] for j in range(nlat)] for k in range(nlevel)])
    QVO=np.array([[[BADPT for i in range(nlon)] for j in range(nlat)] for k in range(nlevel)])
    TN=np.array([[[BADPT for i in range(nlon)] for j in range(nlat)] for k in range(nlevel)])
    TO=np.array([[[BADPT for i in range(nlon)] for j in range(nlat)] for k in range(nlevel)])

    QCN=np.array([[[BADPT for i in range(nlon)] for j in range(nlat)] for k in range(nlevel)])
    QCO=np.array([[[BADPT for i in range(nlon)] for j in range(nlat)] for k in range(nlevel)])
    QIN=np.array([[[BADPT for i in range(nlon)] for j in range(nlat)] for k in range(nlevel)])
    QIO=np.array([[[BADPT for i in range(nlon)] for j in range(nlat)] for k in range(nlevel)])

    '''
for i in range(nlevel):
        for j in range(nlat):
            for k in range(nlon):
                QRN[i,j,k] = QRNEW[i,j,k]
                QRO[i,j,k] = QROLD[i,j,k]
                QSN[i,j,k] = QSNEW[i,j,k]
                QSO[i,j,k] = QSOLD[i,j,k]
                QGN[i,j,k] = QGNEW[i,j,k]
                QGO[i,j,k] = QGOLD[i,j,k]
    '''
    QRN = QRNEW
    QRO = QROLD
    QSN = QSNEW
    QSO = QSOLD
    QGN = QGNEW
    QGO = QGOLD

    UN = UNEW
    UO = UOLD
    VN = VNEW
    VO = VOLD

    TN = TNEW
    TO = TOLD
    QVN = QVNEW
    QVO = QVOLD

    REFN = REFNEW
    REFO = REFOLD

    QCN = QCNEW
    QCO = QCOLD
    QIN = QINEW
    QIO = QIOLD


    print(np.shape(VN))
    print(np.shape(UN))
    print(np.shape(REFN))

    
    #mask = np.where(REFN < 5, 1, 0)
    #UN[:, :, 1:][mask==1] = -999.0
    #VN[:, 1:, :][mask==1] = -999.0

    #mask = np.where(REFO < 5, 1, 0)
    #UO[:, :, 1:][mask==1] = -999.0
    #VO[:, 1:, :][mask==1] = -999.0
     
    radarmask = '/public/home/huozhaoyang/research/SHPAR/case20200716/wrffdda_true_mask/radar_mask/wrffdda_d03_smallmask' 
    #radarmask = '/public/home/huozhaoyang/research/SHPAR/case20200716/wrffdda_true_mask/wrffdda_d03_radar_mask'
    ncfile = Dataset(radarmask)#,'w',format='NETCDF4')
    RM = ncfile.variables['QR_NDG_OLD']
    RM = RM[0,:,:,:]
    #RM = RM[None]
    RM = np.array(RM)
    mask = np.where(RM == -999, 1, 0)
   

    #mask = np.where(QRN == 0, 1, 0)
    QRN[mask==1] = -999
    #mask = np.where(QRO == 0, 1, 0)
    QRO[mask==1] = -999
    #mask = np.where(QSN == 0, 1, 0)
    QSN[mask==1] = -999
    #mask = np.where(QSO == 0, 1, 0)
    QSO[mask==1] = -999
    #mask = np.where(QGN == 0, 1, 0)
    QGN[mask==1] = -999
    #mask = np.where(QGO == 0, 1, 0)
    QGO[mask==1] = -999

    #mask = np.where(UN == 0, 1, 0)
    UN[:,:,1:][mask==1] = -999
    #mask = np.where(UO == 0, 1, 0)
    UO[:,:,1:][mask==1] = -999
    #mask = np.where(VN == 0, 1, 0)
    VN[:,1:,:][mask==1] = -999
    #mask = np.where(VO == 0, 1, 0)
    VO[:,1:,:][mask==1] = -999

    #mask = np.where(TN == 0, 1, 0)
    TN[mask==1] = -999
    #mask = np.where(TO == 0, 1, 0)
    TO[mask==1] = -999
    #mask = np.where(QVN == 0, 1, 0)
    QVN[mask==1] = -999
    #mask = np.where(QVO == 0, 1, 0)
    QVO[mask==1] = -999
    
    REFN[mask==1] = -999
    REFO[mask==1] = -999
    QCN[mask==1] = -999
    QCO[mask==1] = -999
    QIN[mask==1] = -999
    QIO[mask==1] = -999


    write_nc_data(fileoutput,BADPT,TIME,nlevel,nlat,nlon,QRN,QRO,QSN,QSO,QGN,QGO,UN,UO,VN,VO,TN,TO,QVN,QVO,REFN,REFO,QCN,QCO,QIN,QIO)



