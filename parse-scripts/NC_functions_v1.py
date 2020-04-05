#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thurs March 12 17:22:21 2020

@author: Heather Guy

Functions to generate netCDF files. 

"""

def NC_SpecificVariables(fn_nc, var, np):
    """
    Collects specific variables from standard excel file and
    writes to file. 
    
    Parameters:
        fn_nc : NetCDF file to write on. 
        var :   Specific variables DataFrame from standard
                excel file. Read by pd.read_excel()
        np:     numpy
    
    """
    for i in range(0,len(var['Variable'].dropna())):
        if i+1 == len(var['Variable'].dropna()):
            att_end = var.index[-1]
        else: 
            att_end = var['Variable'].dropna().index[i+1] 
           
        att_list = var['Attribute'][var['Variable'].dropna().index[i]:att_end+1].dropna()
        val_list = var['Value'][att_list.index]       
       
        if len(att_list[att_list=='_FillValue'])==1:
            fv=float(val_list[att_list[att_list=='_FillValue'].index])
        else:
            fv=None
       
        varn = fn_nc.createVariable(var['Variable'].dropna().iloc[i], np.float64, val_list.iloc[1].split(', '),fill_value=fv)

        for j in range(0,len(att_list)):
            if att_list.iloc[j][0]=='_':
                continue
            else:
                varn.setncattr(att_list.iloc[j],val_list.iloc[j])
                
    return


def NC_CommonVariables(fn_nc, time_list, np):
    """
    Writes common variables to netCDF file. 
    
    Parameters:
        fn_nc :     NetCDF file to write on. 
        time_list:  Pandas date_range for the time dimension
        np:         numpy
    
    """
    import netCDF4
    #time
    times = fn_nc.createVariable('time', np.float64, ('time',))
    #variable attributes
    times.type = 'float64'
    times.units = 'seconds since 1970-01-01 00:00:00 UTC'
    times.standard_name = 'time'
    times.long_name = 'Time (seconds since 1970-01-01 00:00:00)'
    times.axis = 'T'
    times.valid_min = np.float64(netCDF4.date2num(min(time_list),units=times.units))
    times.valid_max = np.float64(netCDF4.date2num(max(time_list),units=times.units))
    times.calendar = 'standard'
    #write data
    times[:] = np.float64(netCDF4.date2num(time_list.to_list(),units=times.units))
 
    #year
    years = fn_nc.createVariable('year', np.int32, ('time',))
    #variable attributes
    years.type = 'int32'
    years.units = '1'
    years.long_name = 'Year'
    years.valid_min = np.int32(1900)
    years.valid_max = np.int32(2100) 
    #write data
    years[:] = np.int32(time_list.year.to_numpy())

    #month
    months = fn_nc.createVariable('month', np.int32, ('time',))
    #variable attributes
    months.type = 'int32'
    months.units = '1'
    months.long_name = 'Month'
    months.valid_min = np.int32(1)
    months.valid_max = np.int32(12) 
    #write data
    months[:] = np.int32(time_list.month.to_numpy())
   
    #day
    days = fn_nc.createVariable('day', np.int32, ('time',))
    #variable attributes
    days.type = 'int32'
    days.units = '1'
    days.long_name = 'Day'
    days.valid_min = np.int32(1)
    days.valid_max = np.int32(31)
    #write data
    days[:] = np.int32(np.int32(time_list.day.to_numpy()))
   
    #hour
    hours = fn_nc.createVariable('hour', np.int32, ('time',))
    #variable attributes
    hours.type = 'int32'
    hours.units = '1'
    hours.long_name = 'Hour'
    hours.valid_min = np.int32(0)
    hours.valid_max = np.int32(23) 
    #write data
    hours[:] = np.int32(time_list.hour.to_numpy())
   
    #minute
    minutes = fn_nc.createVariable('minute', np.int32, ('time',))
    #variable attributes
    minutes.type = 'int32'
    minutes.units = '1'
    minutes.long_name = 'Minute'
    minutes.valid_min = np.int32(0)
    minutes.valid_max = np.int32(59)  
    #write data
    minutes[:] = np.int32(time_list.minute.to_numpy())
   
    #second
    seconds = fn_nc.createVariable('second', np.float32, ('time',))
    #variable attributes
    seconds.type = 'float32'
    seconds.units = '1'
    seconds.long_name = 'Second'
    seconds.valid_min = np.float32(0)
    seconds.valid_max = np.float32(59) 
    #write data
    seconds[:] = np.int32(time_list.second.to_numpy())
   
    #doy
    doys = fn_nc.createVariable('day_of_year', np.float32, ('time',))
    #variable attributes
    doys.type = 'float32'
    doys.units = '1'
    doys.long_name = 'Day of Year'
    doys.valid_min = np.float32(0)
    doys.valid_max = np.float32(365)
    #write data
    doys[:] = np.float32(np.asarray([time_list[i].timetuple().tm_yday for i in range(0,len(time_list))]))
    
    lats = fn_nc.createVariable('latitude', np.float32, ('latitude',))
    #variable attributes
    lats.type = 'float32'
    lats.units = 'degree_north'
    lats.long_name = 'Latitude'
    lats[:] = [72.572962]
   
    lons = fn_nc.createVariable('longitude', np.float32, ('longitude',))
    #variable attributes
    lons.type = 'float32'
    lons.units = 'degree_east'
    lons.long_name = 'Longitude'
    lons[:] = [-38.470361]
   
    return
      
   
def NC_Dimensions(fn_nc, len_time,index=False):
    """
    Writes dimensions to netCDF file
    
    Parameters:
        fn_nc :     NetCDF file to write on. 
        len_time :  Length of time dimension
        index:      length of optional index dimension
    
    """
    fn_nc.createDimension('time', len_time )
    fn_nc.createDimension('latitude', 1)
    fn_nc.createDimension('longitude', 1) 
    if index:
        index = fn_nc.createDimension('index', index)
    
    return


def NC_Global_Attributes(fn_nc, meta, start_date,end_date):
    """
    Writes global attributes to netCDF file
    
    Parameters:
        fn_nc :     NetCDF file to write on.
        meta :      DataFrame of global attributes from standard
                    excel file. Read by pd.read_excel()
        start_date: Start datetime of data in file
        end_date:   end datetime of data in file

    """
    from datetime import datetime
    import numpy as np
    name = meta.loc[:, 'Name':'Name':1].values
    exp = meta.loc[:, 'Example':'Example':1].values
    pos = exp[34]
    pos = pos[0]
    ix1 = pos.find('N')
    ix2 = pos.find(',')
    ix3 = pos.find('E')
    lat = np.float32(pos[0:ix1])
    lon = np.float32(pos[ix2+1:ix3])
    pos = exp[35]
    pos = pos[0]
    ix1 = pos.find('m')
    base_height = np.float32(pos[0:ix1])
   
    for i in range(40):
       msg1 = np.array(name[i])
       msg2 = np.array(exp[i])
       fn_nc.setncattr(msg1[0], msg2[0])
   
    fn_nc.last_revised_date = datetime.utcnow().isoformat()  
    fn_nc.time_coverage_start = start_date.isoformat()
    fn_nc.time_coverage_end = end_date.isoformat()
   
    return

