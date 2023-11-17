#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# author(s): michael.r.gallagher@noaa.gov

# Heather Guy has adapted some of these functions from Michael so they fit with the 
# ICECAPS-ACE data workflow. 

import os, re, glob, fnmatch, argparse # builtin python libs
import datetime as dt
import xarray   as xr # installed python libs
import numpy    as np
import pandas   as pd

def parse_simba_serialstream(sample_dict, tstamp_dict,num_temp_sensors):
    """
    
    Take the raw data put into a dictionary by extract_serial_samples() and
    parse through all the strings to get meaningful data. We pull out the max
    possible amount of information that might or might not be useful. 
    
    Args:
        sample_dict: (dictionary) timestamped [key] datastreams [value] from extract
        tstamp_dict: (dictionary) line timestamps from serial stream to know when each line came
    
    
    Returns:
        parsed_samples: (dictionary) timestamped [key] dictionary of measurements vars  [value]
    
    """
    simba_header_keywords = {
        'sequence_number' : re.compile('sequen\D*(?P<match>[0-9]+\.*[0-9]*)', re.IGNORECASE), 
        'battery_voltage' : re.compile('diode\D*(?P<match>[0-9]+\.*[0-9]*)', re.IGNORECASE), 
    }
    
    good_file_regex   = re.compile(fnmatch.translate(f"*sd card*"), re.IGNORECASE)
    temps_start_regex = re.compile(fnmatch.translate(f"*chain sensor temperatures*"), re.IGNORECASE)
    
    number_bad_samples = 0
    parsed_samples     = {} # a dict of dicts, keys are first timestamp of entire data sample 
    for tstamp_sample, data_sample in sample_dict.items(): 
        
        # a dictionary of the measurements in this sample, this is used later to fill
        # the xarray dataset, so the variable names must match up
        current_sample = {} # dict that goes in the dict of dicts, so many dicts
        cs = current_sample # shorthand so code is more readable
        
        # keep some timestamp information about this sample to write out
        cs['sample_start'] = tstamp_dict[tstamp_sample][0]
        cs['sample_end']   = tstamp_dict[tstamp_sample][-1]
        cs['sample_span']  = (cs['sample_end']-cs['sample_start'])#.total_seconds()
        
        if len(data_sample)<2:
            print('bad data sample')
            print(data_sample)
            continue
        
        if re.match(good_file_regex, data_sample[0]) and re.match(good_file_regex, data_sample[1]):
            try:
            
                for iline, data_line in enumerate(data_sample):
                    if re.match(temps_start_regex, data_line):
                        raw_temp_data = data_sample[iline+1:iline+1+num_temp_sensors]
                        data_tail     = data_sample[iline+2+num_temp_sensors:-1]
                        data_header   = data_sample[0:iline]
                        break
                
                for hline in data_header:
                    for var, regex in simba_header_keywords.items():
                        m = re.search(regex, hline)
                        if m: cs[var] = m.group('match')
                
                cs['temperature'] = [t.split()[1] for t in raw_temp_data]
                
                #for dt in raw_temp_data: print(dt)
                #for dt in data_tail: print(dt)
                #for dt in data_header: print(dt)
                #for k, v in cs.items(): print(f"{k} : {v}")
                
                parsed_samples[tstamp_sample] = current_sample
                
            except: 
                print(f"!!! sample from {tstamp_sample} failed !!!")
                number_bad_samples+=1
                continue
        else: 
            print(f"!!! there was a bad sample at... {tstamp_sample} !!!")
            number_bad_samples+=1
            
    
    perc_bad   = np.round((number_bad_samples/len(sample_dict))*100, 1)
    date_range = f"{list(sample_dict.keys())[0]} ---> {tstamp_dict[tstamp_sample][-1].replace(microsecond=0)}"
    
    # bad should print to stderr............................................
    if perc_bad>0: print(f"!!! {perc_bad}% of SIMBA samples {date_range} were BAD!!!")
    else:          print(f"... all samples were good for {date_range}")
    
    return parsed_samples # this dict should be sorted before returning

def get_file_list_serial(today, dpath, file_str='simba'):
    """
    Uses the "standardized" naming from the daily serial stream
    script to find corresponding files.

    Args:
        today    : (datetime) the date that you would like to find data for
        dpath    : (str) path to find the data in
        file_str : (str) a prefix for narrowing down file names

    Returns:
        return: list of files
    """
    # gallagher naming convention
    #file_regex = re.compile(fnmatch.translate(f"*{file_str}*{today.strftime('%Y%m%d')}*"), re.IGNORECASE)
    #file_list = [j for j in glob.glob(f"{dpath}/*") if re.match(file_regex, j)]
    
    # guy naming convention
    #if len(file_list) == 0:
    file_regex = re.compile(fnmatch.translate(f"*{today.strftime('%y%m')}*.simba"), re.IGNORECASE)
    file_list = [j for j in glob.glob(f"{dpath}/*") if re.match(file_regex, j)]
        
    # sort file_list by time:
    file_dates = []
    for dfile in file_list: file_dates.append(dt.datetime.strptime(os.path.basename(dfile)[0:11], "%y%m%d_%H%M"))
    sorted_inds = np.argsort(file_dates)
    sorted_files = [file_list[i] for i in sorted_inds]
    
    if len(sorted_files) == 0: raise FileNotFoundError("!!! Didn't find any simba data files matching 'conventions'")
    
    for sfile in sorted_files: print(f" ... {sfile}")
    return sorted_files

def extract_serial_samples(file_list, first_char='['):
    """

    Take a dumped serial file and convert it to a well formatted but-still-raw
    time-keyed dictionary that you might process and convert to say, a pd.DataFrame

    Args:
        file_list: (list(str)) filenames that you want to pull data out of
        first_char: (str)      streaming serial is ugly, this is the character that
                               marks the beginning of a 'proper' line from the serial stream 

    Returns:
        fsamples_dict: dictionary where datetime.datetime of measurement is the key
                       and the value is a dictionary that provides the var name str [key]
                       and the measurement list(float) [value]... a dict of dicts

    """
    alldata_dict = {} # start with kludgy list of strings line-by-line 
    
    # loop through files, open, and clean up blank lines and the first missed few
    # lines that happen cause a raw serial stream isn't the best format
    for curr_file in file_list: 
        
        with open(curr_file,"r", encoding="utf-8", errors="ignore", newline='\n') as f:
            
            # the if==first_char is because minicom fails to timestamp the first lines
            # and occasionally there is garbled stuff, garble garble
            
            file_data = []
            for l in f.readlines():
                try: 
                    #print(repr(l[0]), first_char) # this was really annoying
                    if l[0]==first_char: file_data.append(l)
                except: 
                    print(l)
            
            # try:    file_data = [x for x in f.readlines() if x[0]==first_char]
            # except: raise("Bad file - couldn't read, this shouldn't happen...")
            
        alldata_dict[curr_file] = [item for item in file_data if item] # drop empty lines, put in dict
    
    # loop over all lines and strip out minicom timestamps, aggregate into lines into "samples"
    # aka unique datapoints in timeseries, as long as timestamps of all lines are all in a window of:
    sample_time_window = 3  # minutes
    fsamples_dict      = {} # this gets filled with all the data "samples" and returned 
    
    first_timestamp = dt.datetime(1970,1,1)
    
    for curr_file, file_data in alldata_dict.items():
        
        sample_dict    = {} # this is where the data will live, keys are timestamps of the first line 
        timestamp_dict = {} # but also keep all line timestamps here for debugging
        
        first_failure = False
        for data_line in file_data:
            
            # rip out timestamp from first N characters of the line of data
            # the timestamp format is specific to the minicom timestamp, it's always
            # formatted the same so we can just use some not-flexible indexing instead of regex
            try: 
                if first_char == '[':
                    sample_time = dt.datetime.strptime(data_line[1:20], '%Y-%m-%d %H:%M:%S')
                else:
                    sample_time = dt.datetime.strptime(data_line[0:23], '%Y %m %d %H %M %S.%f')
            
                if (sample_time - first_timestamp) > dt.timedelta(minutes=sample_time_window):
                    first_timestamp                 = sample_time
                    sample_dict[first_timestamp]    = []
                    timestamp_dict[first_timestamp] = []
            except:
                print(f'Cannot read time stamps {curr_file} is broken.. moving on')
                break
            
            # sanity check, hasn't happened yet
            if first_timestamp > sample_time: 
                if first_failure: 
                    print(f'You went backwards in time... wtf. {curr_file} is very broken.. moving on')
                    break
                    
                else:# accounts for the random dropping of old lines into heather's new files.. should only be one line
                    first_failure = True
                    continue 
                    
            data = data_line[26:-1] # all the data that's not part of timestamp
            if data != '':          # and don't keep the line if it's blank... again
                try: sample_dict[first_timestamp].append(data)
                except: 
                    if 'reinit' not in data: 
                        print(f"!!! dropping !!!\n {data}") # sometimes heather's files have meaningless split lines
                        # but it seemsto always be the useless Filesystemreinitialised message, print if not
                    continue;
                    
                if first_char == '[':
                    timestamp_dict[first_timestamp] .append(dt.datetime.strptime(data_line[1:24]+'000',
                                                                               '%Y-%m-%d %H:%M:%S.%f'))
                else:
                    timestamp_dict[first_timestamp] .append(dt.datetime.strptime(data_line[0:23],
                                                                                 '%Y %m %d %H %M %S.%f'))
        if len(sample_dict)>0: 
            fsamples_dict[curr_file] = (sample_dict, timestamp_dict)
        else: print(f"!!! {curr_file} was found empty !!!")
        
    return fsamples_dict