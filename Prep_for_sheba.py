#!/usr/bin/env python

#RMS 2018 
#Prepare waveform files for shear wave splitting with sheba

import obspy as op
import sys
import glob
import pandas as pd
import os
from obspy.taup import TauPyModel
from obspy.geodetics import locations2degrees

#################################################################################
# HELPER FUNCTIONS
#################################################################################

def filter_traces(trace, basename,station_code,stlat,stlon,evlat,evlon,evdep,win_start,
    win_end,cluster_1,cluster_2,filters=[(0.02,0.1),(0.05,0.15),(0.05,0.2),(0.05,0.3),(0.1,0.3)]):

    '''
    Select a trace, then loop over a range of bands and filter it. Then write to file. Each one of these
    filter files can then be investigated with SHEBA
    '''

    #construct the name of the new SAC file to write

    t = trace.copy()

    fnameparts = basename.split('..')
    stanet = fnameparts[0]
    comptime = fnameparts[1].split('__')

    for band in filters:

        b1 = band[0]
        b2 = band[1]

        newname = stanet+'_'+comptime[1]+'.'+str(b1)+'.'+str(b2)+'.SAC.'+comptime[0]

        t.detrend('demean')
        t.detrend('simple')
        t.filter("bandpass",freqmin=b1,freqmax=b2)

        #first write the file. It will then need to be re-opened as SAC
        #this is bad practice because of the IO, but if the files are coming in in mseed
        #format it may be the only option.

        t.write(newname,format='SAC')

        write_sac(newname, station_code, stlat, stlon, evlat, evlon, evdep, win_start, win_end, cluster_1, cluster_2)

def write_sac(newname,station_code,stlat,stlon,evlat,evlon,evdep,win_start,win_end,cluster_1,cluster_2):

    '''
    Read a SAC file and write important components to its header
    '''

    trace = op.read(newname,format='SAC')
    trace = trace[0]
    #fill the sac header
    trace.stats.sac.kstnm = station_code
    trace.stats.sac.stla = stlat
    trace.stats.sac.stlo = stlon
    trace.stats.sac.evla = evlat
    trace.stats.sac.evlo = evlon
    trace.stats.sac.evdp = evdep

    #set the start and end time of the analysis windows
    trace.stats.sac.user0 = win_start
    trace.stats.sac.user1 = cluster_1
    trace.stats.sac.user2 = cluster_2
    trace.stats.sac.user3 = win_end 

    #Set the azimuth and inc information of the components
    if 'BHZ' in newname:
        trace.stats.sac.cmpaz = 0.0
        trace.stats.sac.cmpinc = 0.0
    elif 'BHE' in newname:
        trace.stats.sac.cmpaz = 90
        trace.stats.sac.cmpinc = 90
    elif 'BHN' in newname:
        trace.stats.sac.cmpaz = 0
        trace.stats.sac.cmpinc = 90

    print("Writing %s" %newname)
    trace.write(newname,format='SAC')


#################################################################################
# MAIN PROGRAM
#################################################################################

def main():
  
    #select the earth model
    earth_model = TauPyModel(model='ak135')

    cwd = os.getcwd()
    waveform_dir = 'sheba_test'
    stations_file = 'Stations_2008-08-01T00:00:00.000000Z_2009-08-01T00:00:00.000000Z_38_-110_40_-112_6.0.dat'
    events_file = 'Events_2008-08-01T00:00:00.000000Z_2009-08-01T00:00:00.000000Z_38_-110_40_-112_6.0.dat'

    #Load station and event files as data frames
    stations = pd.read_csv(stations_file,
                       sep=' ',names=['lon','lat','ele','net','code','stime'])

    events = pd.read_csv(events_file,
                       sep=' ',names=['lon','lat','dep','mag','time'])

    os.chdir(waveform_dir)

    #Make a directory where SAC files for use with sheba can be put
    if not os.path.isdir('sheba_files'):
        os.system('mkdir sheba_files')

    #make a list of all the BHZ mseed files
    all_BHZ_waveforms = glob.glob("*BHZ*Z.mseed") 

    #check for three components
    for wfile in all_BHZ_waveforms[:10]:
        fparts = wfile.split('..')
        three_comp = glob.glob('%s*%s' %(fparts[0],fparts[1][3:]))

        #ensure that three component data exists for this event
        if len(three_comp) == 3:

            trace = op.read(wfile,format='mseed')
            starttime_str = str(trace[0].stats.starttime)[:-8]
            station_code = trace[0].stats.station

            #Get the station and event metadata associated with this waveform

            station_info = stations[stations['code'] == station_code]
            event_info = events[events['time'].str.contains(starttime_str)]

            
            evlat = event_info['lat'].values[0]
            evlon = event_info['lon'].values[0]
            evdep = event_info['dep'].values[0]
            evmag = event_info['mag'].values[0]
            stlat = station_info['lat'].values[0]
            stlon = station_info['lon'].values[0]

            ddeg = locations2degrees(evlat,evlon,stlat,stlon)

            #Get the arrival times. We should add options for SKKS and PKS too!
            arrivals = earth_model.get_travel_times(source_depth_in_km=evdep,distance_in_degree=ddeg,phase_list=["SKS"])
            
            if len(arrivals) >= 1:

                sks_arrival = arrivals[0].time

                #These window definitions will probably need editing
                #We might want to calculate them based on some waveform characteristics
                #such as snr in a sliding window around the sks time
                win_start = 35
                win_end = win_start + 20 
                cluster_1 = win_start + 2
                cluster_2 = win_end - 2

                for mseed_file in three_comp:

                    trace = op.read(mseed_file,format='mseed')

                    #Slice the trace around the SKS arrival
                    trace_slice = trace[0].slice(starttime = (trace[0].stats.starttime + sks_arrival - 20), endtime = (trace[0].stats.starttime + sks_arrival + 80))

                    filter_traces(trace_slice, mseed_file, station_code, stlat, stlon, evlat, evlon, evdep, 
                        win_start, win_end, cluster_1, cluster_2)


    print("Moving SAC files to sheba project directory")
    os.system('mv *SAC* sheba_files')



if __name__ == '__main__':

    main()






