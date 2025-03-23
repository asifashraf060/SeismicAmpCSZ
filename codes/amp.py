# Inputs
station    = '30001'
shots      = list(range(67178,70040))

theShot    = '/Users/asifashraf/Documents/SEGY_conversion/GIPP_tools/Project_files/PD16.txt'
theArrival = '/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/Stingray_structures/tlArrival_chk6M35_06-Nov-2024_PPgPn_srModel_06-Nov-2024_CG_SG_softB_5.mat'

## (seconds) This time window determines the range of data after user selects the arrival
timeWindow = .02      

# Import modules
import numpy as np
from obspy.core import read
from obspy.core import UTCDateTime
from obspy.clients.fdsn.client import Client
import matplotlib.pyplot as plt
import os
from obspy import UTCDateTime
import scipy.io as sio
import pandas as pd
from datetime import datetime, timedelta

for i in range(len(shots)):
    try:
        shot = str(shots[i])
        print('shot: '+ shot)
        # load the mseed file for that station
        #st_path = ('/Volumes/LaCie/Ally_Project/' + station + '/' + station + '210620064043.mseed')
        st_path = ('/Volumes/Cascadia_Exp_Backup2/Ally_project/' + station + '/' + station + '210620064043.dp3')
        mseed   = read(st_path)
        tr      = mseed[0]

        # Load the tlArrival structure
        data      = sio.loadmat(theArrival)
        tlArrival = data['tlArrival']
        #Access the structure's fields; use:tlArrival.dtype.names to show the fields
        events   = tlArrival['eventid']
        stations = tlArrival['station']
        times    = tlArrival['time']
        # make arrays for events, stations and times
        events   = events[0][0]
        stations = stations[0][0]
        times    = times[0][0]
        idx_shot    = np.where(events == int(shot))[0]
        idx_station = np.where(stations == (station))[0]
        cmmn        = np.intersect1d(idx_station, idx_shot)
        if len(cmmn) > 0:
                print('Arrival pick for the station-shot exists!')
                print('The station-shot index:' + str(cmmn))
        else:
                print("No station-shot combination found")

        traveltime = times[cmmn]; traveltime = traveltime[0][0]
        print('Picked time: ' + str(traveltime) + ' seconds')

        # Load the shot data
        df          = pd.read_csv(theShot, delim_whitespace=True)
        shotno      = np.array(df.iloc[:,1])                       # All shots
        shottimes   = df.iloc[:,6]                                 # All shot-times
        idx_shotID  = np.where(shotno == int(shot))                # Find the shot index
        shottime    = shottimes[idx_shotID[0][0]]                  # Shottime for the desired shot
        shottime_dt = UTCDateTime(shottime)                        # Convert that shottime to a UTCdatetime format

        arrivalTime = UTCDateTime(shottime_dt + traveltime)                      # Calculate the arrival time for that shot
        before_shot = arrivalTime - 20                              # sec before the shot arrival
        after_shot  = arrivalTime + 20                              # sec after the shot arrival

        tr_slice = tr.slice(starttime=before_shot, endtime=after_shot)

        tmAx = np.linspace(1, 20*2, len(tr_slice.data))

        plt.close('all')
        plt.figure(1)
        plt.plot(tmAx, tr_slice.data)
        plt.savefig("/Users/asifashraf/Documents/code_repo/ampSeismic/plots/shot_" + shot +"_data.png", dpi=300, bbox_inches='tight')
    except:
        pass