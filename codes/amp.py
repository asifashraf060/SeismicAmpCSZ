## Inputs
station      = '30001'
shots        = list(range(67178,67199))

shot_file    = 'PD16.txt'
arrival_file = 'tlArrival_chk6M35_06-Nov-2024_PPgPn_srModel_06-Nov-2024_CG_SG_softB_5.mat'

# shotfile
import os
script_dir  = os.path.dirname(os.path.abspath(__file__))
parent_dir  = os.path.dirname(script_dir)
inp_dir     = os.path.join(parent_dir, 'input')
out_dir     = os.path.join(parent_dir, 'output')

theShot     = inp_dir + '/' + shot_file
theArrival  = inp_dir + '/' + arrival_file
thePlot     = out_dir + '/plots/' + station
theData     = out_dir + '/data/' + station

# make those folders, if not in there
os.makedirs(thePlot, exist_ok=True)
os.makedirs(theData, exist_ok=True)
  
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

# make the output ascii file
out_txt = os.path.join(theData, "shot_rms.txt")
with open(out_txt, "w", encoding="utf-8") as f:

    f.write("shot\tRMS-amp\n")

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
                    print('     Arrival pick for the station-shot exists!')
                    print('     The station-shot index:' + str(cmmn))
            else:
                    print("     No station-shot combination found")

            traveltime = times[cmmn]; traveltime = traveltime[0][0]
            print('         Picked time: ' + str(traveltime) + ' seconds')

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

            tr_slice  = tr.slice(starttime=before_shot, endtime=after_shot)
            tr_slice2 = tr.slice(starttime=arrivalTime, endtime=arrivalTime + 1)

            # Apply a Butterworth bandpass filter to the trace
            tr_slice.filter("bandpass", freqmin=3.0, freqmax=12.0, corners=4, zerophase=True)
            tr_slice2.filter("bandpass", freqmin=3.0, freqmax=12.0, corners=4, zerophase=True)

            tmAx  = np.linspace(-20, 20, len(tr_slice.data))
            tmAx2 = np.linspace(0, 1, len(tr_slice2.data))

            # find where the maximum amplitude occurs
            idx_maxAmp = np.where(tr_slice2.data == np.max(tr_slice2.data))
            idx_maxAmp = idx_maxAmp[0][0]
            tm_maxAmp  = tmAx2[idx_maxAmp]
            tr_slice3  = tr.slice(starttime=(arrivalTime+tm_maxAmp - .15), endtime=(arrivalTime+tm_maxAmp + .15))

            rms_amp    = np.sqrt(np.mean(abs(tr_slice3.data)))

            f.write(f"{int(shot)}\t{rms_amp}\n")

        
            plt.close('all')
            plt.figure(figsize=(10, 10))

            # First subplot: with y-axis limits
            plt.subplot(3, 1, 1)   # (nrows, ncols, index)
            plt.plot(tmAx, tr_slice.data)
            plt.title(f'Filtered Y-axis Data (Shot {shot})')
            plt.xlabel('Time relative to arrival [s]')
            plt.ylabel('Amplitude')
            plt.ylim([-0.015, 0.015])
            plt.axvline(0, color='r', linestyle='--', linewidth=1, label='Arrival Time')
            plt.legend()
            plt.grid()

            # Second subplot: without y-axis limits
            plt.subplot(3, 1, 2)
            plt.plot(tmAx, tr_slice.data)
            plt.title(f'Filtered Y-axis Data (Shot {shot}) - Zoomed In')
            plt.xlabel('Time relative to arrival [s]')
            plt.xlim(-1,3)
            plt.ylabel('Amplitude')
            plt.axvline(0, color='r', linestyle='--', linewidth=1, label='Arrival Time')
            plt.legend()
            plt.grid()

            plt.subplot(3, 1, 3)   # (nrows, ncols, index)
            plt.plot(tmAx2, tr_slice2.data)
            plt.title(f'Filtered Y-axis Data (Shot {shot}) - Further zoomed in')
            plt.xlabel('Time relative to arrival [s]')
            plt.ylabel('Amplitude')
            plt.axvline(0, color='r', linestyle='--', linewidth=1, label='Arrival Time')
            plt.axvline(tmAx2[idx_maxAmp], color='g', linestyle='--', linewidth=1, label='Max Amplitude')
            plt.plot([0+tm_maxAmp - .15, 0+tm_maxAmp + .15], [0, 0], '-g')
            plt.legend()
            plt.grid()

            plt.tight_layout()
            plt.savefig(thePlot + '/shot' + shot + ".png", dpi=300, bbox_inches='tight')

        except:
            pass