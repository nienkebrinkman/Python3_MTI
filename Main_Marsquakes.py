# import matplotlib
# matplotlib.use('Agg')
import obspy
from obspy.geodetics import gps2dist_azimuth
from obspy.core.stream import Stream
from obspy.geodetics import kilometer2degrees

import geographiclib.geodesic as geo

from Get_Parameters_Marsquake_235 import Get_Parameters
from Process_mSEED import Process_mSEED
from Cut_windows import Cut_windows
from MCMC import MCMC
from Plot_waveforms import Plot_waveforms

def main():
    # === Read the mSEED file ===
    get_parameters = Get_Parameters()
    mSEED_path = get_parameters.Get_Path()
    stream = obspy.read(mSEED_path)
    PRIOR = get_parameters.PRIOR(stream,inventory=False)

    # === Process the data ===
    # Process = Process_mSEED(stream)
    # st = Process.remove_instrument() # Remove instrument response
    # stream = Process.automatic_rotate(PRIOR['baz']) # Rotate your Traces

    # === GET PRIOR ===
    PRIOR = get_parameters.Get_ranges(PRIOR) # WHEN RUN IS NOT YET DONE

    sample_path = get_parameters.Start_sample_path(PRIOR)

    # === Cut the BW windows (P&S) ===
    BW_obs = Cut_windows(PRIOR['VELOC_taup'], P_HP=PRIOR['P_HP'], P_LP=PRIOR['P_LP'], S_HP=PRIOR['S_HP'],
                         S_LP=PRIOR['S_LP'],Pre_P= PRIOR['Pre_P'],Pre_S = PRIOR['Pre_S'],Post_P = PRIOR['Post_P'],Post_S = PRIOR['Post_S'])
    if PRIOR['P_pick'] == None or PRIOR['S_pick'] == None:
        BW_obs.Get_bw_windows(stream, PRIOR['epi_s'], PRIOR['depth_s'], PRIOR['origin_time'], npts=PRIOR['npts'])
    else:
        BW_obs.Get_bw_windows_MANUAL(stream, PRIOR['P_pick'], PRIOR['S_pick'], PRIOR['origin_time'], npts=PRIOR['npts'])

    #
    # tt_P = obspy.UTCDateTime(2019, 5, 23, 2, 22, 59.1) # EVENT 173
    # tt_S = obspy.UTCDateTime(2019, 5, 23, 2, 25, 53.8) #  EVENT 173

    # tt_P = obspy.UTCDateTime(2019, 7, 26, 12, 19, 19.3) # EVENT 235
    # tt_S = obspy.UTCDateTime(2019, 7, 26, 12, 22, 3) #  EVENT 235


    # === Start the MCMC ===
    # mcmc = MCMC(PRIOR['origin_time'], PRIOR, sample_path=sample_path)
    # mcmc.start_BW(BW_obs)

    # === Plot waveforms from a previous run ===
    PRIOR['VELOC'] = PRIOR['VELOC']
    path_txt = '/home/nienke/Documents/Master/Data/Mars/S0235b/waveforms/Output/TAYAK_UPDATE_2.txt'
    savedir = '/home/nienke/Documents/Master/Data/Mars/S0235b/waveforms/Output/'
    skiprows = 40 # 26
    plot = Plot_waveforms(BW_obs,path_txt,savedir,PRIOR,PRIOR['origin_time'],skiprows)
    # plot.get_Cut_waveforms()
    plot.get_waveforms()

if __name__ == '__main__':
    main()


