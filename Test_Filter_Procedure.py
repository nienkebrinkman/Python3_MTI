# import matplotlib
# matplotlib.use('Agg')
import obspy
from obspy.geodetics import gps2dist_azimuth
from obspy.core.stream import Stream
from obspy.geodetics import kilometer2degrees
import toml
import matplotlib.pylab as plt

import geographiclib.geodesic as geo

from Get_Parameters import Get_Parameters
# from Process_mSEED import Process_mSEED
from Cut_windows import Cut_windows
# from MCMC import MCMC
from Plot_waveforms import Plot_waveforms
from Grid_Search import Grid_Search

""" DEFINE YOUR INPUT FILE PATH HERE: """
toml_path = '/home/nienke/Documents/Master/Data/Mars/S0235b/waveforms/Input/GS_4.toml'
# toml_path = '/home/nienke/Documents/Master/Data/Mars/S0235b/waveforms/Input/GS_Trial_1_NO_ZP.toml'


def main():
    # === Read input file ===
    param = toml.load(toml_path, _dict=dict)

    # === Read the mSEED file ===
    get_parameters = Get_Parameters()
    mSEED_path = get_parameters.Get_Path(param['directory'], param['mSEED_file'])
    stream = obspy.read(mSEED_path)
    PRIOR = get_parameters.PRIOR(stream, param['PRIOR'])

    import numpy as np
    len_trace = len(stream.traces[0].data)
    test__data = np.hstack((np.zeros(56279-100),np.hanning(100)* 0.2 ,np.zeros(56280)))
    stream.traces[0].data = test__data
    stream.traces[1].data = test__data
    stream.traces[2].data = test__data

    # === Process the data ===
    # Process = Process_mSEED(stream)
    # st = Process.remove_instrument() # Remove instrument response
    # stream = Process.automatic_rotate(PRIOR['baz']) # Rotate your Traces

    # === GET PRIOR ===
    PRIOR = get_parameters.Get_ranges(PRIOR)  # WHEN RUN IS NOT YET DONE

    sample_path = get_parameters.Start_sample_path(PRIOR)

    # === Cut the BW windows (P&S) ===
    BW_obs = Cut_windows(PRIOR['VELOC_taup'], P_HP=PRIOR['P_HP'], P_LP=PRIOR['P_LP'], S_HP=PRIOR['S_HP'],
                         S_LP=PRIOR['S_LP'], Pre_P=PRIOR['Pre_P'], Pre_S=PRIOR['Pre_S'], Post_P=PRIOR['Post_P'],
                         Post_S=PRIOR['Post_S'], zero_phase=PRIOR['Zero_Phase'], Order=PRIOR['Order'],
                         Taper=PRIOR['Taper_obs'])
    if PRIOR['P_pick'] == None or PRIOR['S_pick'] == None:
        BW_obs.Get_bw_windows(stream, PRIOR['epi_s'], PRIOR['depth_s'], PRIOR['origin_time'])
    else:
        P_pick = obspy.UTCDateTime(2019, 7, 26, 12, 46, 8)
        BW_obs.Get_bw_windows_MANUAL(stream, P_pick, P_pick, PRIOR['origin_time'])

    # === Start Grid Search ===
    depth =  58691.9  #           #5000
    epi =    24.7437  #            #PRIOR['epi_s']
    M0 =  135943762646494.86 #      #PRIOR['M0']
    GS = Grid_Search(PRIOR)
    GS.start_GS(BW_obs,depth,epi,M0)

    # === Start the MCMC (Metropolis Hastings) ===
    # mcmc = MCMC(PRIOR['origin_time'], PRIOR, sample_path=sample_path)
    # mcmc.start_BW(BW_obs)

    # === Plot waveforms from a previous run ===
    # PRIOR['VELOC'] = PRIOR['VELOC']
    # path_txt = '/home/nienke/Documents/Master/Data/Mars/S0235b/waveforms/Output/GS_Trial_1_ZP.txt'
    # savedir = '/home/nienke/Documents/Master/Data/Mars/S0235b/waveforms/Output/'
    # skiprows = 48  # 26
    # plot = Plot_waveforms(BW_obs, path_txt, savedir, PRIOR, skiprows)
    # # plot.get_Cut_waveforms()
    # plot.get_waveforms(Norm=False)


if __name__ == '__main__':
    main()
