# import matplotlib
# matplotlib.use('Agg')
import obspy
import toml

from Get_Parameters import Get_Parameters
# from Process_mSEED import Process_mSEED
from Cut_windows import Cut_windows
# from MCMC import MCMC
from Plot_waveforms import Plot_waveforms
from Grid_Search import Grid_Search

""" DEFINE YOUR INPUT FILE PATH HERE: """
toml_path = '/home/nienke/Documents/Master/Data/Mars/S0235b/waveforms/Input/TAYAK.toml'


# toml_path = '/home/nienke/Documents/Master/Data/Mars/S0173a/waveforms/Input/EH45Tcold_173.toml'


def main():
    # === Read input file ===
    param = toml.load(toml_path, _dict=dict)

    # === Read the mSEED file ===
    get_parameters = Get_Parameters()
    mSEED_path = get_parameters.Get_Path(param['directory'], param['mSEED_file'])
    stream = obspy.read(mSEED_path)
    PRIOR = get_parameters.PRIOR(stream, param['PRIOR'])

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
                         Post_S=PRIOR['Post_S'], global_P_shift=PRIOR['Global_P_shift'],
                         global_S_shift=PRIOR['Global_S_shift'], zero_phase=PRIOR['Zero_Phase'], Order=PRIOR['Order'],
                         Taper=PRIOR['Taper_obs'], Taper_len=PRIOR['Taper_len'], Zero_len=PRIOR['Zero_len'])
    if PRIOR['P_pick'] == None or PRIOR['S_pick'] == None:
        BW_obs.Get_bw_windows(stream, PRIOR['epi_s'], PRIOR['depth_s'], PRIOR['origin_time'], MANUAL=False)
    else:
        BW_obs.Get_bw_windows(stream, PRIOR['P_pick'], PRIOR['S_pick'], PRIOR['origin_time'], MANUAL=True)

    # === Start Grid Search ===
    # depth = 58691.9##           #5000
    # epi =    PRIOR['epi_s']#24.7437  #            #PRIOR['epi_s']
    # M0 = PRIOR['M0']# 135943762646494.86 #      #
    # GS = Grid_Search(PRIOR)
    # GS.start_GS(BW_obs,depth,epi,M0)

    # === Start the MCMC (Metropolis Hastings) ===
    # mcmc = MCMC(PRIOR['origin_time'], PRIOR, sample_path=sample_path)
    # mcmc.start_BW(BW_obs)

    # === Plot waveforms from a previous run ===
    PRIOR['VELOC'] = PRIOR['VELOC']
    path_txt = '/home/nienke/Documents/Master/Data/Mars/S0235b/waveforms/Output/TAYAK_Shift_1.txt'
    savedir = '/home/nienke/Documents/Master/Data/Mars/S0235b/waveforms/Output/'
    # path_txt = '/home/nienke/Documents/Master/Data/Mars/S0173a/waveforms/Output/TAYA.txt'
    # savedir = '/home/nienke/Documents/Master/Data/Mars/S0173a/waveforms/Output/'
    skiprows = 56
    plot = Plot_waveforms(BW_obs, path_txt, savedir, PRIOR, skiprows)
    # plot.get_Cut_waveforms()
    plot.get_waveforms(Norm=True)


if __name__ == '__main__':
    main()
