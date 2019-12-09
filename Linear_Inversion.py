"""


"""

# import matplotlib
# matplotlib.use('Agg')
import obspy
import toml
import instaseis

from Get_Parameters import Get_Parameters
from Cut_windows import Cut_windows

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

    # === GET PRIOR ===
    PRIOR = get_parameters.Get_ranges(PRIOR)  # WHEN RUN IS NOT YET DONE

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

    ## === 2. Get Green's Functions ===
    db = instaseis.open_db(PRIOR['VELOC'])
    depth = 58691.9

    a=1

if __name__ == '__main__':
    main()
