# import matplotlib
# matplotlib.use('Agg')
import obspy
from obspy.geodetics import gps2dist_azimuth
from obspy.core.stream import Stream
from obspy.geodetics import kilometer2degrees

import geographiclib.geodesic as geo

from Get_Parameters import Get_Parameters
from Process_mSEED import Process_mSEED
from Cut_windows import Cut_windows
from MCMC import MCMC

def main():
    # === Read the mSEED file ===
    get_parameters = Get_Parameters()
    mSEED_path = get_parameters.Get_Path()
    stream = obspy.read(mSEED_path)
    # stream.plot()
    #

    # Stream for only BFO station
    # st = Stream()
    # st.append(stream[2])
    # st.append(stream[0])
    # st.append(stream[1])

    # Get the inventory
    # Process = Process_mSEED(stream)

    # Process the data
    # st = Process.remove_instrument()

    PRIOR = get_parameters.PRIOR(stream)
    # PARAMETERS = get_parameters.Get_Unknown()

    la_s=  PRIOR['la_s']
    lo_s=  PRIOR['lo_s']
    la_r = PRIOR['la_r']
    lo_r = PRIOR['lo_r']
    dist, az, baz = gps2dist_azimuth(lat1=la_s, lon1=lo_s,
                                     lat2=la_r, lon2=lo_r,a=PRIOR['radius'], f=PRIOR['f'])
    # PRIOR['baz'] = 243
    # PRIOR['az'] = 84

    # st = Process.automatic_rotate(PRIOR['baz'])

    epi = kilometer2degrees(dist, radius=PRIOR['radius'])
    # epi = 86 # Uncertainty is 10 degrees
    # depth = PARAMETERS['depth_s']
    # depth = 36000 # 25000

    PRIOR = get_parameters.Get_ranges(epi,depth,PRIOR)
    sample_path = get_parameters.Start_sample_path(PRIOR)

    # === Cut the BW windows (P&S) ===
    BW_obs = Cut_windows(PRIOR['VELOC_taup'])

    npts = stream.traces[0].stats.npts
    # BW_obs.Get_bw_windows(stream,epi,depth,PARAMETERS['origin_time'], npts = npts)
    tt_P = obspy.UTCDateTime(2019,1,3,15,9,54.9)
    tt_S = obspy.UTCDateTime(2019,1,3,15,18,34.6)
    BW_obs.Get_bw_windows_MANUAL(stream,tt_P,tt_S,PARAMETERS['origin_time'], npts = npts)

    # === Start the MCMC ===
    mcmc = MCMC(PARAMETERS['origin_time'],PRIOR,sample_path=sample_path)
    mcmc.start_BW(BW_obs)



if __name__ == '__main__':
    main()


