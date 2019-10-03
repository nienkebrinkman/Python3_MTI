# import matplotlib
# matplotlib.use('Agg')
import obspy
from obspy.geodetics import gps2dist_azimuth
from obspy.core.stream import Stream
from obspy.geodetics import kilometer2degrees

import geographiclib.geodesic as geo

from Get_Parameters_MSS import Get_Parameters
from Process_mSEED import Process_mSEED
from Cut_windows import Cut_windows
from MCMC import MCMC

def main():
    # === Read the mSEED file ===
    get_parameters = Get_Parameters()
    mSEED_path = get_parameters.Get_Path()
    stream = obspy.read(mSEED_path)
    # stream.plot()

    # Stream for only BFO station
    # st = Stream()
    # st.append(stream[2])
    # st.append(stream[0])
    # st.append(stream[1])

    PRIOR = get_parameters.PRIOR(stream,inventory=False)

    # === Process the data ===
    # Process = Process_mSEED(stream)
    # st = Process.remove_instrument() # Remove instrument response
    # st = Process.automatic_rotate(PRIOR['baz']) # Rotate your Traces

    PRIOR = get_parameters.Get_ranges(PRIOR)
    sample_path = get_parameters.Start_sample_path(PRIOR)

    # === Cut the BW windows (P&S) ===
    BW_obs = Cut_windows(PRIOR['VELOC_taup'])

    npts = stream.traces[0].stats.npts
    # BW_obs.Get_bw_windows(stream, PRIOR['epi_s'], PRIOR['depth_s'], PRIOR['origin_time'], npts = npts)
    tt_P = obspy.UTCDateTime(2019, 1, 3, 15, 9, 54.9)
    tt_S = obspy.UTCDateTime(2019, 1, 3, 15, 18, 34.6)
    BW_obs.Get_bw_windows_MANUAL_OLD(stream, tt_P, tt_S, PRIOR['origin_time'], npts=npts)

    # === Start the MCMC ===
    mcmc = MCMC(PRIOR['origin_time'], PRIOR, sample_path=sample_path)
    mcmc.start_BW(BW_obs)


if __name__ == '__main__':
    main()


