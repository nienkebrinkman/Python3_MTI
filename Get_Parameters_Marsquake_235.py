# === Specify parameters ===

import obspy
import numpy as np
import os
from obspy.geodetics import kilometer2degrees
from obspy.geodetics.base import gps2dist_azimuth
from obspy.core.event.event import Event
from obspy.core.event import read_events
from pyrocko import moment_tensor as mtm



from Create_starting_sample import create_starting_sample

class Get_Parameters:
    def Get_Path(self):
        self.directory = '/home/nienke/Documents/Master/Data/Mars/S0235b/waveforms/'# '/home/nienke/MARSQUAKES/'
        # self.directory = '/home/nienke/Documents/Master/Applied_geophysics/Thesis/Data/Mars/S0235b/waveforms/'# '/home/nienke/MARSQUAKES/'
        self.inv = None #
        mSEED_file = 'waveforms_VBB_ZRT.mseed'#'2018-09-05-mww66-hokkaido-japan-region-5.miniseed'
        # mSEED_file = 'waveforms_VBB.mseed'#'2018-09-05-mww66-hokkaido-japan-region-5.miniseed'
        mSEED_path = self.directory + mSEED_file
        return mSEED_path

    def Start_sample_path(self,PRIOR):
        start_sample_path =  None # '/home/nienke/Documents/Master/Data/MSS/start_sample.txt'


        if start_sample_path == None:
            create = create_starting_sample()
            strike = np.random.uniform(PRIOR['strike']['range_min'], PRIOR['strike']['range_max'])
            dip = np.random.uniform(PRIOR['dip']['range_min'], PRIOR['dip']['range_max'])
            rake = np.random.uniform(PRIOR['rake']['range_min'], PRIOR['rake']['range_max'])
            epi = np.random.uniform(PRIOR['epi']['range_min'], PRIOR['epi']['range_max'])
            depth = np.random.uniform(PRIOR['depth']['range_min'], PRIOR['depth']['range_max'])

            start_sample_path = create.get_sample_manual(epi,depth,strike,dip,rake,PRIOR['M0'],self.directory+'start_sample.txt')
        return start_sample_path


    def PRIOR(self,stream,inventory=False):
        if inventory:
            if self.cat == None:
                raise ValueError('No path for xml file specified!')
            else:
                inv = read_events(self.inv)

        trace = stream.traces[0]
        PRIOR = {}
        PRIOR['PLOT'] = False
        PRIOR['save_name'] = 'TAYAK' + trace.id.replace('.','_')
        PRIOR['save_dir'] = self.directory + 'Output'  #'/home/nienke/MSS'
        if not os.path.exists(PRIOR['save_dir']):
            os.makedirs(PRIOR['save_dir'])

        # = Radius of the body used =
        PRIOR['radius'] = 3389.5 # Mars
        # PRIOR['radius'] = 6378137.0 # Earth

        # = Flattening planet =
        PRIOR['f']= 0 # Mars
        # PRIOR['f'] = 1 / 298.257223563 # Earth

        # = Receiver =
        if inventory == False:
            print('Location of InSight is used!')
            PRIOR['la_r'] = 4.5 # InSight
            PRIOR['lo_r'] = 136 # InSight
        else:
            PRIOR['la_r'] = inv._networks[0].stations[0]._latitude  #4.5 # InSight
            PRIOR['lo_r'] = inv._networks[0].stations[0]._longitude #136 # InSight
        PRIOR['network'] = trace.stats.network
        PRIOR['station'] = trace.stats.channel
        PRIOR['location'] = trace.stats.location
        PRIOR['rec_depth'] = 0#589 # For BFO station

        # = Source =
        ## Catalogue:
        if inventory == False:
            PRIOR['origin_time'] = obspy.UTCDateTime(2019, 7, 26, 12, 16, 15)
            PRIOR['P_pick'] = obspy.UTCDateTime(2019, 7, 26, 12, 19, 19.3) # If not know: None
            PRIOR['S_pick'] = obspy.UTCDateTime(2019, 7, 26, 12, 22, 3) # If not know: None
            PRIOR['depth_s'] = 45000
            PRIOR['la_s'] = 10.99
            PRIOR['lo_s'] = 160.95
            Mw = 3.3
        else:
            PRIOR['la_s'] = inv.events[0].origins[0].latitude
            PRIOR['lo_s'] = inv.events[0].origins[0].longitude
            PRIOR['depth_s'] = inv.events[0].origins[0].depth
            PRIOR['origin_time'] = inv.events[0].origins[0].time
            Mw = None


        PRIOR['M0'] = self.Magnitude2Scalarmoment(Mw)  # Scalar Moment
        exp = mtm.magnitude_to_moment(Mw) # Using Pyrocko package....
        # self.M0 = PRIOR['M0']
        PRIOR['components'] = ["Z", "R", "T"]
        PRIOR['kind'] = 'displacement'


        dist, az, baz = gps2dist_azimuth(lat1=PRIOR['la_s'], lon1=PRIOR['lo_s'],
                                         lat2=PRIOR['la_r'], lon2=PRIOR['lo_r'], a=PRIOR['radius'], f=PRIOR['f'])
        PRIOR['baz'] = baz
        PRIOR['az'] = az
        PRIOR['epi_s'] = kilometer2degrees(dist, radius=PRIOR['radius'])

        # PRIOR['baz'] = 243
        # PRIOR['epi_s'] = 86

        # = Velocity model =

        #   -Mars-
        PRIOR['VELOC'] = 'http://instaseis.ethz.ch/blindtest_1s/TAYAK_1s/'
        # PRIOR['VELOC'] = '/opt/databases/Mars/TAYAK_1s'


        PRIOR['VELOC_taup'] = '/home/nienke/Documents/Master/Data/Database/TAYAK.npz'
        # PRIOR['VELOC_taup'] = '/home/nienke/MARSQUAKES/TAYAK.npz'

        #   -Earth-
        # PRIOR['VELOC'] = 'syngine://iasp91_2s'
        # PRIOR['VELOC_taup'] = 'iasp91'

        # = Noise model =
        # PRIOR['noise_model'] = 'Tcompact' #'STS2' #

        # = Sample information =
        PRIOR['npts'] = 30000
        PRIOR['Temperature'] = 1
        PRIOR['sample_number'] = 100000
        # PRIOR['sampling_rate'] = 20 # [Hz]
        PRIOR['sampling_rate'] = trace.stats.sampling_rate # [Hz] InSight Mission

        PRIOR['directory'] = self.directory

        # = Filter information =
        PRIOR['P_LP'] = 1.0/1.0
        PRIOR['P_HP'] = 1.0 / 8.0 # could also be 1.0/10.0
        PRIOR['S_LP'] = 1.0/1.0
        PRIOR['S_HP'] = 1.0 / 10.
        return PRIOR

    def Get_ranges(self,PRIOR):
        # = Range orientation angles =
        PRIOR['strike'] = {}
        PRIOR['strike']['range_min'] = 0
        PRIOR['strike']['range_max'] = 359.9
        PRIOR['dip']={}
        PRIOR['dip']['range_min'] = 0
        PRIOR['dip']['range_max'] = 89.9
        PRIOR['angle_spread'] = 5
        PRIOR['rake']={}
        PRIOR['rake']['range_min'] = -180
        PRIOR['rake']['range_max'] = 179.9
        PRIOR['rake']['spread'] = 5

        # = Range epi and depth =
        PRIOR['epi']={}
        PRIOR['epi']['range_min'] = PRIOR['epi_s'] - 1
        PRIOR['epi']['range_max'] = PRIOR['epi_s']  + 1
        PRIOR['epi']['spread'] = 1
        if PRIOR['depth_s'] == None:
            PRIOR['depth'] = {}
            PRIOR['depth']['range_min'] = 0
            PRIOR['depth']['range_max'] = 50000
            PRIOR['depth']['spread'] = 1000
        else:
            PRIOR['depth'] = {}
            PRIOR['depth']['range_min'] = 20000# PRIOR['depth_s'] - 10000
            PRIOR['depth']['range_max'] = 46000# PRIOR['depth_s'] + 10000
            PRIOR['depth']['spread'] = 1000

        return PRIOR

    def Magnitude2Scalarmoment(self,Mw):
        M = 10.0 ** ((Mw / 2.0 * 3.0 + 9.1))
        return M

