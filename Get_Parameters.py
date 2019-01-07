# === Specify parameters ===

import obspy
import numpy as np
from obspy.geodetics import kilometer2degrees
from obspy.geodetics.base import gps2dist_azimuth
from obspy.core.event.event import Event

from Create_starting_sample import create_starting_sample

class Get_Parameters:
    def Get_Path(self):
        self.directory ='/home/nienke/Documents/Applied_geophysics/Thesis/BBB_project/Database/MSS/'# '/home/nienke/'#
        mSEED_file = 'mss_event.mseed'#'2018-09-05-mww66-hokkaido-japan-region-5.miniseed'
        mSEED_path = self.directory + mSEED_file
        return mSEED_path

    def Start_sample_path(self,PRIOR):
        start_sample_path = '/home/nienke/Documents/Applied_geophysics/Thesis/anaconda/MSS/start_sample.txt'#'/home/nienke/start_sample.txt'#


        if start_sample_path == None:
            create = create_starting_sample()
            strike = np.random.uniform(PRIOR['strike']['range_min'], PRIOR['strike']['range_max'])
            dip = np.random.uniform(PRIOR['dip']['range_min'], PRIOR['dip']['range_max'])
            rake = np.random.uniform(PRIOR['rake']['range_min'], PRIOR['rake']['range_max'])
            epi = np.random.uniform(PRIOR['epi']['range_min'], PRIOR['epi']['range_max'])
            depth = np.random.uniform(PRIOR['depth']['range_min'], PRIOR['depth']['range_max'])

            start_sample_path = create.get_sample_manual(epi,depth,strike,dip,rake,self.M0,self.directory+'start_sample.txt')
        return start_sample_path


    def PRIOR(self,stream,inventory=None):
        trace = stream.traces[0]
        PRIOR = {}
        PRIOR['PLOT'] = True
        PRIOR['save_name'] = 'MAAK_' + trace.id.replace('.','_')
        PRIOR['save_dir'] = '/home/nienke/Documents/Applied_geophysics/Thesis/anaconda/MSS'#'/home/nienke/MSS'
        # = Radius of the body used =
        PRIOR['radius'] = 3389.5 # Mars
        # PRIOR['radius'] = 6378137.0 # Earth

        # = Flattening planet =
        PRIOR['f']= 0 # Mars
        # PRIOR['f'] = 1 / 298.257223563 # Earth

        # = Receiver =
        if inventory == None:
            print('Location of InSight is used!')
            PRIOR['la_r'] = 4.5 # InSight
            PRIOR['lo_r'] = 136 # InSight
        else:
            PRIOR['la_r'] = inventory._networks[0].stations[0]._latitude  #4.5 # InSight
            PRIOR['lo_r'] = inventory._networks[0].stations[0]._longitude #136 # InSight
        PRIOR['network'] = trace.stats.network
        PRIOR['station'] = trace.stats.channel
        PRIOR['location'] = trace.stats.location
        PRIOR['rec_depth'] = 589 # For BFO station

        # = Source =
        Mw = 3.6
        PRIOR['M0'] = self.Magnitude2Scalarmoment(Mw) # Scalar Moment
        self.M0 = PRIOR['M0']
        PRIOR['components'] = ["Z", "R", "T"]

        PRIOR['kind'] = 'velocity'

        # = Velocity model =

        #   -Mars-
        PRIOR['VELOC'] = 'http://instaseis.ethz.ch/blindtest_1s/MAAK_1s'
        # PRIOR['VELOC'] = 'http://instaseis.ethz.ch/blindtest_1s/MAAK_1s'

        # PRIOR['VELOC'] = '/home/nienke/mnt_databases/databases/blindtestmodels_1s/MAAK_1s'
        # PRIOR['VELOC'] = 'mnt_databases/databases/blindtestmodels_1s/EH45TcoldCrust1'
        # PRIOR['VELOC_taup'] = 'EH45TcoldCrust1b.npz'
        PRIOR['VELOC_taup'] = 'MAAK.npz'

        #   -Earth-
        # PRIOR['VELOC'] = 'syngine://iasp91_2s'
        # PRIOR['VELOC_taup'] = 'iasp91'

        # = Noise model =
        PRIOR['noise_model'] = 'Tcompact' #'STS2' #

        # = Sample information =
        PRIOR['npts'] = 30000
        PRIOR['Temperature'] = 1
        PRIOR['sample_number'] = 5000
        # PRIOR['sampling_rate'] = 20 # [Hz]
        PRIOR['sampling_rate'] = trace.stats.sampling_rate # [Hz] InSight Mission

        PRIOR['directory'] = self.directory
        return PRIOR

    def Get_ranges(self,epi,depth,PRIOR):
        self.epi = epi
        self.depth = depth
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
        PRIOR['epi']['range_min'] = epi - 1
        PRIOR['epi']['range_max'] = epi + 1
        PRIOR['epi']['spread'] = 1
        if self.depth == None:
            PRIOR['depth'] = {}
            PRIOR['depth']['range_min'] = 0
            PRIOR['depth']['range_max'] = 50000
            PRIOR['depth']['spread'] = 1000
        else:
            PRIOR['depth'] = {}
            PRIOR['depth']['range_min'] = depth - 25000
            PRIOR['depth']['range_max'] = depth + 25000
            PRIOR['depth']['spread'] = 1000

        return PRIOR


    def Get_Unknown(self):
        PARAMETERS = {}
        # event = obspy.read_events('/home/nienke/mss_event.xml')

        event = obspy.read_events(
            '/home/nienke/Documents/Applied_geophysics/Thesis/BBB_project/Database/MSS/mss_event.xml')

        magnitude = Event.preferred_magnitude(event.events[0])
        # Mw = magnitude.mag
        source = Event.preferred_origin(event.events[0])
        depth = source.depth
        la_s = source.latitude
        lo_s = source.longitude
        time = source.time

        # Source parameters
        PARAMETERS['la_s'] = -26
        PARAMETERS['lo_s'] = 53

        PARAMETERS['depth_s'] = 36000  # [m]
        PARAMETERS['strike'] = 333  # 79
        PARAMETERS['dip'] = 61  # 50
        PARAMETERS['rake'] = 83  # 20

        PARAMETERS['origin_time'] = obspy.UTCDateTime(2019,1,3,15,00,53)#obspy.UTCDateTime(2018,9,5,18,7,59)
        return PARAMETERS

    def Magnitude2Scalarmoment(self,Mw):
        M = 10.0 ** ((Mw / 2.0 * 3.0 + 9.1))
        return M

