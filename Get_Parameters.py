# === Specify parameters ===

import obspy
import numpy as np
import os
from obspy.geodetics import kilometer2degrees
from obspy.geodetics.base import gps2dist_azimuth
from obspy.core.event.event import Event
from obspy.core.event import read_events

from Create_starting_sample import create_starting_sample

class Get_Parameters:
    def Get_Path(self, directory, mSEED_file):
        self.directory = directory
        mSEED_file = mSEED_file
        mSEED_path = self.directory + mSEED_file
        return mSEED_path

    def Start_sample_path(self,PRIOR):
        start_sample_path =  "None"


        if start_sample_path == "None":
            create = create_starting_sample()
            strike = np.random.uniform(PRIOR['strike']['range_min'], PRIOR['strike']['range_max'])
            dip = np.random.uniform(PRIOR['dip']['range_min'], PRIOR['dip']['range_max'])
            rake = np.random.uniform(PRIOR['rake']['range_min'], PRIOR['rake']['range_max'])
            epi = np.random.uniform(PRIOR['epi']['range_min'], PRIOR['epi']['range_max'])
            depth = np.random.uniform(PRIOR['depth']['range_min'], PRIOR['depth']['range_max'])

            start_sample_path = create.get_sample_manual(epi,depth,strike,dip,rake,PRIOR['M0'],self.directory+'start_sample.txt')
        return start_sample_path


    def PRIOR(self, stream, Param):
        if Param['inventory']:
            if Param['inv_file'] == "None":
                raise ValueError('No path for xml file specified!')
            else:
                inv = read_events(Param['inv_file'])

        trace = stream.traces[0]
        PRIOR = {}
        PRIOR['PLOT'] = Param['PLOT']
        PRIOR['save_name'] = Param['save_name']
        PRIOR['save_dir'] = self.directory + 'Output'  #'/home/nienke/MSS'
        if not os.path.exists(PRIOR['save_dir']):
            os.makedirs(PRIOR['save_dir'])

        # = Radius of the body used =
        PRIOR['radius'] = Param['PLANET']['radius']
        # = Flattening planet =
        PRIOR['f']= Param['PLANET']['f']

        # = Receiver =
        if Param['inventory'] == False:
            print('Location of InSight is used!')
            PRIOR['la_r'] = Param['RECEIVER']['la_r']
            PRIOR['lo_r'] = Param['RECEIVER']['lo_r']
        else:
            PRIOR['la_r'] = inv._networks[0].stations[0]._latitude  #4.5 # InSight
            PRIOR['lo_r'] = inv._networks[0].stations[0]._longitude #136 # InSight
        PRIOR['network'] = trace.stats.network
        PRIOR['station'] = trace.stats.channel
        PRIOR['location'] = trace.stats.location
        PRIOR['rec_depth'] = 0#589 # For BFO station

        # = Source =
        ## Catalogue:
        if Param['inventory'] == False:
            PRIOR['origin_time'] = obspy.UTCDateTime(Param['SOURCE']['origin_time'])
            PRIOR['P_pick'] = obspy.UTCDateTime(Param['SOURCE']['P_pick']) # If not know: None
            PRIOR['S_pick'] = obspy.UTCDateTime(Param['SOURCE']['S_pick']) # If not know: None
            PRIOR['depth_s'] = Param['SOURCE']['depth']
            PRIOR['la_s'] = Param['SOURCE']['la_s']
            PRIOR['lo_s'] = Param['SOURCE']['lo_s']
            Mw = Param['SOURCE']['Mw']
        else:
            PRIOR['la_s'] = inv.events[0].origins[0].latitude
            PRIOR['lo_s'] = inv.events[0].origins[0].longitude
            PRIOR['depth_s'] = inv.events[0].origins[0].depth
            PRIOR['origin_time'] = inv.events[0].origins[0].time
            Mw = "None"


        PRIOR['M0'] = self.Magnitude2Scalarmoment(Mw)  # Scalar Moment
        PRIOR['components'] = Param['components']
        PRIOR['kind'] = Param['kind']

        dist, az, baz = gps2dist_azimuth(lat1=PRIOR['la_s'], lon1=PRIOR['lo_s'],
                                         lat2=PRIOR['la_r'], lon2=PRIOR['lo_r'], a=PRIOR['radius'], f=PRIOR['f'])
        PRIOR['baz'] = baz
        PRIOR['az'] = az
        PRIOR['epi_s'] = kilometer2degrees(dist, radius=PRIOR['radius'])

        # = Velocity model =
        PRIOR['VELOC'] = Param['VELOC']['VELOC']
        PRIOR['VELOC_taup'] = Param['VELOC']['VELOC_taup']

        # = Sample information =
        # PRIOR['npts'] = 15000
        PRIOR['Temperature'] = Param['SAMPLE']['Temperature']
        PRIOR['sample_number'] =  Param['SAMPLE']['sample_number']
        if  Param['SAMPLE']['sampling_rate'] == 'None':
            PRIOR['sampling_rate'] = trace.stats.sampling_rate # [Hz] InSight Mission
        else:
            PRIOR['sampling_rate'] = Param['SAMPLE']['sampling_rate']

        PRIOR['directory'] = self.directory

        # = Filter information =
        PRIOR['P_LP'] = 1.0/ Param['FILTER']['P_LP']
        PRIOR['P_HP'] = 1.0/ Param['FILTER']['P_HP']
        PRIOR['S_LP'] = 1.0/ Param['FILTER']['S_LP']
        PRIOR['S_HP'] = 1.0/ Param['FILTER']['S_HP']

        PRIOR['Taper_obs'] = Param['FILTER']['Taper_obs']
        PRIOR['Taper_syn'] = Param['FILTER']['Taper_syn']
        PRIOR['Zero_Phase'] = Param['FILTER']['Zero_Phase']
        PRIOR['Order'] = Param['FILTER']['Order']

        # = Pick information in seconds =
        PRIOR['Pre_P'] = Param['WINDOW']['Pre_P']
        PRIOR['Pre_S'] = Param['WINDOW']['Pre_S']
        PRIOR['Post_P'] = Param['WINDOW']['Post_P']
        PRIOR['Post_S'] = Param['WINDOW']['Post_S']

        return PRIOR

    def Get_ranges(self,PRIOR):
        # = Range orientation angles =
        PRIOR['strike'] = {}
        PRIOR['strike']['range_min'] = 0
        PRIOR['strike']['range_max'] = 359.9
        PRIOR['dip']={}
        PRIOR['dip']['range_min'] = 0
        PRIOR['dip']['range_max'] = 89.9
        PRIOR['angle_spread'] =8
        PRIOR['rake']={}
        PRIOR['rake']['range_min'] = -180
        PRIOR['rake']['range_max'] = 179.9
        PRIOR['rake']['spread'] = 8

        # = Range epi and depth =
        PRIOR['epi']={}
        PRIOR['epi']['range_min'] = PRIOR['epi_s'] - 1
        PRIOR['epi']['range_max'] = PRIOR['epi_s']  + 1
        PRIOR['epi']['spread'] = 1
        if PRIOR['depth_s'] == None:
            PRIOR['depth'] = {}
            PRIOR['depth']['range_min'] = 20000
            PRIOR['depth']['range_max'] = 100000
            PRIOR['depth']['spread'] = 3000
        else:
            PRIOR['depth'] = {}
            PRIOR['depth']['range_min'] = 20000# PRIOR['depth_s'] - 10000
            PRIOR['depth']['range_max'] = 100000# PRIOR['depth_s'] + 10000
            PRIOR['depth']['spread'] = 3000

        return PRIOR

    def Magnitude2Scalarmoment(self,Mw):
        M = 10.0 ** ((Mw / 2.0 * 3.0 + 9.1))
        return M

