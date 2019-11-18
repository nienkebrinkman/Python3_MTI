import numpy as np
import geographiclib.geodesic as geo
import os
import matplotlib.pylab as plt

from Model_samples import Model_samples
from Cut_windows import Cut_windows
from Get_Seismogram import Get_Seismogram
from Misfit import Misfit


class Grid_Search:
    def __init__(self, PRIOR):
        self.or_time = PRIOR['origin_time']
        self.prior = PRIOR

        self.model = Model_samples(self.prior)
        self.BW_syn = Cut_windows(self.prior['VELOC_taup'], P_HP=PRIOR['P_HP'], P_LP=PRIOR['P_LP'], S_HP=PRIOR['S_HP'],
                                  S_LP=PRIOR['S_LP'], Pre_P=PRIOR['Pre_P'], Pre_S=PRIOR['Pre_S'],
                                  Post_P=PRIOR['Post_P'], Post_S=PRIOR['Post_S'], zero_phase=PRIOR['Zero_Phase'],
                                  Order=PRIOR['Order'], Taper=PRIOR['Taper_syn'])
        self.seis = Get_Seismogram(self.prior)
        self.mis = Misfit()

    def start_GS(self, BW_obs, depth, epi, M0):
        savepath = self.prior['save_dir'] + '/%s.txt' % self.prior['save_name']

        strike_len = np.linspace(0, 360, int(360 / 5) + 1, endpoint=True)
        dip_len = np.linspace(0, 90, int(90 / 5) + 1, endpoint=True)
        rake_len = np.linspace(-180, 180, int(360 / 5) + 1, endpoint=True)


        with open(savepath, 'w') as save_file:
            self.write_par(save_file)
            i = 0
            print('Iteration: %i' % i)
            for i_s, strike in enumerate(strike_len):
                for i_d, dip in enumerate(dip_len):
                    for i_r, rake in enumerate(rake_len):
                        ## Get the synthetic Data:
                        dict = geo.Geodesic(a=self.prior['radius'], f=self.prior['f']).ArcDirect(
                            lat1=self.prior['la_r'],
                            lon1=self.prior['lo_r'],
                            azi1=self.prior['baz'],
                            a12=epi, outmask=1929)
                        st_syn = self.seis.get_seis_manual(la_s=dict['lat2'], lo_s=dict['lon2'], depth=depth,
                                                           strike=strike, dip=dip, rake=rake,
                                                           time=self.or_time, M0=M0)

                        self.BW_syn.Get_bw_windows(st_syn, epi, depth, self.or_time)

                        ## Determine the misfit:
                        Xi_bw, amplitude, time_shift, fig = self.mis.CC_BW(BW_obs, self.BW_syn,
                                                                           self.or_time, self.prior['PLOT'])
                        if self.prior['PLOT'] == True:
                            # self.plot()
                            if not os.path.exists(self.prior['save_dir'] + '/plots/'):
                                os.makedirs(self.prior['save_dir'] + '/plots/')
                            fig.savefig(
                                self.prior['save_dir'] + '/plots/SHIFT_%s_%05i.png' % (self.prior['save_name'], i))
                            plt.close("all")

                        self.write_sample(save_file, epi, depth, strike, dip, rake, M0, Xi_bw, time_shift, i, accept=1)
                        i += 1
                        print('Iteration: %i' % i)
            save_file.close()

    def write_sample(self, file_name, epi, depth, strike, dip, rake, M0, Xi_bw, time_shift, iteration,
                     accept=0):
        s_z = Xi_bw[0]
        s_r = Xi_bw[1]
        s_t = Xi_bw[2]
        p_z = Xi_bw[3]
        p_r = Xi_bw[4]
        Xi = s_z + s_r + s_t + p_z + p_r

        file_name.write("%.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, " % (
            epi, depth, strike, dip, rake, M0, Xi))
        file_name.write("%.4f, %.4f, %.4f, %.4f, %.4f, %.4f, " % (
            p_z, p_r, s_z, s_r, s_t, Xi))
        file_name.write("%i, %i, %i\n\r" % (time_shift[0], time_shift[1], iteration))

    def write_par(self, file_name):

        file_name.write(
            "epi, depth, strike, dip, rake, M0, Total-misfit, p_z, p_r, s_z, s_r, s_t, bw_tot, shift_S, shift-P, Iteration\n\r")
        file_name.write("Velocity Model:%s\n\r" % self.prior['VELOC'])
        file_name.write("Station:%s\n\r" % self.prior['station'])
        file_name.write("Sampling rate:%.2f\n\r" % self.prior['sampling_rate'])
        file_name.write("la_r:%.4f\n\r" % self.prior['la_r'])
        file_name.write("lo_r:%.4f\n\r" % self.prior['lo_r'])
        file_name.write("kind:%s\n\r" % self.prior['kind'])  #
        file_name.write("Pre_P:%.2f\n\r" % self.prior['Pre_P'])  #
        file_name.write("Pre_S:%.2f\n\r" % self.prior['Pre_S'])  #
        file_name.write("Post_P:%.2f\n\r" % self.prior['Post_P'])  #
        file_name.write("Post_S:%.2f\n\r" % self.prior['Post_S'])  #
        file_name.write("P_LP:%.2f\n\r" % self.prior['P_LP'])  #
        file_name.write("P_HP:%.2f\n\r" % self.prior['P_HP'])  #
        file_name.write("S_LP:%.2f\n\r" % self.prior['S_LP'])  #
        file_name.write("S_HP:%.2f\n\r" % self.prior['S_HP'])  #
        file_name.write("Taper_obs:%5i\n\r" % self.prior['Taper_obs'])  #
        file_name.write("Taper_syn:%5i\n\r" % self.prior['Taper_syn'])  #
        file_name.write("Zero_phase:%5i\n\r" % self.prior['Zero_Phase'])  #
        file_name.write("Order:%5i\n\r" % self.prior['Order'])  #
        file_name.write("amount samples:%i\n\r" % self.prior['sample_number'])  #
        file_name.write("Temperature:%i\n\r" % self.prior['Temperature'])  #
        file_name.write("Radius:%.4f\n\r" % self.prior['radius'])  #
        file_name.write("Flattening:%.4f\n\r" % self.prior['f'])  #
        file_name.write("Azimuth:%.4f\n\r" % self.prior['az'])  #
