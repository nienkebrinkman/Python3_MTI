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
                                  Order=PRIOR['Order'], Taper=PRIOR['Taper_syn'],Taper_len=PRIOR['Taper_len'], Zero_len=PRIOR['Zero_len'])
        self.seis = Get_Seismogram(self.prior)
        self.mis = Misfit()

    def start_GS(self, BW_obs, depth, epi, M0):
        savepath = self.prior['save_dir'] + '/%s.txt' % self.prior['save_name']
        color = 'b'
        self.plot_original_vs_filter(BW_obs,color,self.prior['save_dir'],"OBS")

        # strike_len = np.linspace(0, 360, int(360 / 5) + 1, endpoint=True)
        # dip_len = np.linspace(0, 90, int(90 / 5) + 1, endpoint=True)
        # rake_len = np.linspace(-180, 180, int(360 / 5) + 1, endpoint=True)


        strike_len = np.array([315,5,155])
        dip_len = np.array([45,75,50])
        rake_len = np.array([60,-95,-150])


        with open(savepath, 'w') as save_file:
            self.write_par(save_file)
            i = 0
            print('Iteration: %i' % i)
            # for i_s, strike in enumerate(strike_len):
            #     for i_d, dip in enumerate(dip_len):
            #         for i_r, rake in enumerate(rake_len):

            for i in range(0,len(dip_len)):
                strike = strike_len[i]
                dip = dip_len[i]
                rake = rake_len[i]
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
                color = 'r'
                self.plot_original_vs_filter(self.BW_syn,color, self.prior['save_dir'], str(i))

                # ## Determine the misfit:
                Xi_bw, amplitude, time_shift, fig = self.mis.CC_BW(BW_obs, self.BW_syn,
                                                                   self.or_time, self.prior['PLOT'])
                if self.prior['PLOT'] == True:
                    # self.plot()
                    if not os.path.exists(self.prior['save_dir'] + '/plots/'):
                        os.makedirs(self.prior['save_dir'] + '/plots/')
                    fig.savefig(
                        self.prior['save_dir'] + '/plots/%.3f_%.3f_%.3f_%05i.png' % (strike,dip,rake, i))
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

    def plot_original_vs_filter(self, stream, color,savepath, save_name):
        stream.original.trim(self.prior['origin_time'])

        plt.figure(figsize=(8, 10))
        ax1 = plt.subplot(321)
        x = np.arange(len(stream.original.traces[0]))
        x_cut = x[stream.or_P_len - 500:stream.or_P_len + stream.P_len + 500]
        plt.plot(stream.original.traces[0], 'g')
        # plt.plot(x,self.BW_syn.original.traces[0], 'r')
        plt.plot(x_cut, stream.P_stream.traces[0], color)
        ymin, ymax = ax1.get_ylim()
        plt.vlines(stream.or_P_len, ymin=ymin, ymax=ymax, colors=color, linewidth=3, label='Obs_P')
        plt.vlines(stream.or_P_len + self.prior['Pre_P'] * 20 + self.prior['Post_P'] * 20, ymin=ymin, ymax=ymax,
                   colors=color, linewidth=3, label='Obs_P')
        plt.xlim(stream.or_P_len - 100, stream.or_P_len + 200)
        min_ = min(np.hstack((stream.original.traces[0].data[stream.or_P_len - 100:stream.or_P_len + 400],
                              stream.P_stream.traces[0].data)))
        max_ = max(np.hstack((stream.original.traces[0].data[stream.or_P_len - 100:stream.or_P_len + 400],
                              stream.P_stream.traces[0].data)))
        plt.ylim(min_, max_)
        plt.grid(True)
        plt.title('PZ')

        ax2 = plt.subplot(323)
        x = np.arange(len(stream.original.traces[1]))
        x_cut = x[stream.or_P_len - 500:stream.or_P_len + stream.P_len + 500]
        plt.plot(stream.original.traces[1], 'g')
        # plt.plot(x,self.BW_syn.original.traces[0], 'r')
        plt.plot(x_cut, stream.P_stream.traces[1], color)
        ymin, ymax = ax2.get_ylim()
        plt.vlines(stream.or_P_len, ymin=ymin, ymax=ymax, colors=color, linewidth=3, label='Obs_P')
        plt.vlines(stream.or_P_len + self.prior['Pre_P'] * 20 + self.prior['Post_P'] * 20, ymin=ymin, ymax=ymax,
                   colors=color, linewidth=3, label='Obs_P')
        plt.xlim(stream.or_P_len - 100, stream.or_P_len + 200)
        min_ = min(np.hstack((stream.original.traces[1].data[stream.or_P_len - 100:stream.or_P_len + 400],
                              stream.P_stream.traces[1].data)))
        max_ = max(np.hstack((stream.original.traces[1].data[stream.or_P_len - 100:stream.or_P_len + 400],
                              stream.P_stream.traces[1].data)))
        plt.ylim(min_, max_)
        plt.grid(True)
        plt.title('PR')

        ax3 = plt.subplot(322)
        x = np.arange(len(stream.original.traces[0]))
        x_cut = x[stream.or_S_len - 500:stream.or_S_len + stream.S_len + 500]
        plt.plot(stream.original.traces[0],'g')
        # plt.plot(x,self.BW_syn.original.traces[0], 'r')
        plt.plot(x_cut, stream.S_stream.traces[0], color)
        ymin, ymax = ax3.get_ylim()
        plt.vlines(stream.or_S_len, ymin=ymin, ymax=ymax, colors=color, linewidth=3, label='Obs_S')
        plt.vlines(stream.or_S_len + self.prior['Pre_S'] * 20 + self.prior['Post_S'] * 20, ymin=ymin, ymax=ymax,
                   colors=color, linewidth=3, label='Obs_S')
        plt.xlim(stream.or_S_len - 100, stream.or_S_len + 400)
        min_ = min(np.hstack((stream.original.traces[0].data[stream.or_S_len - 100:stream.or_S_len + 400],
                              stream.S_stream.traces[0].data)))
        max_ = max(np.hstack((stream.original.traces[0].data[stream.or_S_len - 100:stream.or_S_len + 400],
                              stream.S_stream.traces[0].data)))
        plt.ylim(min_, max_)
        plt.grid(True)
        plt.title('SZ')

        ax4 = plt.subplot(324)
        x = np.arange(len(stream.original.traces[1]))
        x_cut = x[stream.or_S_len - 500:stream.or_S_len + stream.S_len + 500]
        plt.plot(stream.original.traces[1], 'g')
        # plt.plot(x,self.BW_syn.original.traces[0], 'r')
        plt.plot(x_cut, stream.S_stream.traces[1], color)
        ymin, ymax = ax4.get_ylim()
        plt.vlines(stream.or_S_len, ymin=ymin, ymax=ymax, colors=color, linewidth=3, label='Obs_S')
        plt.vlines(stream.or_S_len + self.prior['Pre_S'] * 20 + self.prior['Post_S'] * 20, ymin=ymin, ymax=ymax,
                   colors=color, linewidth=3, label='Obs_S')
        plt.xlim(stream.or_S_len - 100, stream.or_S_len + 400)
        min_ = min(np.hstack((stream.original.traces[1].data[stream.or_S_len - 100:stream.or_S_len + 400],stream.S_stream.traces[1].data)))
        max_ = max(np.hstack((stream.original.traces[1].data[stream.or_S_len - 100:stream.or_S_len + 400],stream.S_stream.traces[1].data)))
        plt.ylim(min_,max_)
        plt.grid(True)
        plt.title('SR')

        ax5 = plt.subplot(326)
        x = np.arange(len(stream.original.traces[2]))
        x_cut = x[stream.or_S_len - 500:stream.or_S_len + stream.S_len + 500]
        plt.plot(stream.original.traces[2], 'g')
        # plt.plot(x,self.BW_syn.original.traces[0], 'r')
        plt.plot(x_cut, stream.S_stream.traces[2], color)
        ymin, ymax = ax5.get_ylim()
        plt.vlines(stream.or_S_len, ymin=ymin, ymax=ymax, colors=color, linewidth=3, label='Obs_S')
        plt.vlines(stream.or_S_len + self.prior['Pre_S'] * 20 + self.prior['Post_S'] * 20, ymin=ymin, ymax=ymax,
                   colors=color, linewidth=3, label='Obs_S')
        plt.xlim(stream.or_S_len - 100, stream.or_S_len + 400)
        min_ = min(np.hstack((stream.original.traces[2].data[stream.or_S_len - 100:stream.or_S_len + 400],stream.S_stream.traces[2].data)))
        max_ = max(np.hstack((stream.original.traces[2].data[stream.or_S_len - 100:stream.or_S_len + 400],stream.S_stream.traces[2].data)))
        plt.ylim(min_,max_)
        plt.grid(True)
        plt.title('ST')


        plt.tight_layout()
        savepath = savepath + '/plots/'
        if not os.path.exists(savepath):
            os.makedirs(savepath)
        plt.savefig(savepath + '/Waveforms_%s.png' % save_name)
        plt.close("all")


