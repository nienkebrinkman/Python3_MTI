import os
import numpy as np
import pandas as pd
import matplotlib.pylab as plt
import geographiclib.geodesic as geo
from obspy.imaging.beachball import aux_plane
import pylab
from obspy.imaging.beachball import beachball
from matplotlib.patches import Circle
from obspy.imaging.beachball import beach
import matplotlib.image as mpimg
import io
import matplotlib.gridspec as gridspec
import obspy
from obspy.core.stream import Trace,Stream

from Get_Seismogram import Get_Seismogram
from Cut_windows import Cut_windows

class Plot_waveforms:
    def __init__(self,BW_obs,path_txt_inversion,save_directory,PRIOR,otime,skiprows=26):
        self.dir = save_directory
        if not os.path.exists(self.dir):
            os.makedirs(self.dir)


        self.column_names = ["Epi", "Depth", "Strike", "Dip", "Rake","M0","Total_misfit","p_z","p_r","s_z","s_r","s_t",'bw_tot','Shift_S','Shift_P','accept']
        data = np.loadtxt(path_txt_inversion, delimiter=',', skiprows=skiprows)
        self.df = pd.DataFrame(data,columns=self.column_names)
        self.BW_obs = BW_obs
        self.prior = PRIOR
        self.otime = otime

    def get_waveforms(self):
        # fig_bb, ax_bb = plt.subplots(1, 1, figsize=(4, 4))
        #
        # ax_bb.set_xticks([])
        # ax_bb.set_yticks([])
        # ax_bb.axis('off')


        epi = self.df['Epi']
        depth = self.df['Depth']
        strike = self.df['Strike']
        dip = self.df['Dip']
        rake = self.df['Rake']
        M0 = self.df['M0']

        P_shift = self.df['Shift_P']
        S_shift = self.df['Shift_S']


        seis = Get_Seismogram(self.prior)
        BW_syn = Cut_windows(self.prior['VELOC_taup'],P_HP = self.prior['P_HP'], P_LP= self.prior['P_LP'], S_HP =  self.prior['S_HP'], S_LP= self.prior['S_LP'], Pre_P=self.prior['Pre_P'],
                              Pre_S=self.prior['Pre_S'], Post_P=self.prior['Post_P'], Post_S=self.prior['Post_S'])


        fig = plt.figure(figsize=(10, 10))
        delta = self.BW_obs.P_stream.traces[0].meta.delta
        p_time_array = np.arange(len(self.BW_obs.P_stream.traces[0].data)) * delta
        s_time_array = np.arange(len(self.BW_obs.S_stream.traces[0].data)) * delta

        start_P = int((self.BW_obs.start_P.timestamp - self.otime.timestamp - 50) / delta)
        end_P = int((self.BW_obs.start_P.timestamp - self.otime.timestamp + 130) / delta)

        start_S = int((self.BW_obs.start_S.timestamp - self.otime.timestamp - 50) / delta)
        end_S = int((self.BW_obs.start_S.timestamp - self.otime.timestamp + 100) / delta)

        ax1 = plt.subplot2grid((5, 1), (0, 0))
        ax2 = plt.subplot2grid((5, 1), (1, 0))
        ax3 = plt.subplot2grid((5, 1), (2, 0))
        ax4 = plt.subplot2grid((5, 1), (3, 0))
        ax5 = plt.subplot2grid((5, 1), (4, 0))

        n_lowest = 30
        lowest_indices = self.df['Total_misfit'].values.argsort()[0:n_lowest]
        lowest_misfits = self.df['Total_misfit'].values[lowest_indices]
        depths_inds = self.df['Depth'].values.argsort()[-n_lowest:]
        depths = self.df['Depth'].values[depths_inds]


        for i in np.arange(len(epi) - 100, len(epi), 1):
        # for i,v in enumerate(lowest_indices):
        # for i in np.arange(100):
            dict = geo.Geodesic(a=self.prior['radius'], f=self.prior['f']).ArcDirect(lat1=self.prior['la_r'], lon1=self.prior['lo_r'],
                                                                           azi1=self.prior['baz'], a12=epi[i],
                                                                           outmask=1929)

            st_syn = seis.get_seis_manual(la_s=dict['lat2'], lo_s=dict['lon2'], depth=depth[i],
                                          strike=strike[i], dip=dip[i], rake=rake[i],
                                          time=self.otime, M0=M0[i])

            BW_syn.Get_bw_windows(st_syn, epi[i], depth[i], self.otime, self.prior['npts'])

            P_shift_array = self.shift(BW_syn.P_stream.traces[0].data, int(P_shift[i]))
            #
            # ax1.plot(p_time_array[start_P:end_P], BW_syn.P_stream.traces[0].data[start_P:end_P], 'g',
            #          label='Synthetic', linewidth = 0.1)

            ax1.plot(p_time_array[start_P:end_P], self.normalize(P_shift_array[start_P:end_P]), 'r',
                     label='Synthetic', linewidth=0.1)
            # ax1.plot( P_shift_array, 'r',  label='Synthetic', linewidth = 0.1)
            ax1.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
            ax1.tick_params(axis='x', labelsize=18)
            ax1.tick_params(axis='y', labelsize=18)

            # plt.tight_layout()
            # plt.legend(loc='lower left', fontsize=15)

            P_shift_array = self.shift(BW_syn.P_stream.traces[1].data, int(P_shift[i]))
            ax2.plot(p_time_array[start_P:end_P], self.normalize(P_shift_array[start_P:end_P]), 'r',
                     label='Synthetic', linewidth=0.1)
            # ax2.plot(P_shift_array, 'r',label='Synthetic', linewidth = 0.1)
            ax2.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
            ax2.tick_params(axis='x', labelsize=18)
            ax2.tick_params(axis='y', labelsize=18)
            # plt.tight_layout()

            S_shift_array = self.shift(BW_syn.S_stream.traces[0].data, int(S_shift[i]))
            ax3.plot(s_time_array[start_S:end_S], self.normalize(S_shift_array[start_S:end_S]), 'r', linewidth=0.1)
            # ax3.plot(S_shift_array, 'r', linewidth = 0.1)
            ax3.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
            ax3.tick_params(axis='x', labelsize=18)
            ax3.tick_params(axis='y', labelsize=18)
            # plt.tight_layout()

            S_shift_array = self.shift(BW_syn.S_stream.traces[1].data, int(S_shift[i]))
            ax4.plot(s_time_array[start_S:end_S], self.normalize(S_shift_array[start_S:end_S]), 'r', linewidth=0.1)
            # ax4.plot(S_shift_array, 'r', linewidth = 0.1)
            ax4.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
            ax4.tick_params(axis='x', labelsize=18)
            ax4.tick_params(axis='y', labelsize=18)
            # plt.tight_layout()

            S_shift_array = self.shift(BW_syn.S_stream.traces[2].data, int(S_shift[i]))
            ax5.plot(s_time_array[start_S:end_S], self.normalize(S_shift_array[start_S:end_S]), 'r', linewidth=0.1)
            # ax5.plot( S_shift_array, 'r', linewidth = 0.1)
            ax5.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
            ax5.tick_params(axis='x', labelsize=18)
            ax5.tick_params(axis='y', labelsize=18)
            ax5.set_xlabel(self.BW_obs.start_P.strftime('From P arrival: %Y-%m-%dT%H:%M:%S + [sec]'), fontsize=18)

        ax1.plot(p_time_array[start_P:end_P], self.normalize(self.BW_obs.P_stream.traces[0].data[start_P:end_P]), 'b')
        # ax1.plot( BW_obs.P_stream.traces[0].data, 'b', label='Observed', linewidth = 0.1)
        ymin, ymax = ax1.get_ylim()
        xmin, xmax = ax1.get_xlim()
        ax1.text(xmax - 5, ymax / 1.7, "P-Z", fontsize=20, color='b')

        ax2.plot(p_time_array[start_P:end_P], self.normalize(self.BW_obs.P_stream.traces[1].data[start_P:end_P]), 'b')
        # ax2.plot(BW_obs.P_stream.traces[1].data, 'b', linewidth = 0.1)
        ymin, ymax = ax2.get_ylim()
        xmin, xmax = ax2.get_xlim()
        ax2.text(xmax - 5, ymax / 1.7, "P-R", fontsize=20, color='b')

        ax3.plot(s_time_array[start_S:end_S], self.normalize(self.BW_obs.S_stream.traces[0].data[start_S:end_S]), 'b')
        # ax3.plot( BW_obs.S_stream.traces[0].data, 'b', linewidth = 0.1)
        ymin, ymax = ax3.get_ylim()
        xmin, xmax = ax3.get_xlim()
        ax3.text(xmax - 10, ymax / 1.5, "S-Z", fontsize=20, color='b')

        ax4.plot(s_time_array[start_S:end_S], self.normalize(self.BW_obs.S_stream.traces[1].data[start_S:end_S]), 'b')
        # ax4.plot(BW_obs.S_stream.traces[1].data, 'b', linewidth = 0.1)
        ymin, ymax = ax4.get_ylim()
        xmin, xmax = ax4.get_xlim()
        ax4.text(xmax - 10, ymax / 1.5, "S-R", fontsize=20, color='b')

        ax5.plot(s_time_array[start_S:end_S], self.normalize(self.BW_obs.S_stream.traces[2].data[start_S:end_S]), 'b')
        # ax5.plot(BW_obs.S_stream.traces[2].data, 'b', linewidth = 0.1)
        ymin, ymax = ax5.get_ylim()
        xmin, xmax = ax5.get_xlim()
        ax5.text(xmax - 10, ymax / 1.7, "S-T", fontsize=20, color='b')

        plt.tight_layout()

        # plt.show()
        plt.savefig(self.dir + '/Waveforms_Normalized.pdf')
        # plt.show()
        plt.close()

    def get_Cut_waveforms(self):
        epi = self.df['Epi']
        depth = self.df['Depth']
        strike = self.df['Strike']
        dip = self.df['Dip']
        rake = self.df['Rake']
        M0 = self.df['M0']

        P_shift = self.df['Shift_P']
        S_shift = self.df['Shift_S']


        seis = Get_Seismogram(self.prior)

        stream = self.BW_obs.original
        BW_P_stream,BW_S_stream , BW_start_P, BW_start_S= self.cut_stream(stream, Manual= True)


        fig = plt.figure(figsize=(10, 10))
        delta = BW_P_stream.traces[0].meta.delta
        p_time_array = np.arange(len(BW_P_stream.traces[0].data)) * delta
        s_time_array = np.arange(len(BW_S_stream.traces[0].data)) * delta

        # start_P = int((self.BW_obs.start_P.timestamp - self.otime.timestamp - 50) / delta)
        # end_P = int((self.BW_obs.start_P.timestamp - self.otime.timestamp + 130) / delta)
        #
        # start_S = int((self.BW_obs.start_S.timestamp - self.otime.timestamp - 50) / delta)
        # end_S = int((self.BW_obs.start_S.timestamp - self.otime.timestamp + 100) / delta)

        ax1 = plt.subplot2grid((5, 1), (0, 0))
        ax2 = plt.subplot2grid((5, 1), (1, 0))
        ax3 = plt.subplot2grid((5, 1), (2, 0))
        ax4 = plt.subplot2grid((5, 1), (3, 0))
        ax5 = plt.subplot2grid((5, 1), (4, 0))

        for i in np.arange(len(epi) - 2, len(epi), 1):
        # for i in np.arange(100):
            dict = geo.Geodesic(a=self.prior['radius'], f=self.prior['f']).ArcDirect(lat1=self.prior['la_r'], lon1=self.prior['lo_r'],
                                                                           azi1=self.prior['baz'], a12=epi[i],
                                                                           outmask=1929)

            st_syn = seis.get_seis_manual(la_s=dict['lat2'], lo_s=dict['lon2'], depth=depth[i],
                                          strike=strike[i], dip=dip[i], rake=rake[i],
                                          time=self.otime, M0=M0[i])

            Syn_P_stream, Syn_S_stream ,  Syn_start_P, Syn_start_S= self.cut_stream(st_syn, epi = epi[i], depth = depth[i])
            # P_shift_array = self.shift(Syn_P_stream.traces[0].data, int(P_shift[i]))

            ax1.plot(p_time_array, self.normalize(Syn_P_stream.traces[0].data), 'r', alpha = 0.2)
            ax2.plot(p_time_array, self.normalize(Syn_P_stream.traces[1].data), 'r', alpha = 0.2)
            ax3.plot(s_time_array, self.normalize(Syn_S_stream.traces[0].data), 'r', alpha = 0.2)
            ax4.plot(s_time_array, self.normalize(Syn_S_stream.traces[1].data), 'r', alpha = 0.2)
            ax5.plot(s_time_array, self.normalize(Syn_S_stream.traces[2].data), 'r', alpha = 0.2)


        ## PLOT OBSERVED
        ax1.plot(p_time_array, self.normalize(BW_P_stream.traces[0].data), 'b')
        ax2.plot(p_time_array, self.normalize(BW_P_stream.traces[1].data), 'b')
        ax3.plot(s_time_array, self.normalize(BW_S_stream.traces[0].data), 'b')
        ax4.plot(s_time_array, self.normalize(BW_S_stream.traces[1].data), 'b')
        ax5.plot(s_time_array, self.normalize(BW_S_stream.traces[2].data), 'b')


        plt.tight_layout()

        plt.show()
        # plt.savefig(self.dir + '/Waveforms_Normalized.pdf')
        # plt.show()
        plt.close()

        stream.filter('highpass', freq = self.prior['S_HP'] , zerophase = True)
        st_syn.filter('highpass', freq = self.prior['S_HP'] , zerophase = True)


        stream.filter('lowpass', freq = self.prior['S_LP'] , zerophase = True)
        st_syn.filter('lowpass', freq = self.prior['S_LP'] , zerophase = True)

        stream.trim(self.prior['origin_time'])
        st_syn.trim(self.prior['origin_time'])

        ax1 = plt.subplot(111)
        plt.plot(stream.traces[2].data, 'b')
        plt.plot(st_syn.traces[2].data, 'r')
        ymin, ymax = ax1.get_ylim()
        plt.vlines(BW_start_P, ymin=ymin, ymax=ymax, colors='b', linewidth=3, label='Obs_P')
        plt.vlines(Syn_start_P, ymin=ymin, ymax=ymax, colors='r', linewidth=3, label = 'Syn_P %.2f' %depth[i])
        plt.vlines(BW_start_S , ymin=ymin, ymax=ymax, colors='b', linewidth=3, label='Obs_S')
        plt.vlines(Syn_start_S, ymin=ymin, ymax=ymax, colors='r', linewidth=3, label = 'Syn_S %.2f' %depth[i])
        plt.legend()
        plt.show()


    def cut_stream(self,stream, Manual = False, epi = None, depth = None):
        WINDOWS = Cut_windows(self.prior['VELOC_taup'], P_HP=self.prior['P_HP'], P_LP=self.prior['P_LP'],
                              S_HP=self.prior['S_HP'], S_LP=self.prior['S_LP'], Pre_P=self.prior['Pre_P'],
                              Pre_S=self.prior['Pre_S'], Post_P=self.prior['Post_P'], Post_S=self.prior['Post_S'])

        or_time_sec = self.otime.timestamp
        if Manual:
            start_P = obspy.UTCDateTime(self.prior['P_pick'].timestamp - self.prior['Pre_P'])
            or_P_len = int((start_P - self.otime) / stream.traces[0].stats.delta)
            start_S = obspy.UTCDateTime(self.prior['S_pick'].timestamp - self.prior['Pre_S'])
            or_S_len = int((start_S - self.otime) / stream.traces[0].stats.delta)

            self.dt = stream.traces[0].stats.delta

            end_P = obspy.UTCDateTime(self.prior['P_pick'].timestamp + self.prior['Post_P'])
            end_S = obspy.UTCDateTime(self.prior['S_pick'].timestamp + self.prior['Post_S'])
        else:
            tt_P = WINDOWS.get_P(epi, depth)
            tt_S = WINDOWS.get_S(epi, depth)

            start_P = obspy.UTCDateTime(or_time_sec + tt_P - self.prior['Pre_P'])
            or_P_len = int((start_P - self.otime)/ stream.traces[0].stats.delta)
            start_S = obspy.UTCDateTime(or_time_sec + tt_S - self.prior['Pre_S'])
            or_S_len = int((start_S - self.otime) / stream.traces[0].stats.delta)


            end_P = obspy.UTCDateTime(or_time_sec + tt_P + self.prior['Post_P'])
            end_S = obspy.UTCDateTime(or_time_sec + tt_S + self.prior['Post_S'])


        P_stream = Stream()
        S_stream = Stream()
        # BW_stream = Stream()

        for i,trace in enumerate(stream.traces):
            # Filter entire trace for P- and S separately.
            trace.filter('highpass', freq=1. / (end_P - start_P), zerophase=True)
            trace.filter('highpass', freq=1. / (end_S - start_S), zerophase=True)
            if i is not 2:
                P_trace = Trace.slice(trace, start_P, end_P)
                P_trace.taper(0.05, 'hann') # Taper the data
                P_len = len(P_trace)
                npts_p = P_len + 2 * or_P_len
                P_stream += P_trace
                P_stream = WINDOWS.Filter(P_stream,self.prior['P_HP'],self.prior['P_LP'])
            S_trace = Trace.slice(trace, start_S, end_S)
            S_trace.taper(0.05, 'hann')
            S_len = len(S_trace)
            npts_s = S_len + 2 * or_S_len
            S_stream += S_trace
            S_stream = WINDOWS.Filter(S_stream, self.prior['P_HP'], self.prior['P_LP'])
        return P_stream,S_stream, or_P_len, or_S_len



    def normalize(self, v):
        """This function fails when division by zero"""
        normalized_v = v / np.sqrt(np.sum(v ** 2))
        return normalized_v

    def paper_plot(self,real_v):
        params = {'legend.fontsize': 'x-large',
                  'figure.figsize': (15, 15),
                  'axes.labelsize': 25,
                  'axes.titlesize': 'x-large',
                  'xtick.labelsize': 25,
                  'ytick.labelsize': 25}
        pylab.rcParams.update(params)
        df_select = self.df[["Epi", "Depth", "Strike", "Dip", "Rake"]]
        strike,dip,rake = aux_plane(real_v[2],real_v[3],real_v[4])



        fig = plt.figure(figsize=(25,6))

        row = 7
        col = 3
        ax1 = plt.subplot2grid((row,col),(1,0))
        plt.hist(df_select['Strike'],bins=100, alpha = 0.8)

        ymin, ymax = ax1.get_ylim()
        plt.vlines(real_v[2], ymin=ymin, ymax=ymax, colors='g', linewidth=3,label = 'Auxiliary plane')
        plt.vlines(strike, ymin=ymin, ymax=ymax, colors='k', linewidth=3, label='Fault plane')
        ax1.set_title("Density Strike", color='b', fontsize=25)
        ax1.set_xlabel("N=%i" % (len(df_select['Strike'])), fontsize=25)
        ax1.set_ylabel("Posterior marginal" , fontsize=25)
        ax1.tick_params(axis='x', labelsize=20)
        ax1.tick_params(axis='y', labelsize=20)
        ax1.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        ax1.set_xlim(0,360)
        plt.tight_layout()

        ax2 = plt.subplot2grid((row,col),(1,1))
        plt.hist(df_select['Dip'],bins=100, alpha = 0.8)
        ymin, ymax = ax2.get_ylim()
        plt.vlines(real_v[3], ymin=ymin, ymax=ymax, colors='g', linewidth=3)
        plt.vlines(dip, ymin=ymin, ymax=ymax, colors='k', linewidth=3,label='True model')
        ax2.set_title("Density Dip", color='b', fontsize=25)
        ax2.set_xlabel("N=%i" % (len(df_select['Dip'])), fontsize=25)
        ax2.xaxis.set_ticks(np.arange(0, 90, 10))
        ax2.tick_params(axis='x', labelsize=20)
        ax2.tick_params(axis='y', labelsize=20)
        ax2.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        ax2.set_xlim(0, 90)
        plt.tight_layout()

        ax3 = plt.subplot2grid((row,col),(1,2))
        plt.hist(df_select['Rake'],bins=100,alpha = 0.8)
        ymin, ymax = ax3.get_ylim()
        plt.vlines(real_v[4], ymin=ymin, ymax=ymax, colors='g', linewidth=3, label = 'Auxiliary plane')
        plt.vlines(rake, ymin=ymin, ymax=ymax, colors='k', linewidth=3,label='Fault plane')
        ax3.set_title("Density Rake", color='b', fontsize=25)
        ax3.set_xlabel("N=%i" % (len(df_select['Rake'])), fontsize=25)
        ax3.tick_params(axis='x', labelsize=20)
        ax3.tick_params(axis='y', labelsize=20)
        ax3.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        ax3.set_xlim(-180, 180)
        plt.legend(loc = 'upper left',fontsize=20)
        plt.tight_layout()

        ax4 = plt.subplot2grid((row,col),(2,1),rowspan= 5)
        plt.hist(df_select['Depth']/1000,bins=100,alpha=0.8)
        ymin, ymax = ax4.get_ylim()
        xmin, xmax = ax4.get_xlim()
        plt.vlines(real_v[1], ymin=ymin, ymax=ymax, colors='k', linewidth=3, label='True model')
        ax4.set_xlabel("Depth [km]", fontsize=20)
        ax4.set_ylabel("Posterior marginal", fontsize=20)
        ax4.tick_params(axis='x', labelsize=20)
        ax4.tick_params(axis='y', labelsize=20)
        ax4.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        ax4.ticklabel_format(style="sci", axis='x', scilimits=(-2, 2))

        # ax4.set_xlim(real_v[1]-10000, real_v[1]+10000)
        ax4.set_xlim(xmin, xmax)
        # plt.legend( fontsize=20)
        plt.tight_layout()

        ax5 = plt.subplot2grid((row,col),(2,2), rowspan = 5)
        plt.hist(df_select['Epi'],bins=50,alpha=0.8)
        ymin, ymax = ax5.get_ylim()
        plt.vlines(real_v[0], ymin=ymin, ymax=ymax, colors='k', linewidth=3,label = 'True model')
        ax5.set_title("Density Epicentral Distance", color='b', fontsize=20)
        ax5.set_xlabel("N=%i" % (len(df_select['Epi'])), fontsize=20)
        ax5.tick_params(axis='x', labelsize=20)
        ax5.tick_params(axis='y', labelsize=20)
        ax5.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        # plt.legend(fontsize=20)
        plt.tight_layout()

        ax6 = plt.subplot2grid((row,col),(0,0))
        Mw = 2.0 / 3.0 * (np.log10(self.df['M0']) - 9.1)
        plt.hist(Mw,bins=np.arange(4,7,0.05),alpha=0.8)
        ymin, ymax = ax6.get_ylim()
        Mw = self.Scalarmoment2Magnitude(self.prior['M0'])
        plt.vlines(Mw, ymin=ymin, ymax=ymax, colors='k', linewidth=3, label= 'True model')
        ax6.set_title("Density Moment Magnitude", color='b', fontsize=20)
        ax6.set_xlabel("N=%i" % (len(self.df['M0'])), fontsize=18)
        ax6.tick_params(axis='x', labelsize=18)
        ax6.tick_params(axis='y', labelsize=18)
        ax6.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        # ax6.set_xscale('log')
        # ax6.set_xlim(0e18, 0.5e18 )
        # plt.legend(loc='center left', bbox_to_anchor=(1.1, 0.5))
        plt.legend(loc='upper right',fontsize = 20)
        plt.tight_layout()

        ax7 = plt.subplot2grid((row,col),(0,1))
        # b = beach([real_v[2], real_v[3], real_v[4]], size=200, linewidth=2, facecolor='b')
        xmin,xmax = ax7.get_xlim()
        ymin,ymax = ax7.get_ylim()
        b = beach([real_v[2], real_v[3], real_v[4]], width=100, linewidth=0, facecolor='b',xy=(0.5 * xmax, 0.5* ymax), axes=ax7)
        ax7.add_collection(b)
        ax7.set_xticks([])
        ax7.set_yticks([])
        ax7.axis('off')


        ax8 = plt.subplot2grid((row,col),(0,2))
        fig_bb, ax_bb = plt.subplots(1, 1, figsize=(4, 4))

        ax_bb.set_xticks([])
        ax_bb.set_yticks([])
        ax_bb.axis('off')
        img = None
        buf = io.BytesIO()

        strike = self.df['Strike']
        dip = self.df['Dip']
        rake = self.df['Rake']

        # for i in range(0, len(strike)):
        for i in range(0, len(strike),10000):

            b = beach(fm=[strike[i], dip[i], rake[i]],
                      width=200, linewidth=0, facecolor='b',
                      xy=(0, 0), axes=ax_bb, alpha=1, zorder=i)
            ax_bb.add_collection(b)
            ax_bb.set_xlim((-0.1, 0.1))
            ax_bb.set_ylim((-0.1, 0.1))

            p = Circle((0., 0,), 0.065, linewidth=2, edgecolor='k', facecolor='b', zorder=5)
            ax_bb.add_patch(p)

            buf.seek(0)
            fig_bb.savefig(buf, format='png')
            buf.seek(0)
            if img is None:
                img = mpimg.imread(buf)
            else:
                img += mpimg.imread(buf)
        plt.close(fig_bb)
        ax8.imshow(img / np.max(img.flatten()))
        ax8.set_xticks([])
        ax8.set_yticks([])
        ax8.axis('off')
        plt.show()


        ax9 = plt.subplot2grid((row,col),(2,0))
        inner = gridspec.GridSpecFromSubplotSpec(5,1,ax9)

        epi = self.df['Epi']
        depth = self.df['Depth']
        strike = self.df['Strike']
        dip = self.df['Dip']
        rake = self.df['Rake']
        M0 = self.df['M0']

        P_shift = self.df['Shift_P']
        S_shift = self.df['Shift_S']


        seis = Get_Seismogram(self.prior)
        BW_syn = Cut_windows(self.prior['VELOC_taup'])

        delta = self.BW_obs.P_stream.traces[0].meta.delta
        p_time_array = np.arange(len(self.BW_obs.P_stream.traces[0].data)) * delta
        s_time_array = np.arange(len(self.BW_obs.S_stream.traces[0].data)) * delta

        start_P = int((self.BW_obs.start_P.timestamp - self.otime.timestamp - 50) / delta)
        end_P = int((self.BW_obs.start_P.timestamp - self.otime.timestamp + 130) / delta)

        start_S = int((self.BW_obs.start_S.timestamp - self.otime.timestamp - 50) / delta)
        end_S = int((self.BW_obs.start_S.timestamp - self.otime.timestamp + 100) / delta)


        a= plt.subplot2grid((row,col), (2, 0))
        b= plt.subplot2grid((row,col), (3, 0))
        c= plt.subplot2grid((row,col), (4, 0))
        d= plt.subplot2grid((row,col), (5, 0))
        e= plt.subplot2grid((row,col), (6, 0))
        # a = np.where(np.logical_and(epi>=prior['Real_epi'] , epi <= prior['Real_epi']))
        # a = np.where(np.logical_and(epi>=88 -0.5  , epi <= 88 + 0.5))

        # for i in range(len(a[0])):
        # for i in np.arange(len(epi) - 100, len(epi), 1):
        for i in np.arange(2):
            dict = geo.Geodesic(a=self.prior['radius'], f=self.prior['f']).ArcDirect(lat1=self.prior['la_r'], lon1=self.prior['lo_r'],
                                                                           azi1=self.prior['baz'], a12=epi[i],
                                                                           outmask=1929)

            st_syn = seis.get_seis_manual(la_s=dict['lat2'], lo_s=dict['lon2'], depth=depth[i],
                                          strike=strike[i], dip=dip[i], rake=rake[i],
                                          time=self.otime, M0=M0[i])

            BW_syn.Get_bw_windows(st_syn, epi[i], depth[i], self.otime, self.prior['npts'])

            P_shift_array = self.shift(BW_syn.P_stream.traces[0].data, -int(P_shift[i]))
            #
            # ax1.plot(p_time_array[start_P:end_P], BW_syn.P_stream.traces[0].data[start_P:end_P], 'g',
            #          label='Synthetic', linewidth = 0.1)

            a.plot(p_time_array[start_P:end_P], P_shift_array[start_P:end_P], 'r',
                     label='Synthetic', linewidth=0.1)
            # ax1.plot( P_shift_array, 'r',  label='Synthetic', linewidth = 0.1)
            a.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
            a.tick_params(axis='x', labelsize=18)
            a.tick_params(axis='y', labelsize=18)

            # plt.tight_layout()
            # plt.legend(loc='lower left', fontsize=15)

            P_shift_array = self.shift(BW_syn.P_stream.traces[1].data, -int(P_shift[i]))
            b.plot(p_time_array[start_P:end_P], P_shift_array[start_P:end_P], 'r',
                     label='Synthetic', linewidth=0.1)
            # ax2.plot(P_shift_array, 'r',label='Synthetic', linewidth = 0.1)
            b.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
            b.tick_params(axis='x', labelsize=18)
            b.tick_params(axis='y', labelsize=18)
            # plt.tight_layout()

            S_shift_array = self.shift(BW_syn.S_stream.traces[0].data, -int(S_shift[i]))
            c.plot(s_time_array[start_S:end_S], S_shift_array[start_S:end_S], 'r', linewidth=0.1)
            # ax3.plot(S_shift_array, 'r', linewidth = 0.1)
            c.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
            c.tick_params(axis='x', labelsize=18)
            c.tick_params(axis='y', labelsize=18)
            # plt.tight_layout()

            S_shift_array = self.shift(BW_syn.S_stream.traces[1].data, -int(S_shift[i]))
            d.plot(s_time_array[start_S:end_S], S_shift_array[start_S:end_S], 'r', linewidth=0.1)
            # ax4.plot(S_shift_array, 'r', linewidth = 0.1)
            d.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
            d.tick_params(axis='x', labelsize=18)
            d.tick_params(axis='y', labelsize=18)
            # plt.tight_layout()

            S_shift_array = self.shift(BW_syn.S_stream.traces[2].data, -int(S_shift[i]))
            e.plot(s_time_array[start_S:end_S], S_shift_array[start_S:end_S], 'r', linewidth=0.1)
            # ax5.plot( S_shift_array, 'r', linewidth = 0.1)
            e.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
            e.tick_params(axis='x', labelsize=18)
            e.tick_params(axis='y', labelsize=18)
            e.set_xlabel(self.BW_obs.start_P.strftime('From P arrival: %Y-%m-%dT%H:%M:%S + [sec]'), fontsize=18)

        a.plot(p_time_array[start_P:end_P], self.BW_obs.P_stream.traces[0].data[start_P:end_P], 'b')
        # ax1.plot( BW_obs.P_stream.traces[0].data, 'b', label='Observed', linewidth = 0.1)
        ymin, ymax = a.get_ylim()
        xmin, xmax = a.get_xlim()
        a.text(xmax - 5, ymax / 1.7, "P-Z", fontsize=20, color='b')

        b.plot(p_time_array[start_P:end_P], self.BW_obs.P_stream.traces[1].data[start_P:end_P], 'b')
        # ax2.plot(BW_obs.P_stream.traces[1].data, 'b', linewidth = 0.1)
        ymin, ymax = b.get_ylim()
        xmin, xmax = b.get_xlim()
        b.text(xmax - 5, ymax / 1.7, "P-R", fontsize=20, color='b')

        c.plot(s_time_array[start_S:end_S], self.BW_obs.S_stream.traces[0].data[start_S:end_S], 'b')
        # ax3.plot( BW_obs.S_stream.traces[0].data, 'b', linewidth = 0.1)
        ymin, ymax = c.get_ylim()
        xmin, xmax = c.get_xlim()
        c.text(xmax - 10, ymax / 1.5, "S-Z", fontsize=20, color='b')

        d.plot(s_time_array[start_S:end_S], self.BW_obs.S_stream.traces[1].data[start_S:end_S], 'b')
        # ax4.plot(BW_obs.S_stream.traces[1].data, 'b', linewidth = 0.1)
        ymin, ymax = d.get_ylim()
        xmin, xmax = d.get_xlim()
        d.text(xmax - 10, ymax / 1.5, "S-R", fontsize=20, color='b')

        e.plot(s_time_array[start_S:end_S], self.BW_obs.S_stream.traces[2].data[start_S:end_S], 'b')
        # ax5.plot(BW_obs.S_stream.traces[2].data, 'b', linewidth = 0.1)
        ymin, ymax = e.get_ylim()
        xmin, xmax = e.get_xlim()
        e.text(xmax - 10, ymax / 1.7, "S-T", fontsize=20, color='b')

        plt.show()
        plt.savefig(self.dir + '/Total_plot.pdf')



    def Scalarmoment2Magnitude(self,M0):
        Mw = 2.0 / 3.0 * (np.log10(M0) - 9.1)
        return Mw

    def shift(self, np_array, time_shift):
        new_array = np.zeros_like(np_array)
        if time_shift < 0:
            new_array[-time_shift:] = np_array[:time_shift]
        elif time_shift == 0:
            new_array[:] = np_array[:]
        else:
            new_array[:-time_shift] = np_array[time_shift:]
        return new_array