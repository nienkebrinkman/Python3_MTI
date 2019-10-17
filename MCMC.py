import numpy as np
import obspy
import geographiclib.geodesic as geo
import matplotlib.pylab as plt
import pylab
import os

from obspy.core.stream import Stream

from Model_samples import Model_samples
from Cut_windows import Cut_windows
from Get_Seismogram import Get_Seismogram
from Misfit import Misfit

class MCMC:
    def __init__(self,or_time,PRIOR,sample_path=None):
        self.or_time = or_time
        self.prior = PRIOR
        self.sample_path = sample_path

        self.model = Model_samples(self.prior)
        self.BW_syn = Cut_windows(self.prior['VELOC_taup'],P_HP = PRIOR['P_HP'], P_LP= PRIOR['P_LP'], S_HP =  PRIOR['S_HP'], S_LP= PRIOR['S_LP'])
        self.seis = Get_Seismogram(self.prior)
        self.mis = Misfit()

    def start_BW(self,BW_obs):
        self.BW_obs = BW_obs
        self.start_MCMC(self.prior['save_dir'] + '/%s.txt' % self.prior['save_name'],BW=True)


    def start_SW(self):
        raise ValueError('Surface waves Not implemented yet')

    def start_Bw_and_SW(self):
        raise ValueError('Surface waves Not implemented yet')


    def start_MCMC(self,savepath,BW = False,SW=False):
        if BW ==False and SW == False:
            raise ValueError('Run either: start_SW (Surface waves) or start_BW (Body waves)')


        with open(savepath, 'w') as save_file:
            self.write_par(save_file,BW, SW)
            accepted = 0
            rejected = 0
            for self.i in range(self.prior['sample_number']):
                if self.i % 10 == 0:
                    print("proposal: %i, accepted: %i" % (self.i, accepted))
                    print("proposal: %i, rejected: %i" % (self.i, rejected))

                if self.i == 0:
                    if self.sample_path == None:
                        update = None
                        epi, depth = self.model.model_samples()
                        M0_old = self.prior['M0']
                        strike, dip, rake = self.model.model_samples_sdr()
                        moment_old = moment_old = np.array([strike, dip, rake])
                    else:
                        data = np.loadtxt(self.sample_path, delimiter=',')
                        epi = data[0]
                        depth = data[1]
                        strike = data[2]
                        dip = data[3]
                        rake = data[4]
                        M0_old = data[5]

                else:
                    update = np.random.choice(['epi', 'depth', 'moment'], 1)[0]
                    epi, depth = self.model.model_samples(update, epi_old, depth_old)
                    if epi < self.prior['epi']['range_min'] or epi > self.prior['epi']['range_max'] or depth < \
                            self.prior['depth']['range_min'] or depth > self.prior['depth']['range_max']:
                        continue

                    strike, dip, rake = self.model.model_samples_sdr(update,strike_old,dip_old,rake_old)

                self.G_function(epi,depth,M0_old,strike,dip,rake)

                if BW == True:
                    Xi_bw_new, amplitude,time_shift_new, fig = self.mis.CC_BW(self.BW_obs,self.BW_syn,self.or_time,self.prior['PLOT'])

                    s_z_new =  Xi_bw_new[0] # 1 *
                    s_r_new =  Xi_bw_new[1] # 1 *
                    s_t_new =  Xi_bw_new[2] # 5 *
                    p_z_new =  Xi_bw_new[3] # 5 *
                    p_r_new =  Xi_bw_new[4] # 2 *
                    bw_new = s_z_new + s_r_new + s_t_new + p_z_new + p_r_new
                    Xi_new = bw_new
                if SW == True:
                    # sw_new =
                    raise ValueError('Surface waves Not implemented yet')
                if BW == True and SW ==True:
                    raise ValueError('Surface waves Not implemented yet')
                    # Xi_new = bw_new + sw_new

                M0 = M0_old / np.mean(amplitude)

                if self.i == 0:
                    if self.prior['PLOT'] == True and self.i % 1 == 0:
                        # self.plot()
                        if not os.path.exists(self.prior['save_dir'] + '/plots/'):
                            os.makedirs(self.prior['save_dir'] + '/plots/')
                        fig.savefig(self.prior['save_dir'] + '/plots/SHIFT_%s_%05i.png' % (self.prior['save_name'], self.i))
                        plt.close("all")
                    if BW == True:
                        self.s_z_old = s_z_new
                        self.s_r_old = s_r_new
                        self.s_t_old = s_t_new
                        self.p_z_old = p_z_new
                        self.p_r_old = p_r_new
                        self.bw_old = bw_new
                    if SW == True:
                        raise ValueError('Surface waves Not implemented yet')
                    if BW == True and SW == True:
                        raise ValueError('Surface waves Not implemented yet')
                    Xi_old = Xi_new
                    epi_old = epi
                    depth_old = depth
                    strike_old = strike
                    dip_old = dip
                    rake_old = rake
                    M0_old = M0
                    time_shift_old = time_shift_new

                    self.write_sample(save_file, epi_old, depth_old, strike_old, dip_old, rake_old, M0_old, Xi_old, BW,
                                      SW,time_shift_old, accept=1)

                    continue
                random = np.random.random_sample((1,))
                if (Xi_new < Xi_old) or (np.exp((Xi_old - Xi_new) / self.prior['Temperature']) > random):
                    if self.prior['PLOT'] == True and self.i % 1 == 0:
                        # self.plot()
                        fig.savefig(self.prior['save_dir'] + '/plots/SHIFT_%s_%05i.png' % (self.prior['save_name'], self.i))
                        plt.close("all")
                    print(Xi_new)
                    if BW == True:
                        self.s_z_old = s_z_new
                        self.s_r_old = s_r_new
                        self.s_t_old = s_t_new
                        self.p_z_old = p_z_new
                        self.p_r_old = p_r_new
                        self.bw_old = bw_new
                    if SW == True:
                        raise ValueError('Surface waves Not implemented yet')
                    if BW == True and SW == True:
                        raise ValueError('Surface waves Not implemented yet')
                    Xi_old = Xi_new
                    epi_old = epi
                    depth_old = depth
                    strike_old = strike
                    dip_old = dip
                    rake_old = rake
                    M0_old = M0
                    time_shift_old = time_shift_new
                    accepted = accepted + 1
                    if self.i%10 ==0:
                        self.write_sample(save_file, epi_old, depth_old, strike_old,dip_old,rake_old, M0_old, Xi_old,BW,SW,time_shift_old,accept=1)
                else:
                    if self.prior['PLOT']:
                        plt.close("all")
                    rejected += 1
                    if self.i % 10 == 0:
                        self.write_sample(save_file, epi_old, depth_old, strike_old,dip_old,rake_old, M0_old, Xi_old,BW,SW,time_shift_old, accept=0)
            save_file.close()

    def G_function(self, epi, depth, M0, strike, dip, rake ):

        dict = geo.Geodesic(a=self.prior['radius'], f=self.prior['f']).ArcDirect(lat1=self.prior['la_r'],
                                                                                 lon1=self.prior['lo_r'],
                                                                                 azi1=self.prior['baz'],
                                                                                 a12=epi, outmask=1929)



        st_syn = self.seis.get_seis_manual(la_s=dict['lat2'], lo_s=dict['lon2'], depth=depth,
                                                               strike=strike, dip=dip, rake=rake,
                                                               time=self.or_time, M0=M0)

        self.BW_syn.Get_bw_windows(st_syn, epi, depth, self.or_time, self.prior['npts'])
        # ax1 = plt.subplot(111)
        # plt.plot(self.BW_obs.original.traces[0], 'b')
        # plt.plot(self.BW_syn.original.traces[0], 'r')
        # ymin, ymax = ax1.get_ylim()
        # plt.vlines(self.BW_obs.or_P_len , ymin=ymin, ymax=ymax, colors='b', linewidth=3, label='Obs_P')
        # plt.vlines(self.BW_syn.or_P_len, ymin=ymin, ymax=ymax, colors='r', linewidth=3, label = 'Syn_P')
        # plt.vlines(self.BW_obs.or_S_len , ymin=ymin, ymax=ymax, colors='b', linewidth=3, label='Obs_S')
        # plt.vlines(self.BW_syn.or_S_len, ymin=ymin, ymax=ymax, colors='r', linewidth=3, label = 'Syn_S')
        # plt.show()

    def write_par(self,file_name,BW,SW):
        if BW == True and SW == True:
            raise ValueError('Surface waves Not implemented yet')
        elif BW == True:
            file_name.write("epi, depth, strike, dip, rake, M0, Total-misfit, p_z, p_r, s_z, s_r, s_t, bw_tot, shift_S, shift-P, Iteration\n\r")
            file_name.write("Velocity Model:%s\n\r" %self.prior['VELOC'])
            file_name.write("Station:%s\n\r" % self.BW_obs.P_stream.traces[0].stats.station)
            file_name.write("Sampling rate:%.2f\n\r" % self.BW_obs.P_stream.traces[0].stats.sampling_rate)
            file_name.write("la_r:%.4f\n\r" % self.prior['la_r'])
            file_name.write("lo_r:%.4f\n\r" % self.prior['lo_r'])
            file_name.write("kind:%s\n\r" % self.prior['kind'])  #
            file_name.write("network:%s\n\r" % self.prior['network'])  #
            file_name.write("amount samples:%i\n\r" % self.prior['sample_number'])  #
            file_name.write("Temperature:%i\n\r" % self.prior['Temperature'])  #
            file_name.write("Radius:%.4f\n\r" % self.prior['radius'])  #
            file_name.write("Flattening:%.4f\n\r" % self.prior['f'])  #
            file_name.write("Azimuth:%.4f\n\r" % self.prior['az'])  #

        elif SW == True:
            raise ValueError('Surface waves Not implemented yet')

    def write_sample(self,file_name, epi, depth, strike,dip,rake, M0, Xi_old,BW,SW,time_shift,accept=0):
        file_name.write("%.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, " % (epi, depth, strike, dip, rake, M0,Xi_old))
        if BW == True:
            file_name.write("%.4f, %.4f, %.4f, %.4f, %.4f, %.4f, " % (self.p_z_old, self.p_r_old, self.s_z_old, self.s_r_old, self.s_t_old,self.bw_old))
        if SW == True:
            raise ValueError('Surface waves Not implemented yet')
        file_name.write("%i, %i, %i\n\r" % (time_shift[0],time_shift[1],self.i))

    def plot(self):
        stream = Stream()
        stream += self.BW_obs.BW_stream.traces[0]
        # stream += self.BW_obs.S_stream.traces[1]
        # stream += self.BW_obs.S_stream.traces[2]
        #
        #
        stream += self.BW_syn.BW_stream.traces[0]
        # stream += self.BW_syn.S_stream.traces[1]
        # stream += self.BW_syn.S_stream.traces[2]
        #
        stream.plot()

        fig = plt.figure(figsize=(10,12))
        delta = self.BW_obs.P_stream.traces[0].meta.delta
        p_time_array = np.arange(len(self.BW_obs.P_stream.traces[0].data)) * delta
        s_time_array = np.arange(len(self.BW_obs.S_stream.traces[0].data)) * delta

        start_P = int((self.BW_obs.start_P.timestamp - self.or_time.timestamp - 10) /delta )
        end_P = int((self.BW_obs.start_P.timestamp - self.or_time.timestamp + 30) /delta)

        start_S = int((self.BW_obs.start_S.timestamp - self.or_time.timestamp - 20)/delta)
        end_S = int((self.BW_obs.start_S.timestamp - self.or_time.timestamp + 100)/delta)

        ax1 = plt.subplot2grid((5,1),(0,0))
        plt.plot(p_time_array[start_P:end_P],self.BW_obs.P_stream.traces[0].data[start_P:end_P],'b', label = 'Observed')
        plt.plot(p_time_array[start_P:end_P],self.BW_syn.P_stream.traces[0].data[start_P:end_P],'r', label = 'Synthetic')
        ymin, ymax = ax1.get_ylim()
        xmin, xmax = ax1.get_xlim()
        plt.text(xmax-5, ymax/1.7, "P-Z", fontsize=20, color='b')
        ax1.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        ax1.tick_params(axis='x', labelsize=18)
        ax1.tick_params(axis='y', labelsize=18)
        plt.tight_layout()
        plt.legend(loc='lower left',fontsize=15)

        ax2 = plt.subplot2grid((5,1),(1,0))
        plt.plot(p_time_array[start_P:end_P],self.BW_obs.P_stream.traces[1].data[start_P:end_P],'b')
        plt.plot(p_time_array[start_P:end_P],self.BW_syn.P_stream.traces[1].data[start_P:end_P],'r')
        ymin, ymax = ax2.get_ylim()
        xmin, xmax = ax2.get_xlim()
        plt.text(xmax-5, ymax/1.7, "P-R", fontsize=20, color='b')
        ax2.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        ax2.tick_params(axis='x', labelsize=18)
        ax2.tick_params(axis='y', labelsize=18)
        plt.tight_layout()

        ax3 = plt.subplot2grid((5,1),(2,0))
        plt.plot(s_time_array[start_S:end_S],self.BW_obs.S_stream.traces[0].data[start_S:end_S],'b')
        plt.plot(s_time_array[start_S:end_S],self.BW_syn.S_stream.traces[0].data[start_S:end_S],'r')
        ymin, ymax = ax3.get_ylim()
        xmin, xmax = ax3.get_xlim()
        plt.text(xmax-10, ymax/1.7, "S-Z", fontsize=20, color='b')
        ax3.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        ax3.tick_params(axis='x', labelsize=18)
        ax3.tick_params(axis='y', labelsize=18)
        plt.tight_layout()

        ax4 = plt.subplot2grid((5,1),(3,0))
        plt.plot(s_time_array[start_S:end_S],self.BW_obs.S_stream.traces[1].data[start_S:end_S],'b')
        plt.plot(s_time_array[start_S:end_S],self.BW_syn.S_stream.traces[1].data[start_S:end_S],'r')
        ymin, ymax = ax4.get_ylim()
        xmin, xmax = ax4.get_xlim()
        plt.text(xmax-10, ymax/1.7, "S-R", fontsize=20, color='b')
        ax4.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        ax4.tick_params(axis='x', labelsize=18)
        ax4.tick_params(axis='y', labelsize=18)
        plt.tight_layout()

        ax5 = plt.subplot2grid((5,1),(4,0))
        plt.plot(s_time_array[start_S:end_S],self.BW_obs.S_stream.traces[2].data[start_S:end_S],'b')
        plt.plot(s_time_array[start_S:end_S],self.BW_syn.S_stream.traces[2].data[start_S:end_S],'r')
        ymin, ymax = ax5.get_ylim()
        xmin, xmax = ax5.get_xlim()
        plt.text(xmax-10, ymax/1.9, "S-T", fontsize=20, color='b')
        ax5.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        ax5.tick_params(axis='x', labelsize=18)
        ax5.tick_params(axis='y', labelsize=18)
        ax5.set_xlabel(self.BW_obs.start_P.strftime('From P arrival: %Y-%m-%dT%H:%M:%S + [sec]'),fontsize= 18)
        plt.tight_layout()
        # plt.show()
        plt.savefig(self.prior['save_dir']+ '/plots/%s_%i.png' %( self.prior['save_name'],self.i))
        # plt.show()
        plt.close()





