import numpy as np
import obspy.signal.cross_correlation as cc
import matplotlib.pylab as plt
from obspy.core.stream import Stream

class Misfit:
    def L2_BW(self, BW_obs, BW_syn, or_time, var):
        p_obs = BW_obs.P_stream
        p_syn = BW_syn.P_stream
        s_obs = BW_obs.S_stream
        s_syn = BW_syn.S_stream
        dt = s_obs[0].meta.delta
        misfit = np.array([])
        time_shift = np.array([], dtype=int)
        # S - correlations:
        for i in range(len(s_obs)):
            cc_obspy = cc.correlate(s_obs[i].data, s_syn[i].data, int(0.25* len(s_obs[i].data)))
            shift, CC_s = cc.xcorr_max(cc_obspy)

            s_syn_shift = self.shift(s_syn[i].data, -shift)
            time_shift = np.append(time_shift, shift)

            # d_obs_mean = np.mean(s_obs[i].data)
            # var_array = var * d_obs_mean

            var_array = np.var(s_obs[i].data)
            # var_array = var**2

            misfit = np.append(misfit, np.matmul((s_obs[i].data - s_syn_shift).T, (s_obs[i].data - s_syn_shift)) / (
                2 * (var_array)))
            # time = -time_shift * dt  # Relatively, the s_wave arrives now time later or earlier than it originally did

            # plt.plot(s_syn_shift,label='s_shifted',linewidth=0.3)
            # plt.plot(s_syn[i], label='s_syn',linewidth=0.3)
            # plt.plot(s_obs[i], label='s_obs',linewidth=0.3)
            # plt.legend()
            # plt.savefig()
            # plt.close()
        # P- correlation
        for i in range(len(p_obs)):
            cc_obspy = cc.correlate(s_obs[i].data, s_syn[i].data, int(0.25 * len(p_obs[i].data)))
            shift, CC_p = cc.xcorr_max(cc_obspy)

            p_syn_shift = self.shift(p_syn[i].data, -shift)
            time_shift = np.append(time_shift, shift)

            # d_obs_mean = np.mean(p_obs[i].data)
            # var_array = var * d_obs_mean
            var_array = np.var(p_obs[i].data)

            misfit = np.append(misfit, np.matmul((p_obs[i].data - p_syn_shift).T, (p_obs[i].data - p_syn_shift)) / (
                2 * (var_array)))
            # time = -time_shift + len(s_obs)] * dt  # Relatively, the s_wave arrives now time later or earlier than it originally did

            plt.plot(p_syn_shift, label='p_shifted')
            plt.plot(p_syn[i], label='p_syn')
            plt.plot(p_obs[i], label='p_obs')
            plt.legend()
            plt.show()
            # plt.close()
        sum_misfit = np.sum(misfit)
        return sum_misfit, time_shift

    def CC_BW(self,BW_obs,BW_syn,or_time,plot = False):
        # fig = plt.figure(figsize=(12, 10))
        # ax1 = plt.subplot2grid((2, 1), (0, 0))
        # time_array = np.arange(len(BW_obs.BW_stream[0].data)) * BW_obs.BW_stream[0].stats.delta
        # ymin, ymax = ax1.get_ylim()
        # ax1.plot(time_array[0:4000],BW_obs.BW_stream.traces[0][0:4000], alpha = 0.5,c = 'r', label = 'Observed')
        # ax1.plot(time_array[0:4000], BW_syn.BW_stream.traces[0][0:4000], alpha = 0.5,c= 'g', label = 'Synthetic')
        # ax1.tick_params(axis='x', labelsize=20)
        # ax1.tick_params(axis='y', labelsize=20)
        # plt.xlabel(or_time.strftime('Time : %Y-%m-%dT%H:%M:%S + [sec]'),fontsize = 20)
        # plt.legend(loc='upper right')
        #
        #
        # ax2 = plt.subplot2grid((2, 1), (1, 0))
        # plt.vlines(BW_obs.start_P.timestamp-or_time.timestamp, ymin=ymin, ymax=ymax, colors='r', linewidth=3, label = 'Obs_P')
        # plt.vlines(BW_syn.start_P.timestamp-or_time.timestamp, ymin=ymin, ymax=ymax, colors='g', linewidth=3, label = 'Syn_P')
        # plt.vlines(BW_obs.start_S.timestamp-or_time.timestamp, ymin=ymin, ymax=ymax, colors='r', linewidth=3, label = 'Obs_S')
        # plt.vlines(BW_syn.start_S.timestamp-or_time.timestamp, ymin=ymin, ymax=ymax, colors='g', linewidth=3, label = 'Syn_P')
        # ax2.tick_params(axis='x', labelsize=20)
        # ax2.tick_params(axis='y', labelsize=20)
        # plt.legend(loc ='upper center')
        # plt.xlabel(or_time.strftime('Time : %Y-%m-%dT%H:%M:%S + [sec]'),fontsize = 20)
        # plt.tight_layout()
        #
        #
        #
        # plt.show()





        p_obs = BW_obs.P_stream
        p_syn = BW_syn.P_stream
        s_obs = BW_obs.S_stream
        s_syn = BW_syn.S_stream
        p_start_obs = BW_obs.start_P
        p_start_syn = BW_syn.start_P
        s_start_obs = BW_obs.start_S
        s_start_syn = BW_syn.start_S

        dt = s_obs[0].meta.delta
        misfit = np.array([])
        misfit_obs = np.array([])
        time_shift = np.array([], dtype=int)
        amplitude = np.array([])
        #
        # amp_obs = p_obs.copy()
        # amp_obs.trim(p_start_obs, p_start_obs + 30)
        # amp_syn = p_syn.copy()
        # amp_syn.trim(p_start_syn, p_start_syn + 30)

        if plot:
            fig = plt.figure(figsize=(10, 12))
        else:
            fig = 1


        #
        # t1 = s_start_obs
        # t2 = s_start_syn
        #
        # tr1 = BW_obs.S_stream.select(component="Z")[0]
        # tr2 = BW_syn.S_stream.select(component="Z")[0]
        #
        # dt, coeff = cc.xcorr_pick_correction(t1, tr1, t2, tr2, 5, 20, 10, plot=True)

        # S - correlations:
        # Calculate Shift based on T component
        len_S_obs = len(s_obs[2].data)
        cc_obspy = cc.correlate(s_obs[2].data[0:len_S_obs],
                                s_syn[2].data[0:len_S_obs],
                                int(len_S_obs))
        shift_centered, MAX_CC = cc.xcorr_max(cc_obspy, abs_max=False)
        shift = np.argmax(cc_obspy)
        time_shift = np.append(time_shift, shift_centered)

        mu_s = np.array([0.6, 0.6, 0.9])
        sigma_s = np.array([0.3, 0.3, 0.1])

        for i in range(len(s_obs)):
            delta = s_obs[i].stats.delta
            cc_obspy = cc.correlate(s_obs[i].data[0:len_S_obs],
                                    s_syn[i].data[0:len_S_obs],
                                    int(len_S_obs))
            CC_s = cc_obspy[shift]
            #shift, CC_s = cc.xcorr_max(cc_obspy, abs_max=False)

            s_syn_shift_obspy = self.shift(s_syn[i].data, -shift_centered)


            misfit = np.append(misfit, ((CC_s - mu_s[i]) ** 2) / (2 * (sigma_s[i]) ** 2))# + np.abs(shift))

            if plot:
                time_array = np.arange(len(s_obs[i].data)) * delta
                start = int((s_start_obs.timestamp - or_time.timestamp - 100) / delta)
                end =int((s_start_obs.timestamp - or_time.timestamp + 400) / delta)

                ax1 = plt.subplot2grid((5, 1), (i, 0))
                plt.plot(time_array[start:end],s_obs[i][start:end],'b', label = 'Observed')
                plt.plot(time_array[start:end],s_syn[i][start:end],'r', linewidth = 0.3, label = 'Synthetic')
                plt.plot(time_array[start:end],s_syn_shift_obspy[start:end],'g', label = 'Shifted')
                ymin, ymax = ax1.get_ylim()
                xmin, xmax = ax1.get_xlim()
                if i == 0:
                    plt.text(xmax - 30, ymax / 1.7, "S-Z - %.4f - %.4f " % (misfit[i] , CC_s), fontsize=20, color='b')
                elif i == 1:
                    plt.text(xmax - 30, ymax / 1.7, "S-R - %.4f - %.4f " % (misfit[i] , CC_s), fontsize=20, color='b')
                elif i == 2:
                    plt.text(xmax - 30, ymax / 1.7, "S-T - %.4f - %.4f " % (misfit[i] , CC_s), fontsize=20, color='b')
                ax1.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
                ax1.tick_params(axis='x', labelsize=18)
                ax1.tick_params(axis='y', labelsize=18)
                plt.tight_layout()

        len_P_obs = len(p_obs[0].data)
        cc_obspy = cc.correlate(p_obs[0].data[0:len_P_obs],
                                p_syn[0].data[0:len_P_obs],
                                int(len_P_obs))
        shift_centered, MAX_CC = cc.xcorr_max(cc_obspy, abs_max=False)
        shift = np.argmax(cc_obspy)
        time_shift = np.append(time_shift, shift_centered)

        mu_p = np.array([0.95, 0.9])
        sigma_p = np.array([0.1, 0.2])

        # P- correlation
        for i in range(len(p_obs)):
            # cc_obspy = cc.correlate(p_obs[i].data, p_syn[i].data, int(BW_obs.P_len))

            #cc_obspy = cc.correlate(p_obs[i].data, p_syn[i].data, int(len(p_obs[i].data)))
            # cc_obspy = cc.correlate(p_obs[i].data, p_syn[i].data, int(len(p_obs[i].data)* 0.5))
            len_P_obs = len(p_obs[i].data)
            cc_obspy = cc.correlate(p_obs[i].data[0:len_P_obs],
                                    p_syn[i].data[0:len_P_obs],
                                    int(len_P_obs))
            CC_p = cc_obspy[shift]
            # shift, CC_p = cc.xcorr_max(cc_obspy, abs_max=False)

            misfit = np.append(misfit, ((CC_p - mu_p[i] ) ** 2) / (2 * (sigma_p[i]) ** 2) )#+ np.abs(shift))

            p_syn_shift_obspy = self.shift(p_syn[i].data, -shift_centered)
            start = int((p_start_obs.timestamp - or_time.timestamp - 100) / delta)
            end = int((p_start_obs.timestamp - or_time.timestamp + 260) / delta)

            # A = (np.dot(amp_obs.traces[i],amp_syn.traces[i]) / np.dot(amp_obs.traces[i],amp_obs.traces[i]))
            A = (np.dot(p_obs[i].data[start:end], p_syn_shift_obspy[start:end]) /
                 np.dot(p_obs[i].data[start:end], p_obs[i].data[start:end]))
            amplitude = np.append(amplitude,abs(A))

            if plot:
                delta = p_obs[i].stats.delta
                time_array = np.arange(len_P_obs) * delta


                ax1 = plt.subplot2grid((5, 1), (i+3, 0))
                plt.plot(time_array[start:end],p_obs[i][start:end],'b', label = 'Observed')
                plt.plot(time_array[start:end],p_syn[i][start:end],'r', linewidth = 0.3, label = 'Synthetic')
                plt.plot(time_array[start:end],p_syn_shift_obspy[start:end],'g', label = 'Shifted')
                ymin, ymax = ax1.get_ylim()
                xmin, xmax = ax1.get_xlim()
                if i == 0:
                    plt.text(xmax - 30, ymax / 1.7, "P-Z - %.4f - %.4f " % (misfit[i+3] , CC_p), fontsize=20, color='b')
                elif i == 1:
                    plt.text(xmax - 30, ymax / 1.7, "P-R - %.4f - %.4f " % (misfit[i+3] , CC_p), fontsize=20, color='b')
                ax1.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
                ax1.tick_params(axis='x', labelsize=18)
                ax1.tick_params(axis='y', labelsize=18)
                plt.tight_layout()
            # plt.show()
            if i == 1 and plot == True:
                plt.legend(loc='lower left', fontsize=15)

        sum_misfit = np.sum(misfit)
        return misfit, amplitude, time_shift,fig

    def SW_L2(self, SW_env_obs, SW_env_syn, var, amplitude):
        misfit = np.array([])
        for i in range(len(SW_env_obs)):
            dt = SW_env_obs[i].meta.delta

            cc_obspy = cc.correlate(SW_env_obs[i].data, SW_env_syn[i].data,int( 0.25*len(SW_env_syn[i].data)))
            shift, CC_s = cc.xcorr_max(cc_obspy)

            SW_syn_shift = self.shift(SW_env_syn[i].data, -shift)/(np.mean(amplitude))

            misfit = np.append(misfit,
                               np.matmul((SW_env_obs[i].data - SW_syn_shift).T, (SW_env_obs[i].data - SW_syn_shift )) / (
                                   2 * (var*0.1)))
            time = -shift * dt
            # params = {'legend.fontsize': 'x-large',
            #           'figure.figsize': (15, 15),
            #           'axes.labelsize': 25,
            #           'axes.titlesize': 'x-large',
            #           'xtick.labelsize': 25,
            #           'ytick.labelsize': 25}
            # pylab.rcParams.update(params)
            # plt.figure(figsize=(10, 10))
            # Axes = plt.subplot()
            # time_array = np.arange(len(self.zero_to_nan(SW_env_obs.traces[i]))) * SW_env_obs.traces[i].meta.delta + (SW_env_obs.traces[i].meta.starttime - obspy.UTCDateTime(2020, 1, 2, 3, 4, 5))
            # Axes.plot(time_array,self.zero_to_nan(SW_syn_shift),label='Shifted+Scaled Rayleigh wave (syn) ')
            # Axes.plot(time_array,self.zero_to_nan(SW_env_syn.traces[i].data),label = 'Synthetic Rayleigh wave',c='r')
            # Axes.plot(time_array,self.zero_to_nan(SW_env_obs.traces[i].data),label = 'Observed Rayleigh wave',c='k')
            # Axes.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
            # plt.xlabel(obspy.UTCDateTime(2020, 1, 2, 3, 4, 5).strftime('Time : %Y-%m-%dT%H:%M:%S + [sec]'))
            # plt.ylabel("Displacement in Z-component [m]")
            # plt.legend()
            # plt.tight_layout()
            # plt.show()
            # plt.savefig(self.save_dir + '/Shift.pdf' )
            # plt.close()

        sum_misfit = np.sum(misfit)
        return misfit

    def shift(self, np_array, time_shift):
        new_array = np.zeros_like(np_array)
        if time_shift < 0:
            new_array[-time_shift:] = np_array[:time_shift]
        elif time_shift == 0:
            new_array[:] = np_array[:]
        else:
            new_array[:-time_shift] = np_array[time_shift:]
        return new_array

    def zero_to_nan(self,values):
        """Replace every 0 with 'nan' and return a copy."""
        return [float('nan') if x==0 else x for x in values]

