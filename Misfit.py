import numpy as np
import obspy.signal.cross_correlation as cc

import matplotlib.pylab as plt
from obspy.core.stream import Stream

class Misfit:
    def CC_BW(self,BW_obs,BW_syn,or_time,plot = False):
        # fig = plt.figure(figsize=(12, 10))
        # # ax1 = plt.subplot2grid((2, 1), (0, 0))
        # ax1 = plt.subplot(511)
        # ymin, ymax = ax1.get_ylim()
        # ax1.plot(BW_obs.P_stream.traces[0], alpha = 0.5,c = 'r', label = 'Observed')
        # ax1.plot(BW_syn.P_stream.traces[0], alpha = 0.5,c= 'g', label = 'Synthetic')
        # ax1.tick_params(axis='x', labelsize=20)
        # ax1.tick_params(axis='y', labelsize=20)
        # plt.legend(loc='upper right')
        #
        # # ax2 = plt.subplot2grid((2, 1), (1, 0))
        # ax2 = plt.subplot(512)
        # ax2.plot(BW_obs.P_stream.traces[1], alpha = 0.5,c = 'r', label = 'Observed')
        # ax2.plot(BW_syn.P_stream.traces[1], alpha = 0.5,c= 'g', label = 'Synthetic')
        # ax2.tick_params(axis='x', labelsize=20)
        # ax2.tick_params(axis='y', labelsize=20)
        # plt.legend(loc ='upper right')
        # plt.xlabel(or_time.strftime('Time : %Y-%m-%dT%H:%M:%S + [sec]'),fontsize = 20)
        #
        # ax3 = plt.subplot(513)
        # ymin, ymax = ax3.get_ylim()
        # ax3.plot(BW_obs.S_stream.traces[0], alpha = 0.5,c = 'r', label = 'Observed')
        # ax3.plot(BW_syn.S_stream.traces[0], alpha = 0.5,c= 'g', label = 'Synthetic')
        # ax3.tick_params(axis='x', labelsize=20)
        # ax3.tick_params(axis='y', labelsize=20)
        # plt.legend(loc='upper right')
        #
        # ax4 = plt.subplot(514)
        # ymin, ymax = ax1.get_ylim()
        # ax4.plot(BW_obs.S_stream.traces[1], alpha = 0.5,c = 'r', label = 'Observed')
        # ax4.plot(BW_syn.S_stream.traces[1], alpha = 0.5,c= 'g', label = 'Synthetic')
        # ax4.tick_params(axis='x', labelsize=20)
        # ax4.tick_params(axis='y', labelsize=20)
        # plt.legend(loc='upper right')
        #
        # ax5 = plt.subplot(515)
        # ymin, ymax = ax1.get_ylim()
        # ax5.plot(BW_obs.S_stream.traces[2], alpha = 0.5,c = 'r', label = 'Observed')
        # ax5.plot(BW_syn.S_stream.traces[2], alpha = 0.5,c= 'g', label = 'Synthetic')
        # ax5.tick_params(axis='x', labelsize=20)
        # ax5.tick_params(axis='y', labelsize=20)
        # plt.legend(loc='upper right')
        #
        #
        # plt.tight_layout()
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


        if plot:
            fig = plt.figure(figsize=(10, 12))
        else:
            fig = 1

        # S - correlations:
        # Calculate Shift based on T component
        len_S_obs = len(s_obs[2].data)

        cc_obspy = cc.correlate(s_syn[2].data,s_obs[2].data,int(len_S_obs))
        shift_centered, MAX_CC = cc.xcorr_max(cc_obspy, abs_max=False)

        shift = np.argmax(cc_obspy)
        time_shift = np.append(time_shift, shift_centered)

        mu_s = np.array([0.7, 0.7, 0.9]) # S_T can be shifted to 0.95 at some point
        sigma_s = np.array([0.3, 0.3, 0.1])


        for i in range(len(s_obs)):
            delta = s_obs[i].stats.delta
            cc_obspy = cc.correlate(s_syn[i].data, s_obs[i].data, int(len(s_obs[i].data)))
            CC_s = cc_obspy[shift]

            s_syn_shift_obspy = self.shift(s_syn[i].data, shift_centered)

            misfit = np.append(misfit, ((CC_s - mu_s[i]) ** 2) / (2 * (sigma_s[i]) ** 2))# + np.abs(shift))

            if plot:
                time_array = np.arange(len(s_obs[i].data)) * delta
                start = int((s_start_obs.timestamp - or_time.timestamp - 10) / delta)
                end =int((s_start_obs.timestamp - or_time.timestamp + 30) / delta)

                ax1 = plt.subplot2grid((5, 1), (i+2, 0))

                plt.plot(time_array[start:end],self.normalize(s_obs[i][start:end]),'b', label = 'Observed')
                # plt.plot(time_array[start:end],self.normalize(s_syn[i][start:end]),'r', linewidth = 0.3, label = 'Synthetic')
                plt.plot(time_array[start:end],self.normalize(s_syn_shift_obspy[start:end]),'g', label = 'Shifted')
                ymin, ymax = ax1.get_ylim()
                xmin, xmax = ax1.get_xlim()
                # if i == 0:
                #     plt.text(xmax - 30, ymax / 1.7, "S-Z - %.4f, %.4f " % (misfit[i] , CC_s), fontsize=20, color='b')
                # elif i == 1:
                #     plt.text(xmax - 30, ymax / 1.7, "S-R - %.4f, %.4f " % (misfit[i] , CC_s), fontsize=20, color='b')
                # elif i == 2:
                #     plt.text(xmax - 30, ymax / 1.7, "S-T - %.4f, %.4f " % (misfit[i] , CC_s), fontsize=20, color='b')
                ax1.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
                ax1.tick_params(axis='x', labelsize=18)
                ax1.tick_params(axis='y', labelsize=18)
                plt.tight_layout()


        ## P - correlation based on Z component
        len_P_obs = len(p_obs[0].data)
        cc_obspy = cc.correlate(p_syn[0].data[0:len_P_obs],p_obs[0].data[0:len_P_obs],int(len_P_obs))
        shift_centered, MAX_CC = cc.xcorr_max(cc_obspy, abs_max=False)

        shift = np.argmax(cc_obspy)
        time_shift = np.append(time_shift, shift_centered)

        mu_p = np.array([0.95, 0.9])
        sigma_p = np.array([0.1, 0.2])



        # P- correlation
        for i in range(len(p_obs)):
            len_P_obs = len(p_obs[i].data)

            cc_obspy = cc.correlate(p_syn[i].data,p_syn[i].data,int(len_P_obs))

            CC_p = cc_obspy[shift]

            misfit = np.append(misfit, ((CC_p - mu_p[i] ) ** 2) / (2 * (sigma_p[i]) ** 2) )#+ np.abs(shift))

            p_syn_shift_obspy = self.shift(p_syn[i].data, shift_centered)
            start = int((p_start_obs.timestamp - or_time.timestamp - 10) / delta)
            end = int((p_start_obs.timestamp - or_time.timestamp + 30) / delta)

            A = (np.dot(p_obs[i].data, p_syn_shift_obspy) / np.dot(p_obs[i].data, p_obs[i].data))
            amplitude = np.append(amplitude,abs(A))

            if plot:
                delta = p_obs[i].stats.delta
                time_array = np.arange(len_P_obs) * delta

                ax1 = plt.subplot2grid((5, 1), (i, 0))
                plt.plot(time_array[start:end],self.normalize(p_obs[i][start:end]),'b', label = 'Observed')
                # plt.plot(time_array[start:end],self.normalize(p_syn[i][start:end]),'r', linewidth = 0.3, label = 'Synthetic')
                plt.plot(time_array[start:end],self.normalize(p_syn_shift_obspy[start:end]),'g', label = 'Shifted')
                ymin, ymax = ax1.get_ylim()
                xmin, xmax = ax1.get_xlim()
                # if i == 0:
                #     plt.text(xmax - 30, ymax / 1.7, "P-Z - %.4f - %.4f " % (misfit[i+3] , CC_p), fontsize=20, color='b')
                # elif i == 1:
                #     plt.text(xmax - 30, ymax / 1.7, "P-R - %.4f - %.4f " % (misfit[i+3] , CC_p), fontsize=20, color='b')
                ax1.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
                ax1.tick_params(axis='x', labelsize=18)
                ax1.tick_params(axis='y', labelsize=18)
                plt.tight_layout()


            if i == 0 and plot == True:
                plt.title('Depth: 90000 (m)')
                plt.legend(loc='lower left', fontsize=15)
        # plt.show()
        sum_misfit = np.sum(misfit)
        return misfit, amplitude, time_shift,fig


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

    def normalize(self, v):
        """This function fails when division by zero"""
        normalized_v = v / np.sqrt(np.sum(v ** 2))
        return normalized_v


