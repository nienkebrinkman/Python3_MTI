from obspy.taup import TauPyModel
from obspy.core.stream import Stream
from obspy.core.trace import Trace
import numpy as np
import obspy

class Cut_windows:
    def __init__(self, veloc_model_taup, P_HP, P_LP, S_HP, S_LP, Pre_P, Pre_S, Post_P, Post_S, global_P_shift, global_S_shift,zero_phase = True, Order = 4, Taper = True, Taper_len = 1.0, Zero_len = 500):
        self.veloc_model = veloc_model_taup
        self.P_HP = P_HP
        self.P_LP = P_LP
        self.S_HP = S_HP
        self.S_LP = S_LP
        self.Pre_P = Pre_P
        self.Pre_S = Pre_S
        self.Post_P = Post_P
        self.Post_S = Post_S
        self.zero_phase = zero_phase
        self.Taper = Taper
        self.Order = Order
        self.Zero_len = Zero_len
        self.Taper_Len = Taper_len # Length of taper that is added to your window
        self.global_P_shift = global_P_shift # Global P shift to the right (so left would be negative) in seconds
        self.global_S_shift = global_S_shift

    def get_P(self, epi, depth_m):
        model = TauPyModel(model=self.veloc_model)
        tt = model.get_travel_times(source_depth_in_km=depth_m / 1000, distance_in_degree=epi,phase_list=['P'])
        tt = tt[0].time + self.global_P_shift
        return tt

    def get_S(self, epi, depth_m):
        model = TauPyModel(model=self.veloc_model)
        tt = model.get_travel_times(source_depth_in_km=depth_m / 1000, distance_in_degree=epi,phase_list=['S'])
        tt = tt[0].time + self.global_S_shift
        return tt

    def get_pp(self, epi, depth_m):
        model = TauPyModel(model=self.veloc_model)
        tt = model.get_travel_times(source_depth_in_km=depth_m / 1000, distance_in_degree=epi,phase_list=['pP'])
        return tt[0].time

    def get_sp(self, epi, depth_m):
        model = TauPyModel(model=self.veloc_model)
        tt = model.get_travel_times(source_depth_in_km=depth_m / 1000, distance_in_degree=epi,phase_list=['sP'])
        return tt[0].time

    def Get_bw_windows(self, stream, UNKNOWN_1, UNKNOWN_2, or_time, MANUAL = False):
        ## YOU can do EITHER MANUAL or not:
        """   if MANUAL = False:
                UNKNOWN_1 = EPI
                UNKNOWN_2 = DEPTH
              if MANUAL = True:
                UNKNOWN_1 = tt_P
                UNKNOWN_2 = tt_S

        """

        self.dt = stream.traces[0].stats.delta
        self.original = stream.copy()
        or_time_sec = or_time.timestamp


        if MANUAL == False:
            epi = UNKNOWN_1
            depth = UNKNOWN_2
            tt_P = self.get_P(epi, depth)
            tt_S = self.get_S(epi, depth)
            self.start_P = obspy.UTCDateTime(or_time_sec + tt_P - self.Pre_P)
            self.or_P_len = int((self.start_P - or_time)/ stream.traces[0].stats.delta)
            self.start_S = obspy.UTCDateTime(or_time_sec + tt_S - self.Pre_S)
            self.or_S_len = int((self.start_S - or_time) / stream.traces[0].stats.delta)
            # print(self.start_P)
            # print(self.start_S)

            end_P = obspy.UTCDateTime(or_time_sec + tt_P + self.Post_P)
            end_S = obspy.UTCDateTime(or_time_sec + tt_S + self.Post_S)

        else:
            tt_P = UNKNOWN_1
            tt_S = UNKNOWN_2
            self.start_P = obspy.UTCDateTime(tt_P.timestamp - self.Pre_P)
            self.or_P_len = int((self.start_P - or_time) / stream.traces[0].stats.delta)
            self.start_S = obspy.UTCDateTime(tt_S.timestamp - self.Pre_S)
            self.or_S_len = int((self.start_S - or_time) / stream.traces[0].stats.delta)

            end_P = obspy.UTCDateTime(tt_P.timestamp + self.Post_P)
            end_S = obspy.UTCDateTime(tt_S.timestamp + self.Post_S)


        self.S_original = self.original.copy()
        self.P_original = self.original.copy()
        if self.Taper == True:
            self.P_original.taper(0.025, 'hann', self.start_P - or_time - 10) # Making sure the P and S wave are not effected
            self.S_original.taper(0.025, 'hann', self.start_P - or_time - 10)

            # CHECK IN OBSPY FUNCTION IF TAPER IS GOOD:
            # import matplotlib.pylab as plt
            # plt.close()
            # x = np.arange(len(taper))
            # plt.plot(x, taper, 'r', label='Hann Taper')
            # plt.plot(3606, 1, 'bx', label='P-arrival')
            # plt.plot(6806, 1, 'kx', label='S-arrival')
            # plt.legend()
            # plt.show()

        if self.zero_phase:
            self.P_original.filter('highpass', freq=1. / (end_P - self.start_P), zerophase=True, corners=self.Order)
            self.P_original.filter('highpass', freq=1. / (end_S - self.start_S), zerophase=True, corners=self.Order)
            self.S_original.filter('highpass', freq=1. / (end_P - self.start_P), zerophase=True, corners=self.Order)
            self.S_original.filter('highpass', freq=1. / (end_S - self.start_S), zerophase=True, corners=self.Order)
        else:
            self.P_original.filter('highpass', freq=1. / (end_P - self.start_P) , corners=self.Order)
            self.P_original.filter('highpass', freq=1. / (end_S - self.start_S) , corners=self.Order)
            self.S_original.filter('highpass', freq=1. / (end_P - self.start_P) , corners=self.Order)
            self.S_original.filter('highpass', freq=1. / (end_S - self.start_S) , corners=self.Order)

        self.P_original = self.Filter(self.P_original, HP=self.P_HP, LP=self.P_LP)
        self.S_original = self.Filter(self.S_original, HP=self.S_HP, LP=self.S_LP)

        wlen_seconds = self.Taper_Len
        zero_len = self.Zero_len
        # wlen = int(wlen_seconds / self.dt)

        P_stream = Stream()
        S_stream = Stream()
        for i in range(0,len(stream.traces)):
            trace_P = self.P_original.traces[i].copy()
            trace_S = self.S_original.traces[i].copy()
            dt = trace_P.meta.delta

            P_trace = Trace.slice(trace_P, self.start_P - wlen_seconds, end_P + wlen_seconds)
            self.P_len = len(P_trace)
            npts_p = self.P_len + 2 * zero_len
            start_p = dt * zero_len
            S_trace = Trace.slice(trace_S, self.start_S - wlen_seconds, end_S + wlen_seconds)
            self.S_len = len(S_trace)
            npts_s = self.S_len + 2 * zero_len
            start_s = dt * zero_len



            if i == 2:
                total_s_trace = Trace(np.zeros(npts_s),
                                      header={"starttime": self.start_S - start_s - wlen_seconds, 'delta': trace_S.stats.delta,
                                              "station": trace_S.stats.station,
                                              "network": trace_S.stats.network, "location": trace_S.stats.location,
                                              "channel": trace_S.stats.channel}).__add__(S_trace, method=0,
                                                                                       interpolation_samples=0,
                                                                                       fill_value=S_trace.data,
                                                                                       sanity_checks=True)

            else:
                total_p_trace = Trace(np.zeros(npts_p),
                                      header={"starttime": self.start_P - start_p - wlen_seconds, 'delta': trace_P.stats.delta,
                                              "station": trace_P.stats.station,
                                              "network": trace_P.stats.network, "location": trace_P.stats.location,
                                              "channel": trace_P.stats.channel}).__add__(P_trace, method=0,
                                                                                       interpolation_samples=0,
                                                                                       fill_value=P_trace.data,
                                                                                       sanity_checks=True)
                total_s_trace = Trace(np.zeros(npts_s),
                                      header={"starttime": self.start_S - start_s - wlen_seconds, 'delta': trace_S.stats.delta,
                                              "station": trace_S.stats.station,
                                              "network": trace_S.stats.network, "location": trace_S.stats.location,
                                              "channel": trace_S.stats.channel}).__add__(S_trace, method=0,
                                                                                       interpolation_samples=0,
                                                                                       fill_value=S_trace.data,
                                                                                       sanity_checks=True)
                P_stream.append(total_p_trace)
            S_stream.append(total_s_trace)
            self.S_stream = S_stream
            self.P_stream = P_stream

            ## Create the taper:
            # Taper_S = self.Create_Taper(self.S_len,wlen, zero_len)
            # Taper_P = self.Create_Taper(self.P_len,wlen,zero_len)
            #
            # # Apply the Taper:
            # if i != 2:
            #     self.P_stream.traces[i].data *= Taper_P
            # self.S_stream.traces[i].data *= Taper_S

            # import matplotlib.pylab as plt
            # plt.figure()
            # plt.subplot(211)
            # plt.plot(self.P_stream.traces[0].data, 'b' , label = 'Non-Tapered P-window')
            # plt.plot(self.P_stream.traces[0].data * Taper_P, 'r' , label = 'Tapered P-window')
            # plt.legend()
            # plt.subplot(212)
            # plt.plot(self.S_stream.traces[0].data, 'b' , label = 'Non-Tapered S-window')
            # plt.plot(self.S_stream.traces[0].data * Taper_S, 'r' , label = 'Tapered S-window')
            # plt.legend()
            # plt.show()
            #
            # a=1



            # Trim the original waveforms after the S arrival so that the waveforms are not too long
            # self.original.trim(endtime=end_S + 25)
            # self.P_original.trim(endtime=end_S + 25)
            # self.S_original.trim(endtime=end_S + 25)

    def Create_Taper(self, Trace_len,wlen, zero_len,):
        hann = np.hanning(wlen * 2)
        # if 2 * wlen == int(wlen_seconds / dt):
        #     hann = np.hanning(wlen * 2)
        # else:
        #     hann = np.hanning(wlen * 2 + 1)
        window = np.hstack((np.zeros(zero_len), hann[:wlen], np.ones(Trace_len - 2 * wlen), hann[len(hann) - wlen:],
                            np.zeros(zero_len)))
        return window

    def Filter(self, stream, HP, LP):
        if self.zero_phase:
            stream.filter('highpass', freq= HP, zerophase=True, corners = self.Order)
            stream.filter('lowpass' , freq=LP, zerophase=True, corners = self.Order)
        else:
            stream.filter('highpass', freq= HP, corners = self.Order)
            stream.filter('lowpass' , freq=LP, corners = self.Order)
        return stream








