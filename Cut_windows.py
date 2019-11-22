from obspy.taup import TauPyModel
from obspy.core.stream import Stream
from obspy.core.trace import Trace
import numpy as np
import obspy

class Cut_windows:
    def __init__(self, veloc_model_taup, P_HP, P_LP, S_HP, S_LP,Pre_P,Pre_S,Post_P,Post_S, zero_phase = True,Order = 4, Taper = True):
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

    def get_P(self, epi, depth_m):
        model = TauPyModel(model=self.veloc_model)
        tt = model.get_travel_times(source_depth_in_km=depth_m / 1000, distance_in_degree=epi,phase_list=['P'])
        return tt[0].time

    def get_S(self, epi, depth_m):
        model = TauPyModel(model=self.veloc_model)
        tt = model.get_travel_times(source_depth_in_km=depth_m / 1000, distance_in_degree=epi,phase_list=['S'])
        return tt[0].time

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
            self.P_original.taper(0.05, 'hann')
            self.S_original.taper(0.05, 'hann')

        if self.zero_phase:
            self.P_original.filter('highpass', freq=1. / (end_P - self.start_P), zerophase=True, corners=self.Order)
            self.P_original.filter('highpass', freq=1. / (end_S - self.start_S), zerophase=True, corners=self.Order)
            self.S_original.filter('highpass', freq=1. / (end_P - self.start_P), zerophase=True, corners=self.Order)
            self.S_original.filter('highpass', freq=1. / (end_S - self.start_S), zerophase=True, corners=self.Order)
        else:
            self.P_original.filter('highpass', freq=1. / (end_P - self.start_P) , corners=self.Order)
            self.P_original.filter('highpass', freq=1. / (end_S - self.start_S) , corners=self.Order)
            self.S_original.filter('highpass', freq=1. / (end_P - self.start_P), corners=self.Order)
            self.S_original.filter('highpass', freq=1. / (end_S - self.start_S), corners=self.Order)

        self.P_original = self.Filter(self.P_original, HP=self.P_HP, LP=self.P_LP)
        self.S_original = self.Filter(self.S_original, HP=self.S_HP, LP=self.S_LP)

        P_stream = Stream()
        S_stream = Stream()
        for i in range(0,len(stream.traces)):
            trace_P = self.P_original.traces[i].copy()
            trace_S = self.S_original.traces[i].copy()
            dt = trace_P.meta.delta


            P_trace = Trace.slice(trace_P, self.start_P, end_P)
            self.P_len = len(P_trace)
            npts_p = self.P_len + 2 * 500
            start_p = dt * 500
            S_trace = Trace.slice(trace_S, self.start_S, end_S)
            self.S_len = len(S_trace)
            npts_s = self.S_len + 2 * 500
            start_s = dt * 500

            if i == 2:
                total_s_trace = Trace(np.zeros(npts_s),
                                      header={"starttime": self.start_S - start_s, 'delta': trace_S.stats.delta,
                                              "station": trace_S.stats.station,
                                              "network": trace_S.stats.network, "location": trace_S.stats.location,
                                              "channel": trace_S.stats.channel}).__add__(S_trace, method=0,
                                                                                       interpolation_samples=0,
                                                                                       fill_value=S_trace.data,
                                                                                       sanity_checks=True)

            else:
                total_p_trace = Trace(np.zeros(npts_p),
                                      header={"starttime": self.start_P - start_p, 'delta': trace_P.stats.delta,
                                              "station": trace_P.stats.station,
                                              "network": trace_P.stats.network, "location": trace_P.stats.location,
                                              "channel": trace_P.stats.channel}).__add__(P_trace, method=0,
                                                                                       interpolation_samples=0,
                                                                                       fill_value=P_trace.data,
                                                                                       sanity_checks=True)
                total_s_trace = Trace(np.zeros(npts_s),
                                      header={"starttime": self.start_S - start_s, 'delta': trace_S.stats.delta,
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


    def Filter(self, stream, HP, LP):
        if self.zero_phase:
            stream.filter('highpass', freq= HP, zerophase=True, corners = self.Order)
            stream.filter('lowpass' , freq=LP, zerophase=True, corners = self.Order)
        else:
            stream.filter('highpass', freq= HP, corners = self.Order)
            stream.filter('lowpass' , freq=LP, corners = self.Order)
        return stream








