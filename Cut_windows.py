from obspy.taup import TauPyModel
from obspy.core.stream import Stream
from obspy.core.trace import Trace
import numpy as np
import obspy

class Cut_windows:
    def __init__(self, veloc_model_taup, P_HP, P_LP, S_HP, S_LP,Pre_P,Pre_S,Post_P,Post_S):
        self.veloc_model = veloc_model_taup
        self.P_HP = P_HP
        self.P_LP = P_LP
        self.S_HP = S_HP
        self.S_LP = S_LP
        self.Pre_P = Pre_P
        self.Pre_S = Pre_S
        self.Post_P = Post_P
        self.Post_S = Post_S

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

    def Get_bw_windows(self, stream, epi, depth, or_time):
        self.original = stream
        or_time_sec = or_time.timestamp
        tt_P = self.get_P(epi, depth)
        tt_S = self.get_S(epi, depth)

        self.start_P = obspy.UTCDateTime(or_time_sec + tt_P - self.Pre_P)
        self.or_P_len = int((self.start_P - or_time)/ stream.traces[0].stats.delta)
        self.start_S = obspy.UTCDateTime(or_time_sec + tt_S - self.Pre_S)
        self.or_S_len = int((self.start_S - or_time) / stream.traces[0].stats.delta)


        end_P = obspy.UTCDateTime(or_time_sec + tt_P + self.Post_P)
        end_S = obspy.UTCDateTime(or_time_sec + tt_S + self.Post_S)


        P_stream = Stream()
        S_stream = Stream()
        # BW_stream = Stream()

        for i, trace in enumerate(stream.traces):
            dt = trace.meta.delta
            trace.filter('highpass', freq=1. / (end_P - self.start_P), zerophase=True)
            trace.filter('highpass', freq=1. / (end_S - self.start_S), zerophase=True)

            # Apply a filter with the length of the window

            P_trace = Trace.slice(trace, self.start_P, end_P)
            self.P_len = len(P_trace)
            npts_p = self.P_len + 2 * 500
            start_p = dt * 500
            S_trace = Trace.slice(trace, self.start_S, end_S)
            self.S_len = len(S_trace)
            npts_s = self.S_len + 2 * 500
            start_s = dt * 500

            # === Taper the data ===
            S_trace.taper(0.05, 'hann')
            P_trace.taper(0.05, 'hann')

            if 'T' in trace.stats.channel:
                total_s_trace = Trace(np.zeros(npts_s),
                                      header={"starttime": self.start_S - start_s, 'delta': trace.stats.delta,
                                              "station": trace.stats.station,
                                              "network": trace.stats.network, "location": trace.stats.location,
                                              "channel": trace.stats.channel}).__add__(S_trace, method=0,
                                                                                       interpolation_samples=0,
                                                                                       fill_value=S_trace.data,
                                                                                       sanity_checks=True)

            else:
                total_p_trace = Trace(np.zeros(npts_p),
                                      header={"starttime": self.start_P - start_p, 'delta': trace.stats.delta,
                                              "station": trace.stats.station,
                                              "network": trace.stats.network, "location": trace.stats.location,
                                              "channel": trace.stats.channel}).__add__(P_trace, method=0,
                                                                                       interpolation_samples=0,
                                                                                       fill_value=P_trace.data,
                                                                                       sanity_checks=True)
                total_s_trace = Trace(np.zeros(npts_s),
                                      header={"starttime": self.start_S - start_s, 'delta': trace.stats.delta,
                                              "station": trace.stats.station,
                                              "network": trace.stats.network, "location": trace.stats.location,
                                              "channel": trace.stats.channel}).__add__(S_trace, method=0,
                                                                                       interpolation_samples=0,
                                                                                       fill_value=S_trace.data,
                                                                                       sanity_checks=True)
                P_stream.append(total_p_trace)
            S_stream.append(total_s_trace)
            # === Apply filters ===
            self.S_stream = self.Filter(S_stream, HP=self.S_HP, LP=self.S_LP)
            self.P_stream = self.Filter(P_stream, HP=self.P_HP, LP=self.P_LP)

    def Get_bw_windows_MANUAL(self, stream, tt_P, tt_S, or_time):
        self.original = stream
        self.start_P = obspy.UTCDateTime(tt_P.timestamp - self.Pre_P)
        self.or_P_len = int((self.start_P - or_time) / stream.traces[0].stats.delta)
        self.start_S = obspy.UTCDateTime(tt_S.timestamp - self.Pre_S)
        self.or_S_len = int((self.start_S - or_time) / stream.traces[0].stats.delta)

        self.dt = stream.traces[0].stats.delta

        end_P = obspy.UTCDateTime(tt_P.timestamp + self.Post_P)
        end_S = obspy.UTCDateTime(tt_S.timestamp + self.Post_S)

        P_stream = Stream()
        S_stream = Stream()

        for i, trace in enumerate(stream.traces):
            dt = trace.meta.delta
            trace.filter('highpass', freq=1. / (end_P - self.start_P), zerophase=True)
            trace.filter('highpass', freq=1. / (end_S - self.start_S), zerophase=True)

            # Apply a filter with the length of the window

            P_trace = Trace.slice(trace, self.start_P, end_P)
            self.P_len = len(P_trace)
            npts_p = self.P_len + 2 * 500
            start_p = dt * 500
            S_trace = Trace.slice(trace, self.start_S, end_S)
            self.S_len = len(S_trace)
            npts_s = self.S_len + 2 * 500
            start_s = dt * 500

            # === Taper the data ===
            S_trace.taper(0.05,'hann')
            P_trace.taper(0.05,'hann')

            if 'T' in trace.stats.channel:
                total_s_trace = Trace(np.zeros(npts_s),
                                      header={"starttime": self.start_S - start_s, 'delta': trace.stats.delta,
                                              "station": trace.stats.station,
                                              "network": trace.stats.network, "location": trace.stats.location,
                                              "channel": trace.stats.channel}).__add__(S_trace, method=0,
                                                                                       interpolation_samples=0,
                                                                                       fill_value=S_trace.data,
                                                                                       sanity_checks=True)

            else:
                total_p_trace = Trace(np.zeros(npts_p),
                                      header={"starttime": self.start_P - start_p, 'delta': trace.stats.delta,
                                              "station": trace.stats.station,
                                              "network": trace.stats.network, "location": trace.stats.location,
                                              "channel": trace.stats.channel}).__add__(P_trace, method=0,
                                                                                       interpolation_samples=0,
                                                                                       fill_value=P_trace.data,
                                                                                       sanity_checks=True)
                total_s_trace = Trace(np.zeros(npts_s),
                                      header={"starttime": self.start_S - start_s, 'delta': trace.stats.delta,
                                              "station": trace.stats.station,
                                              "network": trace.stats.network, "location": trace.stats.location,
                                              "channel": trace.stats.channel}).__add__(S_trace, method=0,
                                                                                       interpolation_samples=0,
                                                                                       fill_value=S_trace.data,
                                                                                       sanity_checks=True)
                P_stream.append(total_p_trace)
            S_stream.append(total_s_trace)

            # === Apply filters ===
            self.S_stream = self.Filter(S_stream, HP=self.S_HP, LP=self.S_LP)
            self.P_stream = self.Filter(P_stream, HP=self.P_HP, LP=self.P_LP)

    def Filter(self, stream, HP, LP):
        stream.filter('highpass', freq= HP, zerophase=True)
        stream.filter('lowpass' , freq=LP, zerophase=True)
        return stream








