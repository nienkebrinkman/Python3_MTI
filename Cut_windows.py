from obspy.taup import TauPyModel
from obspy.core.stream import Stream
from obspy.core.trace import Trace
import numpy as np
import obspy

class Cut_windows:
    def __init__(self, veloc_model_taup):
        self.veloc_model = veloc_model_taup

    def get_P(self, epi, depth_m):
        model = TauPyModel(model=self.veloc_model)
        tt = model.get_travel_times(source_depth_in_km=depth_m / 1000, distance_in_degree=epi,
                                    phase_list=['P'])

        return tt[0].time

    def get_S(self, epi, depth_m):
        model = TauPyModel(model=self.veloc_model)
        tt = model.get_travel_times(source_depth_in_km=depth_m / 1000, distance_in_degree=epi,
                                    phase_list=['S'])
        return tt[0].time

    def Get_bw_windows(self, stream, epi, depth, or_time,npts):
        self.original = stream
        or_time_sec = or_time.timestamp
        tt_P = self.get_P(epi, depth)
        tt_S = self.get_S(epi, depth)



        self.start_P = obspy.UTCDateTime(or_time_sec + tt_P - 10)
        self.or_P_len = int((self.start_P - or_time)/ stream.traces[0].stats.delta)
        self.start_S = obspy.UTCDateTime(or_time_sec + tt_S - 15)
        self.or_S_len = int((self.start_S - or_time) / stream.traces[0].stats.delta)


        end_P = obspy.UTCDateTime(or_time_sec + tt_P + 20)
        end_S = obspy.UTCDateTime(or_time_sec + tt_S + 50)


        P_stream = Stream()
        S_stream = Stream()
        # BW_stream = Stream()

        for i,trace in enumerate(stream.traces):
            trace.filter('highpass', freq=1.0 / 50.0, zerophase=True)

            P_trace = Trace.slice(trace, self.start_P, end_P)
            self.P_len = len(P_trace)
            npts_p = self.P_len + 2 * self.or_P_len
            S_trace = Trace.slice(trace, self.start_S, end_S)
            self.S_len = len(S_trace)
            npts_s = self.S_len + 2 * self.or_S_len

            zero_trace = Trace(np.zeros(npts),
                               header={"starttime": or_time, 'delta': trace.stats.delta,
                                       "station": trace.stats.station,
                                       "network": trace.stats.network, "location": trace.stats.location,
                                       "channel": trace.stats.channel})

            if 'T' in trace.stats.channel:
                # total_trace = zero_trace.__add__(S_trace, method=0, interpolation_samples=0, fill_value=S_trace.data,
                #                                  sanity_checks=True)
                total_s_trace = Trace(np.zeros(npts_s),
                      header={"starttime": or_time, 'delta': trace.stats.delta,
                              "station": trace.stats.station,
                              "network": trace.stats.network, "location": trace.stats.location,
                              "channel": trace.stats.channel}).__add__(S_trace, method=0, interpolation_samples=0,
                                                                       fill_value=S_trace.data,
                                                                       sanity_checks=True)
            else:
                # P_and_S = P_trace.__add__(S_trace, fill_value=0, sanity_checks=True)
                # total_trace = zero_trace.__add__(P_and_S, method=0, interpolation_samples=0,
                #                                  fill_value=P_and_S.data, sanity_checks=True)
                total_p_trace = Trace(np.zeros(npts_p),
                      header={"starttime": or_time, 'delta': trace.stats.delta,
                              "station": trace.stats.station,
                              "network": trace.stats.network, "location": trace.stats.location,
                              "channel": trace.stats.channel}).__add__(P_trace, method=0, interpolation_samples=0,
                                                                       fill_value=P_trace.data,
                                                                       sanity_checks=True)
                total_s_trace = Trace(np.zeros(npts_s),
                      header={"starttime": or_time, 'delta': trace.stats.delta,
                              "station": trace.stats.station,
                              "network": trace.stats.network, "location": trace.stats.location,
                              "channel": trace.stats.channel}).__add__(S_trace, method=0, interpolation_samples=0,
                                                                       fill_value=S_trace.data,
                                                                       sanity_checks=True)
                P_stream.append(total_p_trace)
            S_stream.append(total_s_trace)
            # BW_stream.append(total_trace)
            # === Apply filters ===
            self.S_stream = self.S_filter(S_stream)
            self.P_stream = self.BW_filter(P_stream)
            # self.BW_stream = self.BW_filter(BW_stream)



    def Get_bw_windows_MANUAL(self, stream, tt_P,tt_S, or_time, npts):
        self.original = stream
        self.start_P = obspy.UTCDateTime(tt_P.timestamp - 10)
        # self.start_P = obspy.UTCDateTime(tt_P.timestamp - 5)
        self.or_P_len = int((self.start_P - or_time)/ stream.traces[0].stats.delta)
        self.start_S = obspy.UTCDateTime(tt_S.timestamp - 15)
        self.or_S_len = int((self.start_S - or_time) / stream.traces[0].stats.delta)

        end_P = obspy.UTCDateTime(tt_P.timestamp + 20)
        end_S = obspy.UTCDateTime(tt_S.timestamp + 50)
        # end_S = obspy.UTCDateTime(tt_S.timestamp + 35)

        P_stream = Stream()
        S_stream = Stream()
        BW_stream = Stream()

        for i,trace in enumerate(stream.traces):
            trace.filter('highpass', freq=1.0 / 50.0, zerophase=True)

            P_trace = Trace.slice(trace, self.start_P, end_P)
            self.P_len = len(P_trace)
            npts_p = self.P_len + 2 * self.or_P_len
            S_trace = Trace.slice(trace, self.start_S, end_S)
            self.S_len = len(S_trace)
            npts_s = self.S_len + 2 * self.or_S_len

            zero_trace = Trace(np.zeros(npts),
                               header={"starttime": or_time, 'delta': trace.stats.delta,
                                       "station": trace.stats.station,
                                       "network": trace.stats.network, "location": trace.stats.location,
                                       "channel": trace.stats.channel})

            if 'T' in trace.stats.channel:
                total_trace = zero_trace.__add__(S_trace, method=0, interpolation_samples=0, fill_value=S_trace.data,
                                                 sanity_checks=True)
                total_s_trace = Trace(np.zeros(npts_s),
                      header={"starttime": or_time, 'delta': trace.stats.delta,
                              "station": trace.stats.station,
                              "network": trace.stats.network, "location": trace.stats.location,
                              "channel": trace.stats.channel}).__add__(S_trace, method=0, interpolation_samples=0,
                                                                       fill_value=S_trace.data,
                                                                       sanity_checks=True)
            else:
                P_and_S = P_trace.__add__(S_trace, fill_value=0, sanity_checks=True)
                total_trace = zero_trace.__add__(P_and_S, method=0, interpolation_samples=0,
                                                 fill_value=P_and_S.data, sanity_checks=True)
                total_p_trace = Trace(np.zeros(npts_p),
                      header={"starttime": or_time, 'delta': trace.stats.delta,
                              "station": trace.stats.station,
                              "network": trace.stats.network, "location": trace.stats.location,
                              "channel": trace.stats.channel}).__add__(P_trace, method=0, interpolation_samples=0,
                                                                       fill_value=P_trace.data,
                                                                       sanity_checks=True)
                total_s_trace = Trace(np.zeros(npts_s),
                      header={"starttime": or_time, 'delta': trace.stats.delta,
                              "station": trace.stats.station,
                              "network": trace.stats.network, "location": trace.stats.location,
                              "channel": trace.stats.channel}).__add__(S_trace, method=0, interpolation_samples=0,
                                                                       fill_value=S_trace.data,
                                                                       sanity_checks=True)
                P_stream.append(total_p_trace)
            S_stream.append(total_s_trace)
            BW_stream.append(total_trace)
            # === Apply filters ===
            self.S_stream = self.S_filter(S_stream)
            self.P_stream = self.BW_filter(P_stream)
            self.BW_stream = self.BW_filter(BW_stream)


    def BW_filter(self, stream):
        # stream.filter('highpass', freq=0.5)
        # stream.filter('lowpass', freq=0.1)


        # stream.filter('highpass', freq=1.0 / 10.0, zerophase=True) #mars
        stream.filter('highpass', freq=1.0/20.0, zerophase=True) #earth
        stream.filter('lowpass', freq=1.0/7.0, zerophase=True)#earth
        # stream.filter('lowpass', freq=1.0/5.0, zerophase=True)#mars
        return stream

    def S_filter(self, stream):
        # stream.filter('highpass', freq=0.5)
        # stream.filter('lowpass', freq=0.1)

        # stream.filter('highpass', freq=1.0 / 30.0)#mars
        stream.filter('highpass', freq=1.0 / 30.0, zerophase=True)#earth
        # stream.filter('highpass', freq=0.05)#mars
        stream.filter('lowpass', freq=1.0 / 7.0, zerophase=True)# earth
        return stream






