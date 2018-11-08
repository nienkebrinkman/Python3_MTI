import instaseis
import numpy as np


class Get_Seismogram():
    def __init__(self, PRIOR):
        self.prior = PRIOR
        self.db = instaseis.open_db(PRIOR['VELOC'])

    def get_receiver(self):
        receiver = instaseis.Receiver(latitude=self.prior['la_r'], longitude=self.prior['lo_r'],
                                      network=self.prior['network'], station=self.prior['station'],
                                      location=self.prior['location'])
        return receiver

    def get_source(self, la_s, lo_s, depth, strike, dip, rake, time, M0):
        source = instaseis.Source.from_strike_dip_rake(latitude=la_s, longitude=lo_s,
                                                       depth_in_m=depth,
                                                       strike=strike, dip=dip,
                                                       rake=rake, M0=M0,
                                                       origin_time=time)
        return source

    def get_seis_manual(self, la_s, lo_s, depth, strike, dip, rake, time, M0):
        source = self.get_source(la_s=la_s, lo_s=lo_s, depth=depth, strike=strike, dip=dip, rake=rake, time=time, M0=M0)
        receiver = self.get_receiver()
        traces = self.db.get_seismograms(source=source, receiver=receiver, components=self.prior['components'],kind=self.prior['kind'])
        traces.interpolate(sampling_rate=self.prior['sampling_rate'])
        traces.traces[0].data = np.float64(traces.traces[0].data)
        traces.traces[1].data = np.float64(traces.traces[1].data)
        traces.traces[2].data = np.float64(traces.traces[2].data)
        return traces
