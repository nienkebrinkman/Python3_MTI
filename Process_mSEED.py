from obspy.clients.fdsn.client import Client

class Process_mSEED:
    def __init__(self,stream):
        self.stream = stream
        client = Client("IRIS")
        self.inv = client.get_stations(network=stream.traces[0].stats.network, station=stream.traces[0].stats.station,
                                  level='response')
    def remove_instrument(self):
        # === Remove instrument response ===
        self.stream.remove_response(self.inv,pre_filt=(0.005, 0.01, 0.4, 0.5))
        return self.stream

    def automatic_rotate(self,baz):

        if 'N' in self.stream.traces[0].id or 'N' in self.stream.traces[1].id or 'N' in self.stream.traces[2].id:
            self.stream = self.rotate_NEZ_RTZ(baz)
        elif '1' in self.stream.traces[0].id or '1' in self.stream.traces[1].id or '1' in self.stream.traces[2].id:
            self.stream = self.rotate_Z12_ZRT(baz)
        else:
            raise ValueError('Channels cant be rotated')
        return self.stream


    def rotate_NEZ_RTZ(self,baz):
        # === Convert the NEZ to RTZ coordinates ===
        self.stream.rotate(method='NE->RT', back_azimuth=baz)
        return self.stream

    def rotate_Z12_ZRT(self,baz):
        # === Convert the NEZ to RTZ coordinates ===
        self.stream._rotate_to_zne(self.inv)
        self.stream.rotate(method='NE->RT', back_azimuth=baz)
        return self.stream




