import instaseis
import numpy as np

class Greens_Functions:
    def __init__(self, PRIOR):
        self.prior = PRIOR
        self.db = instaseis.open_db(PRIOR['VELOC'])

    def Get_GF(self, epi, depth):
        gf = self.db.get_greens_function(epicentral_distance_in_degree=epi, source_depth_in_m=depth,
                                          origin_time=self.prior['origin_time'], kind=self.prior['kind'],
                                          dt=1.0 / self.prior['sampling_rate'])

        tss = gf.traces[0].data
        zss = gf.traces[1].data
        rss = gf.traces[2].data
        tds = gf.traces[3].data
        zds = gf.traces[4].data
        rds = gf.traces[5].data
        zdd = gf.traces[6].data
        rdd = gf.traces[7].data
        zep = gf.traces[8].data
        rep = gf.traces[9].data

        G_z = gf.traces[0].meta['npts']
        G_r = gf.traces[0].meta['npts'] * 2
        G_t = gf.traces[0].meta['npts'] * 3
        G = np.ones((G_t, 5))
        G[0:G_z, 0] = zss * (0.5) * np.cos(2 * np.deg2rad(self.prior['az'])) - zdd * 0.5
        G[0:G_z, 1] = - zdd * 0.5 - zss * (0.5) * np.cos(2 * np.deg2rad(self.prior['az']))
        # G[0:G_z, 1] =  zdd * (1/3) + zep * (1/3)
        G[0:G_z, 2] = zss * np.sin(2 * np.deg2rad(self.prior['az']))
        G[0:G_z, 3] = -zds * np.cos(np.deg2rad(self.prior['az']))
        G[0:G_z, 4] = -zds * np.sin(np.deg2rad(self.prior['az']))

        G[G_z:G_r, 0] = rss * (0.5) * np.cos(2 * np.deg2rad(self.prior['az'])) - rdd * 0.5
        G[G_z:G_r, 1] = -0.5 * rdd - rss * (0.5) * np.cos(2 * np.deg2rad(self.prior['az']))
        # G[G_z:G_r, 1] =  rdd * (1/3) + rep * (1/3)
        G[G_z:G_r, 2] = rss * np.sin(2 * np.deg2rad(self.prior['az']))
        G[G_z:G_r, 3] = -rds * np.cos(np.deg2rad(self.prior['az']))
        G[G_z:G_r, 4] = -rds * np.sin(np.deg2rad(self.prior['az']))

        G[G_r:G_t, 0] = -tss * (0.5) * np.sin(2 * np.deg2rad(self.prior['az']))
        G[G_r:G_t, 1] = tss * (0.5) * np.sin(2 * np.deg2rad(self.prior['az']))
        # G[G_r:G_t, 1] =   0
        G[G_r:G_t, 2] = tss * np.cos(2 * np.deg2rad(self.prior['az']))
        G[G_r:G_t, 3] = tds * np.sin(np.deg2rad(self.prior['az']))
        G[G_r:G_t, 4] = -tds * np.cos(np.deg2rad(self.prior['az']))
        return G
