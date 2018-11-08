import numpy as np

class Model_samples:
    def __init__(self,prior):
        self.prior = prior
    def model_samples(self, update = None,epi_old=None, depth_old=None):
        if update == None:
            epi_sample = np.random.uniform(self.prior['epi']['range_min'], self.prior['epi']['range_max'])
            depth_sample = np.around(
                np.random.uniform(self.prior['depth']['range_min'], self.prior['depth']['range_max']), decimals=1)
        else:
            if update == 'epi':
                epi_sample = np.random.normal(epi_old, self.prior['epi']['spread'])
                depth_sample = depth_old
            elif update == 'depth':
                epi_sample = epi_old
                depth_sample = np.around(np.random.normal(depth_old, self.prior['depth']['spread']), decimals=1)
            else:
                epi_sample = epi_old
                depth_sample = depth_old
        return epi_sample, depth_sample

    def model_samples_sdr(self, update = None,strike_old=None, dip_old=None, rake_old=None):
        if update == None:
            strike = np.random.uniform(self.prior['strike']['range_min'], self.prior['strike']['range_max'])
            dip = np.random.uniform(self.prior['dip']['range_min'], self.prior['dip']['range_max'])
            rake = np.random.uniform(self.prior['rake']['range_min'], self.prior['rake']['range_max'])

        else:
            if update == 'moment':
                # Change radian to degree
                dip_rad = np.deg2rad(dip_old)
                strike_rad = np.deg2rad(strike_old)

                # Calculate normal vector of Fault geometry using strike and dip:
                n = np.array(
                    [-np.sin(dip_rad) * np.sin(strike_rad), -np.sin(dip_rad) * np.cos(strike_rad), np.cos(dip_rad)])

                # X,Y,Z coordinate of the Northpole:
                north_coor = np.array([0, 0, 1])

                # Rotation Axis of from Northpole to Old_sample
                R = np.cross(north_coor, n)
                R_norm = R / (np.sqrt(np.sum(np.square(R), axis=0)))

                # New proposal angle which depends on a spread specified:
                random_phi = np.abs(np.random.normal(0, self.prior['angle_spread']))
                phi = np.deg2rad(random_phi)

                # Theta will be choosen from a point on a circle all with epicentral radius of: Phi
                random_theta = np.random.choice(
                    np.around(np.linspace(0, 360, 361),
                              decimals=1))
                theta = np.deg2rad(random_theta)

                # X,Y,Z coordinates of the new_sample, BUT looking from the northpole (SO STILL NEEDS ROTATION)
                new_coor = np.array([-np.sin(phi) * np.sin(theta), -np.sin(phi) * np.cos(theta), np.cos(phi)])

                # X,Y,Z coordinates with rotation included --> using: Rodrigues' Rotation Formula
                beta = dip_rad
                R_new_coor = np.cos(beta) * new_coor + np.sin(beta) * (np.cross(R_norm, new_coor)) + (np.inner(R_norm,
                                                                                                               new_coor)) * (
                                                                                                         1 - np.cos(
                                                                                                             beta)) * R_norm

                # Determine from a normal distribution a new rake: mean = old rake and SD = spread
                rake = np.random.normal(rake_old, self.prior['rake']['spread'])

                if rake < -180:
                    rake = (180 + rake) % 180
                if rake == 180:
                    rake = -rake
                if rake > 179:
                    rake = (180 + rake) % -180

                phi_normal = np.rad2deg(np.arctan2(R_new_coor[1], R_new_coor[0]))
                wrap = phi_normal % 360
                st = 360.0 - 90.0 - wrap
                strike = st % 360.0

                dip = np.rad2deg(np.arctan2(np.sqrt(R_new_coor[0] ** 2 + R_new_coor[1] ** 2), R_new_coor[2]))
                # dip = 90.0
                if dip >= 90 or dip < 0:
                    if dip == 90:
                        dip = 89.0
                        strike = (strike + 180.0) % 360.0
                        if rake == -180:
                            pass
                        else:
                            rake = -rake
                    else:
                        d_new = 90.0 - dip
                        dip = d_new % 90.0
                        strike = (strike + 180.0) % 360.0
                        if rake == -180:
                            pass
                        else:
                            rake = -rake
            else:
                strike = strike_old
                dip = dip_old
                rake = rake_old
        return strike, dip, rake