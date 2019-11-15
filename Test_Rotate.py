import obspy
import numpy as np

def _dip_azimuth2zne_base_vector(dip, azimuth):
    """
    Helper function converting a vector described with azimuth and dip of unit
    length to a vector in the ZNE (Vertical, North, East) base.

    The definition of azimuth and dip is according to the SEED reference
    manual.
    """
    dip = np.deg2rad(dip)
    azimuth = np.deg2rad(azimuth)

    return np.array([-np.sin(dip),
                     np.cos(azimuth) * np.cos(dip),
                     np.sin(azimuth) * np.cos(dip)])

inverse = True
data_1 = np.array([1,1,1,1,1,1,1])
data_2 = np.array([2,2,2,2,2,2,2])
data_3 = np.array([8,8,8,8,8,8,8])

dip_1 = -89.9
dip_2 = 0
dip_3 = 0

azimuth_1 = 285
azimuth_2 = 105.2
azimuth_3 = 345.3

# Define the base vectors of the old base in terms of the new base vectors.
base_vector_1 = _dip_azimuth2zne_base_vector(dip_1, azimuth_1)
base_vector_2 = _dip_azimuth2zne_base_vector(dip_2, azimuth_2)
base_vector_3 = _dip_azimuth2zne_base_vector(dip_3, azimuth_3)

m = np.array([base_vector_1,
              base_vector_2,
              base_vector_3])

if not inverse:
    m = np.linalg.inv(m)

z, n, e = np.dot(m, [data_1, data_2, data_3])

