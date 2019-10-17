import obspy
from obspy.core.trace import Trace
import numpy as np

# Creating trace in obspy
npts = 10000
trace_data = np.random.rand(npts)

or_time = obspy.UTCDateTime(2019,10,15,10,0,0) # Origin Time
dt = 1/20.

Trace = Trace(trace_data, header={"starttime": or_time, 'delta': dt})

start = obspy.UTCDateTime(2019,10,15,10,1,0)
end = obspy.UTCDateTime(2019,10,15,10,1,40)

Trace.trim(starttime=start, endtime=end, pad=True, nearest_sample=True, fill_value=0.0)

a=1