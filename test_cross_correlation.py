import numpy as np
import obspy.signal.cross_correlation as cc

import matplotlib.pylab as plt


def SHIFT(np_array, time_shift):
    new_array = np.zeros_like(np_array)
    if time_shift < 0:
        new_array[-time_shift:] = np_array[:time_shift]
    elif time_shift == 0:
        new_array[:] = np_array[:]
    else:
        new_array[:-time_shift] = np_array[time_shift:]
    return new_array

a = np.ones(5)
trace_a = np.hstack((np.zeros(100),a,np.zeros(100)))
b = np.ones(5)
trace_b = np.hstack((np.zeros(100),b,np.zeros(100)))
trace_b = np.roll(trace_b, 5)

cc_obspy = cc.correlate(trace_b,trace_a,5)
shift_centered, MAX_CC = cc.xcorr_max(cc_obspy, abs_max=False)

shift = np.argmax(cc_obspy)
b_shift = SHIFT(trace_b, shift_centered)

plt.figure()
plt.subplot(111)
plt.plot(trace_a, label = 'Observed')
plt.plot(trace_b, label = 'Synthetic')
plt.plot(b_shift, LineStyle = ':', label = 'Synthetic shifted')
plt.legend()
plt.show()

a=1

