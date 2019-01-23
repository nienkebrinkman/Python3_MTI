import numpy as np
from obspy.core.stream import Stream


class Get_SNR():
    def getSNR(self, X, TSNR, t1, tracenum=0):
        '''
        Function to calculate SNR of certain part of seismogram relative to
        given time (onset) out of given noise and signal windows. A safety gap
        between noise and signal part can be set. Returns SNR and SNR [dB] and
        noiselevel.
        :param: X, time series (seismogram)
        :type:  `~obspy.core.stream.Stream`
        :param: TSNR, length of time windows [s] around t1 (onset) used to determine SNR
        :type: tuple (T_noise, T_gap, T_signal)
        :param: t1, initial time (onset) from which noise and signal windows are calculated
        :type: float
        '''

        assert isinstance(X, Stream), "%s is not a stream object" % str(X)

        SNR = None
        SNRdB = None
        noiselevel = None

        x = X[tracenum].data
        npts = X[tracenum].stats.npts
        sr = X[tracenum].stats.sampling_rate
        dt = X[tracenum].stats.delta
        t = np.arange(0, npts / sr, dt)

        # get noise window
        inoise = self.getnoisewin(t, t1, TSNR[0], TSNR[1])

        # get signal window
        isignal = self.getsignalwin(t, t1, TSNR[2])
        if np.size(inoise) < 1:
            print("getSNR: Empty array inoise, check noise window!")
            return SNR, SNRdB, noiselevel

        # demean over entire waveform
        x = x - np.mean(x[inoise])

        # calculate ratios
        noiselevel = np.sqrt(np.mean(np.square(x[inoise])))
        # signallevel = np.sqrt(np.mean(np.square(x[isignal])))

        if np.size(isignal) < 1:
            print("getSNR: Empty array isignal, check signal window!")
            return SNR, SNRdB, noiselevel

        # noiselevel = np.abs(x[inoise]).max()
        signallevel = np.abs(x[isignal]).max()

        SNR = signallevel / noiselevel
        SNRdB = 10 * np.log10(SNR)

        return SNR, SNRdB, noiselevel

    def getnoisewin(self,t, t1, tnoise, tgap):
        '''
        Function to extract indeces of data out of time series for noise calculation.
        Returns an array of indeces.
        :param: t, array of time stamps
        :type:  numpy array
        :param: t1, time from which relativ to it noise window is extracted
        :type: float
        :param: tnoise, length of time window [s] for noise part extraction
        :type: float
        :param: tgap, safety gap between t1 (onset) and noise window to
                ensure, that noise window contains no signal
        :type: float
        '''

        # get noise window
        inoise, = np.where((t <= max([t1 - tgap, 0])) \
                           & (t >= max([t1 - tnoise - tgap, 0])))
        if np.size(inoise) < 1:
            inoise, = np.where((t >= t[0]) & (t <= t1))
            if np.size(inoise) < 1:
                print("getnoisewin: Empty array inoise, check noise window!")

        return inoise

    def getsignalwin(self,t, t1, tsignal):
        '''
        Function to extract data out of time series for signal level calculation.
        Returns an array of indeces.
        :param: t, array of time stamps
        :type:  numpy array
        :param: t1, time from which relativ to it signal window is extracted
        :type: float
        :param: tsignal, length of time window [s] for signal level calculation
        :type: float
        '''

        # get signal window
        isignal, = np.where((t <= min([t1 + tsignal, t[-1]])) \
                            & (t >= t1))
        if np.size(isignal) < 1:
            print("getsignalwin: Empty array isignal, check signal window!")

        return isignal


