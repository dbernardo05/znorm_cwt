

__author__ = "Danilo Bernardo"
__license__ = "MIT"
__version__ = "0.1"

import numpy as np
import sklearn.mixture
import wavelet
from scipy.stats import iqr
from scipy.optimize import curve_fit
from scipy.stats.mstats import mquantiles

def znorm_cwt(wave):
    z_wave = np.zeros_like(wave, dtype=np.complex_)
    for s, scale in enumerate(wave):
        tf_re_og = np.real(scale)
        tf_im_og = np.imag(scale)
         
        iqr2575_re = iqr(tf_re_og)
        iqr2575_im = iqr(tf_im_og)
        iqr_re = mquantiles(tf_re_og)
        iqr_im = mquantiles(tf_im_og)
        tf_re = tf_re_og[np.bitwise_and(tf_re_og>=iqr_re[0]- iqr2575_re*1.5, tf_re_og<=iqr_re[2] + iqr2575_re*1.5)]
        tf_im = tf_im_og[np.bitwise_and(tf_im_og>=iqr_im[0] - iqr2575_im*1.5, tf_im_og<=iqr_im[2] + iqr2575_im*1.5)]

        gmm_re = sklearn.mixture.GMM()
        gmm_im = sklearn.mixture.GMM()

        # Fit gaussian to real
        tf_re_gauss = gmm_re.fit(tf_re[:, np.newaxis]) # GMM requires 2D data as of sklearn version 0.16
        avg_tf_re = tf_re_gauss.means_[0, 0]
        std_tf_re = np.sqrt(tf_re_gauss.covars_[0, 0])

        # Fit gaussian to imaginary
        tf_im_gauss = gmm_im.fit(tf_im[:, np.newaxis]) # GMM requires 2D data as of sklearn version 0.16
        avg_tf_im = tf_im_gauss.means_[0, 0]
        std_tf_im = np.sqrt(tf_im_gauss.covars_[0, 0])        
        
        zero_avg_re = tf_re_og - avg_tf_re
        zero_avg_im = tf_im_og - avg_tf_im
        
        z_re = np.divide(zero_avg_re, std_tf_re)
        z_im = np.divide(zero_avg_im, std_tf_im)

        for n, z_re_n in enumerate(z_re):
            z_wave[s, n] = z_re_n + z_im[n]*1j
         
    power = (np.abs(z_wave))**2
    return power



if __name__ == "__main__":

    # Z-normalize spectra using wavelet co-efficients
    variance = np.std(data_samp)**2
    mean=np.mean(data_samp)
    data_samp = (data_samp - np.mean(data_samp))/np.sqrt(variance)
    
    # Set wavelet parameters
    mother = 'Morlet'
    param = 6 # Default 6
    dt = 0.0005 # Default for 485 Hz = 0.0025
    dj = 0.125/8 #0.03125 # Default = 0.125 // for 485 Hz: range of ~10-80 Hz //
    pad = 0
    s0 = -1
    j1 = -1
    lag1 = 0.72  # lag-1 autocorrelation for red noise background

    # Wavelet transform
    wave,period,scale,coi = wavelet(data_samp,dt,pad,dj,s0,j1,mother);
    z_power = znorm_cwt(wave):
