# znorm_cwt
Z-normalization of CWT

This function enhances visual inspection of continuous wavelet transform in brain activity, by rectifying the 1/f spectral behavior of the EEG. First, a Gaussian distribution is fitted on the Morlet coefficients at each scale of the CWT. Then, the normalized, z-scaled wavelet power spectrum was defined as the absolute value-squared of the z-scored wavelet coefficients. This method was first described by Roehri et al. and is similar to the implementation used in Bernardo et al.

Example of usage is shown here: https://github.com/dbernardo05/TF_CNN_IEEG/blob/master/prep_data_v1_20170629.ipynb

Bernardo, Danilo, et al. "Visual and semi-automatic non-invasive detection of interictal fast ripples: A potential biomarker of epilepsy in children with tuberous sclerosis complex." Clinical Neurophysiology 129.7 (2018): 1458-1466.

Roehri N, Lina J-M, Mosher JC, Bartolomei F, Benar C-G. Time-frequency strategies for increasing high-frequency oscillation detectability in intracerebral EEG. IEEE Trans Biomed Eng 2016;63:2595–606.
