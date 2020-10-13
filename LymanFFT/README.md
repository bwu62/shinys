# Lyman-break galaxies

### Overview of method

 1. Low-pass filter via Fast Fourier Transform used to smooth spectra
 1. Convolution with template, used to:
   1. correct for red-shift (horizontal translation of log-wavelength)
   1. compute normalized area-under-peak to assess quality-of-match
 1. Partial-quantile rescaling to match flux of spectra at computed red-shift
 1. (Rescaled) Kolmogorov-Smirnov statistic used as second criterion for quality-of-match

### Results

For detailed results, see output file `res.txt`. Below is plot of top 25 scoring spectra.

![Plot of top 25 scoring spectra](https://raw.githubusercontent.com/bwu62/shinys/e5d4f633b2b0c902b368ef9bcb512d08b90eb2aa/LymanFFT/top25.svg)

### To do next

Currently awaiting further feedback on this new method.

