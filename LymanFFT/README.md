# Lyman-break galaxies

## Summary of method

1. Removal of edge artifacts (Gibb's phenomena)
1. Low-pass smoothing via (fast) Fourier Transform
1. Difference with loess-smoothed spectrum used to enhance absorption troughs.
1. Convolution with template, used to:
   1. correct for red-shift (horizontal translation of log-wavelength)
   1. compute normalized area-under-peak to assess quality-of-match
1. Regression of spectra quantiles to match spectra at computed red-shift
1. (Rescaled) Kolmogorov-Smirnov statistic used as second criterion for quality-of-match

## Results

### Initial run of 100 test spectra:

![Plot of top 25 scoring spectra](https://raw.githubusercontent.com/bwu62/shinys/061c9cb70adaf650cb5c30fee44d9f411bb9a139/LymanFFT/top25.svg)

### Actual dataset of 2.459 million spectra:

[Link to results summary page](file:///C:/Users/bi/Desktop/run4/index.html)

## To do next

Currently awaiting further feedback. May investigate better discrimination of quasars by also matching a second template.

