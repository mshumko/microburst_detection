# Introduction
This repo contains a few algorithms to detect electron microbursts. The three 
implemented methods are:
- Paul O'Brien's burst parameter (e.g. O'Brien et al., 2003; Douma et al., 2017; Shumko et al., 2020)
- Wavelet-based filtering 
- Signal to background (Blum et al., 2015; Greeley et al., 2019)

All three methods involve some sort of count or standard deviation threshold
and must be tuned to your data. 

After I thought about it for a minute, the Paul O'Brien's burst parameter and 
the signal to background methods are similar if not identical. The only 
difference is over how much time you average over to estimate your background. 

# Dependencies:
Easiest method to install the numpy and pandas depedencies is with pip and 
the requirements.txt file. 
Example: ```pip3 install -r requirements.txt```

# Example Uses:

# Bibliography

Blum, L., X. Li, and M. Denton (2015), Rapid MeV electron precipitation as observed by SAMPEX/HILT during high-speed stream-driven storms, J. Geophys. Res. Space Physics, 120, 3783–3794, doi:10.1002/2014JA020633.

Douma, E., C. J. Rodger, L. W. Blum, and M. A. Clilverd (2017), Occurrence characteristics of relativistic electron microbursts from SAMPEX observations, J. Geophys. Res. Space Physics, 122, doi:10.1002/2017JA024067.

Greeley, A. D., Kanekal, S. G., Baker, D. N., Klecker, B., & Schiller, Q. (2019). Quantifying the contribution of microbursts to global electron loss in the radiation belts. Journal of Geophysical Research: Space Physics, 124. https://doi.org/10.1029/2018JA026368

O’Brien, T. P., K. R. Lorentzen, I. R. Mann, N. P. Meredith, J. B. Blake, J. F. Fennell, M. D. Looper, D. K. Milling, and R. R. Anderson, Energization of relativistic electrons in the presence of ULF power and MeV microbursts: Evidence for dual ULF and VLF acceleration, J. Geophys. Res., 108(A8), 1329, doi:10.1029/2002JA009784, 2003.

Shumko, M., Johnson, A. T., Sample, J. G., Griffith, B. A., Turner, D. L., O'Brien, T. P., et al. (2020). Electron microburst size distribution derived with AeroCube-6. Journal of Geophysical Research: Space Physics, 125, e2019JA027651. https://doi.org/10.1029/2019JA027651