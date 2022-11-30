The MTM (Multiple Taper Method) toolkit produces estimates of power spectral
density from time series. The multiple taper (also called multiple window) 
method produces good spectral estimates on a variety of time series, including
those with wide dynamic range, short length, and high signal-to-noise ratio.

The method, developed by Thomson (1982), applies a set of windows to the time
series to produce several spectral estimates. These estimates can then be
averaged (linearly or adaptively) to yield a spectral estimate with good
resolution and low bias and variance. The windows which are used are based
on discrete prolate spheroidal sequences. These have the properties that
1) they are mutually orthogonal, so that individual spectral estimates are
approximately independent; and 2) the windows are optimally leakage-resistant
in a specified frequency range.

For a hands-off spectral estimate, try [S,c]=mtm(x); where x is the time series
of interest. This will return a spectral estimate S, with confidence limits c.
Try this for a short series first, as the calculation of the dpss can be slow.

This toolkit was written by:
 
  Eric Breitenberger, Geophysical Institute,University of Alaska - Fairbanks.

Comments and suggestions are welcome and can be directed to eric@gi.alaska.edu. 
 
-----------------------------------------------------------------------------
The following functions are included:

MTM   	 - Simple driver routine to call various MTM and DPSS routines.
  
MTMLIN   - Computes a multiple window power spectrum estimate with confidence
           intervals, using linear weighting.   

MTMADAP  - Computes a multiple window power spectrum estimate with confidence
           intervals, using adaptive weighting. 

DPSS     - Driver routine for getting the DPSS by calculation, interpolation,
           or loading from disk. 

DPSSCALC - Calculates the discrete prolate spheroidal sequences using the 
           tridiagonal method of Slepian. DPSSCALC employs spectrum slicing 
           and other tricks to be *much* faster than the previous version
           for large N. The old DPSSCALC has been renamed DPSSCAL2.

DPSSCAL2 - Directly calculates the discrete prolate spheroidal sequences
           using the tridiagonal method of Slepian. This is faster than 
           DPSSCALC for fairly small N (N~<200).

DPSSINT  - Interpolates DPSS of length N from previously calculated sequences
           of a different length.  

DPSSSAVE - Utility routine to calculate DPSS and save them on disk in the form
           used by dpssint.m and dpss.m.

WILHIL   - Wilson-Hilferty approximation to chi-squared distribution. 
           Used by MTMLIN and MTMADAP in computing confidence intervals.

TRISLICE - Spectrum slicing for a tridiagonal matrix. Used by DPSSCALC.
------------------------------------------------------------------------------
Version history:

First released 10/3/95.
2/6/96: Added a much faster routine for the calculation of dpss for large N.

------------------------------------------------------------------------------

References:

Chave, A. D., D. J. Thomson, and M. E. Ander,1987, On the robust estimation of 
    power spectra, coherences, and transfer functions, Journal of Geophysical
    Research 92 (B1), 633-648.

Park, J., C. R. Lindberg, and F. L. Vernon III, 1987, Multitaper spectral analysis
    of high-frequency seismograms, Journal of Geophysical Research 92, 12,675-12,684.

Percival, D. B., and A. T. Walden, 1993, Spectral Analysis for Physical Applications:
    Multitaper and Conventional Univariate Techniques, Cambridge University Press, 
    Cambridge, 583 pp.

Thomson, D. J., 1982, Spectrum estimation and harmonic analysis, Proceedings of the
    IEEE 70, 1055-1096.

Thomson, D. J., 1995, The seasons, global temperature, and precession,
    Science 268, 59-68.

Thomson, D. J., and A. D. Chave, 1991, Jackknifed error estimates for spectra, 
    coherences, and transfer functions, Advances in Spectrum Analysis and Array
    Processing, Vol. 1, ed. by S. Haykin, Prentice-Hall, 58-113.

Vernon, F. L., J. Fletcher, L. Carroll, A. Chave, and E. Sembera, 1991, Coherence
    of seismic body waves from local events as measured by a small-aperture array,
    Journal of Geophysical Research, 96 (B7), 11981-11996.
