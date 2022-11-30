function freqs = fft_freqs( N, rate )

freqs = rate*(-floor(N/2):ceil(N/2-1))/N;
