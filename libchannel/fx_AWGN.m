function noise = fx_AWGN(Es, EsN0dB, sz)

N0 = Es * 10^(-EsN0dB/10);
noise = sqrt(N0/2) * (randn(sz) + 1i * randn(sz));