function [Y,sigma_n] = fx_AWGN(X, EsN0dB)

Es =  X(1:end) * X(1:end)'/numel(X);

N0 = Es * 10^(-EsN0dB/10);
noise = sqrt(N0/2) * (randn(X) + 1i * randn(X));
sigma_n = 10^(-EsN0dB/20);

Y = X + noise;
end