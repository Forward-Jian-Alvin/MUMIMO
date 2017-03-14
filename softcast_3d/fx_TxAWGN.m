function [Y, sigma_n] = fx_TxAWGN(X, EsN0dB)
Es = X(1:end) * X(1:end)'/numel(X);

sigma_n = Es*10^(-EsN0dB/20);

if iscell(X)
    Y = X;
    for ii = 1:numel(X)
        if isempty(X{ii})
            Y{ii} = [];
        else
            noise = sigma_n .* randn(size(X{ii}));
            Y{ii} = X{ii} + noise;
        end
    end
else
    noise = sigma_n .* randn(size(X));
    Y = X + noise;
end