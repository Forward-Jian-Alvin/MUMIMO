function [Yc, Hc, N0] = fx_MIMO_Nt4(X, Tc, EsN0dB, Nr)

    Nt = 4;
    Es = X(1:end) * X(1:end)' / numel(X);

    Ncw  = size(X, 2);
    Nc = ceil(Ncw / Tc);
    Yc = cell(Nc, 1);
    Hc = cell(Nc, 1);

    offset = 0;
    for ii = 1:Nc
        if ii*Tc <= Ncw
            Xc = X(:,offset+(1:Tc));
        else
            Xc = X(:,offset+1:end);
        end
        offset = offset + Tc;
        % Rayleigh fading
        H = fx_MIMOchannelLoad(Nr);
        Hc{ii} = H;
        % 4 time slots code word
        Yc{ii} = blkdiag(H, H, H, H) * Xc + fx_AWGN(Nt*Es, EsN0dB, [Nr/Nt*size(Xc,1) size(Xc,2)]);
    end
    
    N0 = Nt * Es * 10^(-EsN0dB/10) / 2; % 2 coefficients are mapped to one complex symbol
end