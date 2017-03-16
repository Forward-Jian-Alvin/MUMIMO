function [Yc, Hc, N0, EsN0dB, sgmtSNR, NrHc] = fx_soraMIMO_Nt4(X, Nr, cfg)

    Nt = 4;
    Tpkt = 152;
    Tstbc = 4;
    Es = X(1:end) * X(1:end)' / numel(X);
    
    Tc = Tpkt / Tstbc;

    Ncw  = size(X, 2);
    Nc = ceil(Ncw / Tc);
    Yc = cell(Nc, 1);
    Hc = cell(Nc, 1);

    NrHc = zeros(Nc,1);
    Noise = -inf(4, Nc*Tpkt);
    sgmtSNR = zeros(Nc,1);
    for ii = 1:Nc
        [H, N] = fx_soraMIMOchannelLoad(cfg.localCSIpath, cfg.SORAscen, cfg.maxNr);
        scale = 1./sqrt(mean(abs(H).^2, 2));
        H = diag(scale) * H;
        N = diag(scale) * N;
        
        NrHc(ii) = min(Nr, size(H,1));
        Hc{ii} = H(1:NrHc(ii),:);
        
        Noise(1:size(H,1),(ii-1)*Tpkt+(1:Tpkt)) = sqrt(Es)*N;
        
        sgmtSNR(ii) = 10*log10(Nt/mean(abs(N(:)).^2));
    end
    
    ind = find(Noise ~= -inf);
    N0 = mean(abs(Noise(ind)).^2);
    EsN0dB = 10*log10(Nt*Es / N0);
    N0 = N0 / 2; % 2 coefficients are mapped to one complex symbol
    
    rand('seed', 1);
    Noise(ind(randperm(numel(ind)))) = Noise(ind);
%     Noise = reshape(Noise(1:Nr,:), 4*Nr, []);
    
    offset = 0;
    for ii = 1:Nc
        tmp = reshape(Noise(1:NrHc(ii),(ii-1)*Tpkt+(1:Tpkt)), 4*NrHc(ii), []);
        if ii*Tc <= Ncw
            Xc = X(:,offset+(1:Tc));
            N = tmp;
        else
            Xc = X(:,offset+1:Ncw);
            N = tmp(:,1:Ncw-offset);
        end
        offset = offset + Tc;
        % 4 time slots code word
        Yc{ii} = blkdiag(Hc{ii}, Hc{ii}, Hc{ii}, Hc{ii}) * Xc + N;
    end
end