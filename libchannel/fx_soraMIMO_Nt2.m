function [Yc, Hc, N0, EsN0dB, sgmtSNR] = fx_soraMIMO_Nt2(X, Nr, cfg)

    Nt = 2;
    Tpkt = 152;
    Tstbc = 2;
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
        
        if size(H,1) == 4
            NrHc(ii) = Nr;
            indNr = Nr;
        else
            NrHc(ii) = min(Nr, size(H,1));
            indNr = 0;
            %Hc{ii} = H(1:NrHc(ii),1:Nt);
        end
        Hc{ii} = H(indNr+(1:NrHc(ii)),1:Nt);
        scale2 = 1./sqrt(mean(abs(Hc{ii}).^2, 2));
        Hc{ii} = diag(scale2) * Hc{ii};
        
        Noise(1:size(H,1),(ii-1)*Tpkt+(1:Tpkt)) = sqrt(Es)*N;
        
        sgmtSNR(ii) = 10*log10(Nt/mean(abs(N(:)).^2));
%         Hc{ii} = H(1:Nr,1:Nt);
%         Noise = [Noise sqrt(Es)*N];
%         sgmtSNR = [sgmtSNR 10*log10(Nt/mean(abs(N(:)).^2))];
    end
    
    ind = find(Noise ~= -inf);
    N0 = mean(abs(Noise(ind)).^2);
%     N0 = mean(abs(Noise(:)).^2);
    EsN0dB = 10*log10(Nt*Es / N0);
    N0 = N0 / 2; % 2 coefficients are mapped to one complex symbol
    
    rand('seed', 1);
    Noise(ind(randperm(numel(ind)))) = Noise(ind);
%     Noise(randperm(numel(Noise))) = Noise(1:end);
%     Noise = reshape(Noise(1:Nr,:), 2*Nr, []);
    
    offset = 0;
    for ii = 1:Nc
        tmp = reshape(Noise(indNr+(1:NrHc(ii)),(ii-1)*Tpkt+(1:Tpkt)), 2*NrHc(ii), []);
        if ii*Tc <= Ncw
            Xc = X(:,offset+(1:Tc));
            N = tmp;
%             N = Noise(:,offset+(1:Tc));
        else
            Xc = X(:,offset+1:Ncw);
            N = tmp(:,1:Ncw-offset);
%             N = Noise(:,offset+1:Ncw);
        end
        offset = offset + Tc;
        % 4 time slots code word
        Yc{ii} = blkdiag(Hc{ii}, Hc{ii}) * Xc + N;
    end
end