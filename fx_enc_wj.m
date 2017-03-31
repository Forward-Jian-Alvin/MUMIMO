function [ enc,lambda,g ] = fx_enc_wj( C,imginfo )
%FX_ENC_WJ Summary of this function goes here
%   Detailed explanation goes here
%% 划分成64*8行
    x=[];
    for kk = 1:imginfo.GOP
        for tt = 1:8
            for jj=1:8
                currentBlock = C((tt-1)*imginfo.bh+1:tt*imginfo.bh,(jj-1)*imginfo.bw+1:jj*imginfo.bw,kk);
                x = [x;reshape(currentBlock,1,imginfo.bh*imginfo.bw)];
            end
        end
    end
    P1 = 1;% mean(mean(x.*x));
    P = P1*8*8*8;%total power constraint
    lambda = mean((x.*x)');
    lambda = lambda';
    g = sqrt(P/sum(sqrt(lambda)))./sqrt(sqrt(lambda));
    G = diag(g);
    enc=G*x;
end

