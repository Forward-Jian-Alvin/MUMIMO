clear;
filepath='C:\Users\wu\Desktop\Thesixth semester\mimo\3DVSource\3Dtest\3DV≤‚ ‘–Ú¡–\';
filename = sprintf('balloons3.yuv');
fid = fopen([filepath filename],'rb');
TotalFrameNum = 32;
GOP = 8;
FrameWidth = 768;
FrameHeight = 1024;
NumOfGop = TotalFrameNum/GOP;
SNR = [-5 0 5 10 15 20 25];
psnr = zeros(NumOfGop, numel(SNR));
BandwidthRatio = 1;
BlockNum = 64;
%%
imginfo.H       = FrameHeight;
imginfo.W       = FrameWidth;
imginfo.cH      = imginfo.H / 2;
imginfo.cW      = imginfo.W / 2;
imginfo.blkH    = 16;
imginfo.blkW    = 16;
imginfo.Ysz     = imginfo.H * imginfo.W;
imginfo.Usz     = imginfo.cH * imginfo.cW;
imginfo.Vsz     = imginfo.cH * imginfo.cW;
%%
for indGOP = 1:NumOfGop
    Pics    = fx_LoadNFrm (fid, imginfo, GOP);
    C       = DCT3(Pics - 128);
    if BandwidthRatio == 1
        L = VarEqualChunk(C,BlockNum);
    else
        avgDropN = round((1-BandwidthRatio) * GOP * imginfo.W * imginfo.H);
        [L, DroppedN] = fx_DropLeastE64Coefficient(C, avgDropN);
        if DroppedN < avgDropN
            disp('dropping 3d-DCT coefficients error');
        end
    end
    
    G = EnergyAllocationMethod_3D_Y(1, L);

    S = G.*C;
    
    for ii = 1:numel(SNR)
        EsN0dB = SNR(ii);
        %% transmission
        [Y, sigma_n] = fx_TxAWGN(S, EsN0dB);
        
        %% denoise and power de-allocation
        Chat = fx_dec_LLSE(Y, G, L, [], sigma_n, 0);
        rxPics = IDCT3(Chat) + 128;
        rxPics = round(rxPics);
        rxPics(rxPics > 255) = 255;
        rxPics(rxPics < 0)   = 0;
        %% calculate psnr
        psnr(indGOP, ii) = fx_CalcPSNR(Pics, rxPics);
        disp([indGOP EsN0dB psnr(indGOP, ii)]);
    end
end
fclose('all');