close all;
clear
datestr(now)
addpath('./libchannel');
filepath='C:\Users\wu\Desktop\Thesixth semester\mimo\3DVSource\3Dtest\3DV测试序列\';
ind=[1 3 5];
FrameHeight = 1024;
FrameWidth = 768;
TotalFrameNum = 8;
GOP = 8;
NumOfGop = TotalFrameNum/GOP;
SNR = [-5 0 5 10 15 20 25];
psnr = zeros(NumOfGop, numel(SNR));
BandwidthRatio = 1;
BlockNum = 64;
%%
imginfo.H       = FrameHeight;
imginfo.W       = FrameWidth;
imginfo.bh = FrameHeight/8;
imginfo.bw = FrameWidth/8;
imginfo.cH      = imginfo.H / 2;
imginfo.cW      = imginfo.W / 2;
imginfo.Ysz     = imginfo.H * imginfo.W;
imginfo.Usz     = imginfo.cH * imginfo.cW;
imginfo.Vsz     = imginfo.cH * imginfo.cW;
%%
Tx_line = [];
for ii=1:2
    fid = fopen([filepath 'balloons' num2str(ind(ii)) '.yuv'],'rb');
    for indGOP = 1:NumOfGop
        Pics    = fx_LoadNFrm (fid, imginfo, GOP);
        C       = DCT3(Pics - 128);
        [enc,lambda,g]  =fx_enc_wj(C,imginfo);
        meta.lambda=[meta.lambda;lambda];
        meta.g=[meta.g;g];
     %% zigag  
        currentBlock_complex=[];
        for jj=1:size(enc,1)
            warp_temp=reshape(enc(jj,:),imginfo.bh,imginfo.bw);
            warp_temp_zig=zigzag(warp_temp);
            warp_temp_comp=warp_temp_zig(1:2:end)+warp_temp_zig(2:2:end)*1i;
            currentBlock_complex=[currentBlock_complex warp_temp_comp];
        end
        %% two-way alamounti stbc
        Tx_line=[Tx_line;currentBlock_complex(1:2:end);currentBlock_complex(2:2:end)];     
    end
end
% stbc
T0 = eye(2);
T1 = [-1 0; 0 1];
X = [ T0*    [Tx_line(1,:);Tx_line(4,:)]; ...
                       T1*  [conj(Tx_line(2,:));conj(Tx_line(3,:))] ];% 保证平均信号功率为1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Init channel fading
Nt =2;Nr=2;
Es = X(1:end) * X(1:end)' / numel(X);
%% channel
for indCoherent=1
	for ii = 1:numel(SNR)
        EsN0dB = SNR(ii);
       %% MIMO simulation
        fx_MIMOchannelReset(indCoherent);
        % Rayleigh fading
        H = fx_MIMOchannelLoad(Nr);
        H = H(:,1:Nt);
        Yc = blkdiag(H, H) * X + fx_AWGN(Nt*Es, EsN0dB, [Nr/Nt*size(X,1) size(X,2)]);
        N0 = Nt * Es * 10^(-EsN0dB/10) / 2; % 2 coefficients are mapped to one complex symbol
        meta.N0=N0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %% stbc decoder
        Xr = fx_Decwj(Yc, H, meta);
        rxLuma = fx_Recwj(Xr,bh,bw);
       %% calculate psnr   
        MSE_imgref_warp = reshape(mean(mean((rxLuma - img{3}).^2, 1), 2), 1, 1);
        psnr_imgref_warp = [psnr_imgref_warp 10*log10(255^2 ./ MSE_imgref_warp)]; 
	end
end
imshow(rxLuma,[]);
disp(psnr_imgref_warp);
       
       
       

