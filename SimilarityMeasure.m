close all;
clear
datestr(now)
addpath('modelspecific');
addpath('multigs');
addpath('MatlabFns\Projective');

cd vlfeat-0.9.14/toolbox;
feval('vl_setup');
cd ../..;

cd multigs;
if exist('computeIntersection','file')~=3
    mex computeIntersection.c;
end
cd ..;
%%
filepath='C:\Users\wu\Desktop\Thesixth semester\mimo\3DVSource\3Dtest\3DV测试序列\';
ind=[1 3 5];
FrameHeight = 1024;
FrameWidth = 768;
TotalFrameNum = 32;
imginfo.GOP = 8;% GOP = 8;
NumOfGop = TotalFrameNum/imginfo.GOP;
SNR = -5:5:20;
psnr = zeros(NumOfGop, numel(SNR));
simCoefNum = zeros(NumOfGop,imginfo.GOP);
BandwidthRatio = 1;
imginfo.BlockNum = 64;
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
fid=zeros(2,1);
for indVideo=1:2
    fid(indVideo) = fopen([filepath 'balloons' num2str(2*indVideo-1) '.yuv'],'rb');
end
simCoef=zeros(imginfo.H,imginfo.W,imginfo.GOP*NumOfGop);
for indGOP = 1:NumOfGop
    for ii=1:2
        Pics    = fx_LoadNFrm (fid(ii), imginfo, imginfo.GOP);
        C       = DCT3(Pics - 128);
        if ii==1
            C_balloons1=C;
        end
    end
    for indFram=1:imginfo.GOP
      %% 对img2做warp/对每一帧做warp
%        img{1}=C_balloons1(:,:,indFram);
%        img{2}=C(:,:,indFram);
%        [ H1 ] = rmatrix( img );
%        [imgref_warp,imref_invwarp,x,y,xii,yii]= imTransD_wj(img{1}, H1, size(img{1}));
       %% 数据处理
       test1=C_balloons1(:,:,indFram);
       test2=C(:,:,indFram);
       LikelyhoodCoef=abs(test1-test2)./abs(test1);
       LikelyhoodCoef(LikelyhoodCoef>=0.1)=1;
       LikelyhoodCoef(LikelyhoodCoef<0.1)=0;
       fprintf('similar coef. GOP: %d, Frame:%d, Num : %d !\n',indGOP,indFram,imginfo.H*imginfo.W-sum(LikelyhoodCoef(:)));
       simCoefNum(indGOP,indFram)=imginfo.H*imginfo.W-sum(LikelyhoodCoef(:));
       simCoef(:,:,(indGOP-1)*imginfo.GOP+indFram)=LikelyhoodCoef;
    end
end





