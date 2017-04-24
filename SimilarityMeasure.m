close all;
clear
datestr(now)
% clear global;
global fitfn resfn degenfn psize numpar

addpath('modelspecific');
addpath('multigs');
addpath('MatlabFns\Projective');
addpath('./libchannel');

cd vlfeat-0.9.14/toolbox;
feval('vl_setup');
cd ../..;

cd multigs;
if exist('computeIntersection','file')~=3
    mex computeIntersection.c;
end
cd ..;


% options :0 invalid;1 valid;
DcComponent=0;
MannedSource=0;
ROI=1;
%%
filepath='C:\Users\wu\Desktop\Thesixth semester\mimo\3DVSource\3Dtest\3DV测试序列\';
ind=[1 3 5];
FrameHeight = 1024;
FrameWidth = 768;
FrameHeightROI = 992;
FrameWidthROI = 760;
TotalFrameNum = 32;
imginfo.GOP = 8;% GOP = 8;
NumOfGop = TotalFrameNum/imginfo.GOP;
SNR = -5:5:15;
psnr1 = zeros(NumOfGop, numel(SNR));
psnr2 = zeros(NumOfGop, numel(SNR));
simCoefNum = zeros(NumOfGop,imginfo.GOP);
BandwidthRatio = 1;
imginfo.BlockNum = 64;
%%
imginfo.H       = FrameHeight;
imginfo.W       = FrameWidth;
imginfo.cH      = FrameHeight / 2;
imginfo.cW      = FrameWidth / 2;
imginfo.Ysz     = imginfo.H * imginfo.W;
imginfo.Usz     = imginfo.cH * imginfo.cW;
imginfo.Vsz     = imginfo.cH * imginfo.cW;

imginfo.HROI       = FrameHeightROI;
imginfo.WROI       = FrameWidthROI;
imginfo.bh = FrameHeightROI/8;
imginfo.bw = FrameWidthROI/8;

%%    
fid=zeros(2,1);
for indVideo=1:2
    fid(indVideo) = fopen([filepath 'balloons' num2str(2*indVideo-1) '.yuv'],'rb');
% test 输入相同source
%     fid(indVideo) = fopen([filepath 'balloons3.yuv'],'rb');
end
%*******************************Source Alignment***************************%
Pics_srcAlign1    = fx_LoadNFrm (fid(1), imginfo, imginfo.GOP);
Pics_srcAlign2    = fx_LoadNFrm (fid(2), imginfo, imginfo.GOP);
imgref_warpSet=zeros(imginfo.H, imginfo.W,imginfo.GOP);
for FrmId=1:8
    Frm1=Pics_srcAlign1(:,:,FrmId);
    Frm2=Pics_srcAlign2(:,:,FrmId);
    figure(1);subplot(2,2,1);imshow(Frm1,[]);title('balloons1');
    figure(1);subplot(2,2,2);imshow(Frm2,[]);title('balloons3');
    if FrmId==1
        fitfn = 'homography_fit';
        resfn = 'homography_res';
        degenfn = 'homography_degen';
        psize   = 4;
        numpar  = 9;
    %% feature-based image alignment
        M = 500;
        thr   = 0.1;  % RANSAC threshold.
        [kp1,ds1] = vl_sift(single(Frm1),'PeakThresh', 0,'edgethresh',50); %edit by yhj; original is 500
        [kp2,ds2] = vl_sift(single(Frm2),'PeakThresh', 0,'edgethresh',50); %edit by yhj; original is 500
    %     figure(1);subplot(2,2,1);imshow(img{2},[]);title('img_ori image');
    %     figure(1);subplot(2,2,2);imshow(img{1},[]);title('img_ref image');
        matches   = vl_ubcmatch(ds1,ds2);
        data_orig = [ kp1(1:2,matches(1,:)) ; ones(1,size(matches,2)) ; kp2(1:2,matches(2,:)) ; ones(1,size(matches,2)) ];

        [ dat_norm_img1,T1 ] = normalise2dpts(data_orig(1:3,:));
        [ dat_norm_img2,T2 ] = normalise2dpts(data_orig(4:6,:));
        data_norm = [ dat_norm_img1 ; dat_norm_img2 ];

        [ par,res,inx,tim ] = multigsSampling(100,data_norm,M,10);
        con = sum(res<=thr);
        [ ~, maxinx ] = max(con);
        inliers = find(res(:,maxinx)<=thr);
        pt1_inlier = data_orig(1:2,inliers);
        pt2_inlier = data_orig(4:5,inliers);

        [h,A,D1,D2] = feval(fitfn,data_norm(:,inliers));
        H = T2\(reshape(h,3,3)*T1);
        H1 = pinv(H);
        disp(H1);
    end
    [imgref_warp,imref_invwarp,x,y,xii,yii]= imTransD_wj(Frm2, H1, size(Frm2));
    % 取整可能带来的误差
    imgref_warp=round(imgref_warp);
    imgref_warp(imgref_warp>255)=255;
    imgref_warp(imgref_warp<0)=0;
    imgref_warpSet(:,:,FrmId)=imgref_warp;
    fprintf('Source Aligned!');
end

input('请输入一个数字，并按回车继续：');
%*******************************Encoder***********************************%
simCoef=zeros(64,imginfo.HROI*imginfo.WROI/64,imginfo.GOP*NumOfGop);

for indGOP = 1%:NumOfGop
    Tx=[];
    TxLambda=[];
    TxG=[];
%     Test_Rec=zeros(1024,12288);
    for ii=1:2
        Pics = fx_LoadNFrm (fid(ii), imginfo, imginfo.GOP);
        if ROI && ii==1
            Pics = Pics(33:end,1:760,:);
        elseif ROI && ii==2 
            Pics = imgref_warpSet(33:end,1:760,:);
        end
        C       = DCT3(Pics - 128);
        if ii==1
            Pics_balloons1=Pics;
            C_balloons1=C;
        end
        if MannedSource && ii==2
            C(1:imginfo.bh,1:imginfo.bw,1)=C_balloons1(1:imginfo.bh,1:imginfo.bw,1);
            Pics=IDCT3(C)+128;
        end
        
        % SOFTCAST Encoder
        x=[];
        for kk = 1:imginfo.GOP
            for tt1 = 1:8
                for jj=1:8
                    currentBlock = C((tt1-1)*imginfo.bh+1:tt1*imginfo.bh,(jj-1)*imginfo.bw+1:jj*imginfo.bw,kk);
                    x = [x;reshape(currentBlock,1,imginfo.bh*imginfo.bw)];
                end
            end
        end
        if MannedSource && ii==2
            x(1,:)=reshape(C_balloons1(1:imginfo.bh,1:imginfo.bw,1),1,imginfo.bh*imginfo.bw);
        end
        
        %******************test for Rec
%         Test_Rec(end/2*(ii-1)+1:end/2*ii,:)= x;
        %********************
        P1 = 1;% mean(mean(x.*x));
        P = P1*8*8*8;%total power constraint
        lambda = mean((x.*x)');
        lambda = lambda';
        g = sqrt(P/sum(sqrt(lambda)))./sqrt(sqrt(lambda));
        G = diag(g);
        enc=G*x;
        TxLambda=[TxLambda;lambda];
        TxG=[TxG;g];
        Tx=[Tx;enc];
    end
    Threshold=0.15;
    TxData1=Tx(1:end/2,:);
    TxData2=Tx(end/2+1:end,:);
    LikelyhoodCoef=abs(TxData1-TxData2)./abs(TxData1);
    LikelyhoodCoef(LikelyhoodCoef>=Threshold)=1;
    LikelyhoodCoef(LikelyhoodCoef<Threshold)=0;% valid
    LikelyhoodCoef=1-LikelyhoodCoef;
    fprintf('similar coef. Num : %d !\n',sum(LikelyhoodCoef(:)));

%     for indFram=1:imginfo.GOP
%       %% 对img2做warp/对每一帧做warp
% %        img{1}=C_balloons1(:,:,indFram);
% %        img{2}=C(:,:,indFram);
% %        [ H1 ] = rmatrix( img );
% %        [imgref_warp,imref_invwarp,x,y,xii,yii]= imTransD_wj(img{1}, H1, size(img{1}));
%        %% 数据处理
%        test1=Tx(64*(indFram-1)+1:64*indFram,:);
%        test2=Tx(64*(indFram-1)+1+8*64:64*indFram+8*64,:);
%        LikelyhoodCoefFrm=abs(test1-test2)./abs(test1);
%        LikelyhoodCoefFrm(LikelyhoodCoefFrm>=0.1)=1;
%        LikelyhoodCoefFrm(LikelyhoodCoefFrm<0.1)=0;
%        fprintf('similar coef. GOP: %d, Frame:%d, Num : %d !\n',indGOP,indFram,imginfo.H*imginfo.W-sum(LikelyhoodCoefFrm(:)));
%        simCoefNum(indGOP,indFram)=imginfo.H*imginfo.W-sum(LikelyhoodCoefFrm(:));
%        simCoef(:,:,(indGOP-1)*imginfo.GOP+indFram)=LikelyhoodCoefFrm;
%     end
%需要发送的数据包括：1 元数据：TxLambda& TxG 1024;2 相似矩阵：LikelyhoodCoef 512*12288;3 模拟数据：Tx:2*64*8*12288
% 组帧 1 valid else invalid
% 补足成4的倍数：复数对，stbc组码
    InvalidSet1=zeros(length(find(LikelyhoodCoef==0))+4-mod(length(find(LikelyhoodCoef==0)),4),1);
    validSet1=zeros(length(find(LikelyhoodCoef==1))+4-mod(length(find(LikelyhoodCoef==1)),4),1);
    InvalidSet2=zeros(length(find(LikelyhoodCoef==0))+4-mod(length(find(LikelyhoodCoef==0)),4),1);
    validSet2=zeros(length(find(LikelyhoodCoef==1))+4-mod(length(find(LikelyhoodCoef==1)),4),1);
    indGvalid=sum(LikelyhoodCoef,2);indGvalid(end)=4-mod(length(find(LikelyhoodCoef==1)),4)+indGvalid(end);
%     indGvalidsum = zeros(size(indGvalid));
%     indGInvalidsum = zeros(size(indGvalid));
%     for mm=1:length(indGvalid)
%         tmp=mm;
%         while(tmp~=0)
%             indGvalidsum(mm)=indGvalidsum(mm)+indGvalid(tmp);
%             tmp=tmp-1;
%         end
% %         indGInvalidsum(mm)=12288*mm-indGvalidsum(mm);
%     end
    indGInvalid=sum(1-LikelyhoodCoef,2);indGInvalid(end)=indGInvalid(end)+4-mod(length(find(LikelyhoodCoef==0)),4);
%     indGInvalidsum(end)=4-mod(length(find(LikelyhoodCoef==1)),4)+indGInvalidsum(end)+4-mod(length(find(LikelyhoodCoef==0)),4);
    
    %
    indexValid=1;
    indexinValid=1;
    for sfx=1:size(Tx,1)/2
        for sfy=1:size(Tx,2)
            if(LikelyhoodCoef(sfx,sfy))
                validSet1(indexValid)=Tx(sfx,sfy);
                validSet2(indexValid)=Tx(sfx+end/2,sfy);
                indexValid=indexValid+1;
            else
                InvalidSet1(indexinValid)=Tx(sfx,sfy);
                InvalidSet2(indexinValid)=Tx(sfx+end/2,sfy);
                indexinValid=indexinValid+1;
            end
        end
    end
    % test
    %Es_test1=(sum(abs(validSet1).^2)+sum(abs(validSet2).^2)+sum(abs(InvalidSet1).^2)+sum(abs(InvalidSet2).^2))/numel(Tx);
    % 组成复数
    validSetComp1 = validSet1(1:2:end)+1i* validSet1(2:2:end);
    InvalidSetComp1 = InvalidSet1(1:2:end)+1i* InvalidSet1(2:2:end);
    validSetComp2 = validSet2(1:2:end)+1i* validSet2(2:2:end);
    InvalidSetComp2 = InvalidSet2(1:2:end)+1i* InvalidSet2(2:2:end);
    % test
    %Es_test2=(sum(abs(validSetComp1).^2)+sum(abs(InvalidSetComp1).^2)+sum(abs(validSetComp2).^2)+sum(abs(InvalidSetComp2).^2))/numel(Tx);
    % validSetComp组成stbcodes
    % step1 ：对于前512个数据，x2取共轭取反,对于后512个数据，x1取共轭，x1,x2互换位置
    line1=validSetComp1(1:2:end).';
    line2=-conj(validSetComp1(2:2:end)).';
    line3=conj(validSetComp2(1:2:end)).';
    line4=validSetComp2(2:2:end).';
    % step2 ：
    line1_invalid=InvalidSetComp1(1:2:end).';
    line2_invalid=InvalidSetComp1(2:2:end).';
    line3_invalid=InvalidSetComp2(1:2:end).';
    line4_invalid=InvalidSetComp2(2:2:end).';
    % test
    %Es_test3=(sum(abs(line1).^2)+sum(abs(line2).^2)+sum(abs(line3).^2)+sum(abs(line4).^2)+...
    %    sum(abs(line1_invalid).^2)+sum(abs(line2_invalid).^2)+sum(abs(line3_invalid).^2)+sum(abs(line4_invalid).^2))/numel(Tx);
   TxData=[[line1 line1_invalid];[line4 line4_invalid];[line2 line2_invalid];[line3 line3_invalid]];
%*******************************Channel***********************************%
fprintf('**Channel**\n');
%% Init channel fading
    Nt =2;Nr=2;
    Es = mean(mean(abs(TxData(:).^2)));
    for indCoherent=1
        for ii = 1:numel(SNR)
            EsN0dB = SNR(ii);
           %% MIMO simulation
            fx_MIMOchannelReset(indCoherent);
            % Rayleigh fading
            maxNc1GOP=1024;
            fx_MIMOchannelInit(1, maxNc1GOP);
            H = fx_MIMOchannelLoad(Nr);
            H = H(:,1:Nt);
            Yc = blkdiag(H, H) * TxData + fx_AWGN(Nt*Es, EsN0dB, [Nr/Nt*size(TxData,1) size(TxData,2)]);
            N0 = Nt * Es * 10^(-EsN0dB/10) / Nr;
%             Yc = blkdiag(H, H) * TxData ;
%             N0 = 0;
%*******************************Decoder***********************************%
fprintf('**Channel SNR:%d,Channel State:%d**\n',EsN0dB,indCoherent); 
            % step1: 区分出valid 和 Invalid data
            ValidDec=Yc(:,1:length(line1));
            InValidDec=Yc(:,length(line1)+1:end);
            % step2: 解valid:得到stbc译码结果的Yr后，LLSE译码
            T1 = [0 1; -1 0];
            Yhat=0;
            Hhat=0;
            for jj = 1:2
                Hdiv = [H(jj,:); conj(H(jj,:))*T1'];
                Yhat = Yhat + Hdiv' * [ValidDec(jj,:); conj(ValidDec(jj+2,:))];
                Hhat = Hhat + Hdiv' * Hdiv;
            end
            Pchan = H(:)'*H(:);
            Hhat = real(Hhat) / Pchan;
            Yhat = Yhat / Pchan;
            Ne = N0 / Pchan;
            Yr = [real(Yhat); imag(Yhat)];
            Xrv1 = zeros(4,size(Yr,2)); Xrv2 = zeros(4,size(Yr,2));
            gValid1=zeros(1,4);gValid2=zeros(1,4);
            lambdaValid1=zeros(1,4);lambdaValid2=zeros(1,4);
            
            ValidGset1=[];ValidGset2=[];lambdaValidset1=[];lambdaValidset2=[];
            for ss=1:512
                ValidGset1=[ValidGset1 TxG(ss)*ones(1,indGvalid(ss))];
                ValidGset2=[ValidGset2 TxG(ss+512)*ones(1,indGvalid(ss))];
                lambdaValidset1=[lambdaValidset1 TxLambda(ss)*ones(1,indGvalid(ss))];
                lambdaValidset2=[lambdaValidset2 TxLambda(ss+512)*ones(1,indGvalid(ss))];
            end
            for tt1=1:size(Yr,2)
                % search g
                indg=4*(tt1-1)+(1:4);
                gValid1(:)=ValidGset1(indg);
                gValid2(:)=ValidGset2(indg);
                lambdaValid1(:)=lambdaValidset1(indg);
                lambdaValid2(:)=lambdaValidset2(indg);
%                 for ind4=1:4
%                     indg=4*(tt-1)+ind4;
%                     for ss=length(TxG)/2:-1:1
%                         if(indg<=indGvalidsum(ss))
%                            gValid1(ind4)=TxG(ss);
%                            gValid2(ind4)=TxG(ss+512);
%                            lambdaValid1(ind4)=TxLambda(ss);
%                            lambdaValid2(ind4)=TxLambda(ss+512);
%                         end
%                     end
%                 end
                D1 = blkdiag(Hhat, Hhat) * diag([gValid1(1) gValid1(2) gValid1(3) gValid1(4)]);
                L1=[lambdaValid1(1) lambdaValid1(2) lambdaValid1(3) lambdaValid1(4)];
                DecMat1 = diag(L1) * D1' * pinv( D1 * diag(L1) * D1'  + Ne * eye(size(D1,1)));
                Xrv1(:,tt1)= DecMat1 * Yr([1 3 2 4],tt1);
                
                D2 = blkdiag(Hhat, Hhat) * diag([gValid2(1) gValid2(2) gValid2(3) gValid2(4)]);
                L2=[lambdaValid2(1) lambdaValid2(2) lambdaValid2(3) lambdaValid2(4)];
                DecMat2 = diag(L2) * D2' * pinv( D2 * diag(L2) * D2'  + Ne * eye(size(D2,1)));
                Xrv2(:,tt1)= DecMat2 * Yr([1 3 2 4],tt1);
            end
            XrValid1=Xrv1(:);XrValid2=Xrv2(:);
            fprintf('**Valid Data Decoding**\n');
            % step3: 解Invalid
            % 1 把数据映射回去，2 矩阵求解 
            
            Yrwj=[real(InValidDec);imag(InValidDec)];
            Hr = [real(H) -imag(H); imag(H) real(H)];
            XrInv = zeros(8,size(InValidDec,2));
            gInValid=zeros(1,8);
            lambdaInValid=zeros(1,8);
            InvalidGset1=[];InvalidGset2=[];InvalidLambdaset1=[];InvalidLambdaset2=[];
            for ss=1:512
                InvalidGset1=[InvalidGset1 TxG(ss)*ones(1,indGInvalid(ss))];
                InvalidGset2=[InvalidGset2 TxG(ss+512)*ones(1,indGInvalid(ss))];
                InvalidLambdaset1=[InvalidLambdaset1 TxLambda(ss)*ones(1,indGInvalid(ss))];
                InvalidLambdaset2=[InvalidLambdaset2 TxLambda(ss+512)*ones(1,indGInvalid(ss))];
            end
            for tt1=1:size(InValidDec,2)
                indg=4*(tt1-1)+(1:4);
                gInValid(1:4)=InvalidGset1(indg);
                gInValid(5:8)=InvalidGset2(indg);
                lambdaInValid(1:4)=InvalidLambdaset1(indg);
                lambdaInValid(5:8)=InvalidLambdaset2(indg);
                % search g
%                 for ind4=1:4
%                     gInValid(ind4)=TxG(ss);
%                     for ss=length(TxG)/2:-1:1
%                         if(indg<=indGInvalidsum(ss))
%                            gInValid(ind4)=TxG(ss);
%                            gInValid(ind4+4)=TxG(ss);
%                            lambdaInValid(ind4)=TxLambda(ss);
%                            lambdaInValid(ind4+4)=TxLambda(ss);
%                         else
%                             
%                            break;
%                         end
%                     end
%                 end
%                 D_temp= diag([gInValid(1) gInValid(7) gInValid(2) gInValid(8) gInValid(3) gInValid(5) gInValid(4) gInValid(6)]);
                D1 = Hr * diag([gInValid(1) gInValid(7) gInValid(2) gInValid(8)]);
                D2 = Hr * diag([gInValid(3) gInValid(5) gInValid(4) gInValid(6)]);
                D_temp=blkdiag(D1,D2);
                L=[lambdaInValid(1) lambdaInValid(7) lambdaInValid(2) lambdaInValid(8) lambdaInValid(3) lambdaInValid(5) lambdaInValid(4) lambdaInValid(6)];
                DecMat=diag(L) * D_temp' * pinv(D_temp * diag(L) * D_temp' + blkdiag(N0 * eye(size(D_temp,1))));
                XrInv(:,tt1)= DecMat * Yrwj([1:2 5:6 3:4 7:8],tt1);
            end
            order=[1 3 5 7 6 8 2 4];
            Xtmp=XrInv(order,:);
%             Xtmp(3,:)=-Xtmp(3,:);
%             Xtmp(6,:)=-Xtmp(6,:);
            XrVideo1=Xtmp(1:4,:);XrVideo2=Xtmp(5:end,:);
            XrinValid1=XrVideo1(:);XrinValid2=XrVideo2(:);
            fprintf('**InValid Data Decoding**\n');
            % step4: 合并解的数据 1024*12288
            Rec=zeros(size(Tx,1),size(Tx,2));
            tmpv=1;tmpinv=1;
            for iirec=1:size(Tx,1)/2
                for jjrec=1:size(Tx,2)
                    if LikelyhoodCoef(iirec,jjrec)==1
                        Rec(iirec,jjrec)=XrValid1(tmpv);
                        Rec(iirec+512,jjrec)=XrValid2(tmpv);
                        tmpv=tmpv+1;
                    else
                        Rec(iirec,jjrec)=XrinValid1(tmpinv);
                        Rec(iirec+512,jjrec)=XrinValid2(tmpinv);
                        tmpinv=tmpinv+1;
                    end
                end
            end
%*******************************Reconstruct*******************************%
            % reshape 
            RecVideo1 = [];RecVideo2 = [];
            tt1 = [];tt2=[];
            for recii=1:imginfo.GOP
                for recjj = 1:imginfo.BlockNum
                    currentBlock1=reshape(Rec(recjj+(recii-1)*64,:),imginfo.bh,imginfo.bw);
                    currentBlock2=reshape(Rec(recjj+(recii-1)*64+512,:),imginfo.bh,imginfo.bw);
                    tt1 = [tt1 currentBlock1];
                    tt2 = [tt2 currentBlock2];
                    if mod(recjj,8) == 0 %调整宽度
                        RecVideo1 = [RecVideo1;tt1];
                        RecVideo2 = [RecVideo2;tt2];
                        tt1 = [];
                        tt2 = [];
                    end
                end
            end
            
            Rectmp1=zeros(FrameHeightROI,FrameWidthROI,imginfo.GOP);
            Rectmp2=zeros(FrameHeightROI,FrameWidthROI,imginfo.GOP);
            for Recmm=1:imginfo.GOP
                Rectmp1(:,:,Recmm)=RecVideo1(1+(Recmm-1)*end/imginfo.GOP:Recmm*end/imginfo.GOP,:);
                Rectmp2(:,:,Recmm)=RecVideo2(1+(Recmm-1)*end/imginfo.GOP:Recmm*end/imginfo.GOP,:);
            end
            if DcComponent
                Rectmp1(1,1,1)=C_balloons1(1,1,1);
                Rectmp2(1,1,1)=C(1,1,1);
            end
            rxLuma1 = IDCT3(Rectmp1) + 128;
            rxLuma2 = IDCT3(Rectmp2) + 128;
            rxLuma1 = round(rxLuma1);
            rxLuma1(rxLuma1 > 255) = 255;
            rxLuma1(rxLuma1 < 0)   = 0;
            rxLuma2 = round(rxLuma2);
            rxLuma2(rxLuma2 > 255) = 255;
            rxLuma2(rxLuma2 < 0)   = 0;
            %% calculate psnr   
            psnr1(indGOP, ii) = fx_CalcPSNR(Pics_balloons1, rxLuma1);
            disp([indGOP EsN0dB psnr1(indGOP, ii)]);
            psnr2(indGOP, ii) = fx_CalcPSNR(Pics, rxLuma2);
            disp([indGOP EsN0dB psnr2(indGOP, ii)]);
        end
    end
end
datestr(now)





