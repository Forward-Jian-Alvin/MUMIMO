close all;
clear
datestr(now)
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
%*******************************Encoder***********************************%
simCoef=zeros(64,imginfo.H*imginfo.W/64,imginfo.GOP*NumOfGop);

for indGOP = 1%:NumOfGop
    Tx=[];
    TxLambda=[];
    TxG=[];
    for ii=1:2
        Pics    = fx_LoadNFrm (fid(ii), imginfo, imginfo.GOP);
        C       = DCT3(Pics - 128);
        if ii==1
            C_balloons1=C;
        end
        % SOFTCAST Encoder
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
        TxLambda=[TxLambda;lambda];
        TxG=[TxG;g];
        Tx=[Tx;enc];
    end
    Threshold=0.1;
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
    indGvalidsum = zeros(size(indGvalid));
    indGInvalidsum = zeros(size(indGvalid));
    for mm=1:length(indGvalid)
        tmp=mm;
        while(tmp~=0)
            indGvalidsum(mm)=indGvalidsum(mm)+indGvalid(tmp);
            tmp=tmp-1;
        end
        indGInvalidsum(mm)=12288*mm-indGvalidsum(mm);
    end
    indGInvalidsum(end)=4-mod(length(find(LikelyhoodCoef==1)),4)+indGInvalidsum(end);
    
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
            for tt=1:size(Yr,2)
                % search g
                for ind4=1:4
                    indg=4*(tt-1)+ind4;
                    for ss=length(TxG)/2:-1:1
                        if(indg<=indGvalidsum(ss))
                           gValid1(ind4)=TxG(ss);
                           gValid2(ind4)=TxG(ss+512);
                           lambdaValid1(ind4)=TxLambda(ss);
                           lambdaValid2(ind4)=TxLambda(ss+512);
                        end
                    end
                end
                D1 = blkdiag(Hhat, Hhat) * diag([gValid1(1) gValid1(2) gValid1(3) gValid1(4)]);
                L1=[lambdaValid1(1) lambdaValid1(2) lambdaValid1(3) lambdaValid1(4)];
                DecMat1 = diag(L1) * D1' * pinv( D1 * diag(L1) * D1'  + Ne * eye(size(D1,1)));
                Xrv1(:,tt)= DecMat1 * Yr([1 3 2 4],tt);
                
                D2 = blkdiag(Hhat, Hhat) * diag([gValid2(1) gValid2(2) gValid2(3) gValid2(4)]);
                L2=[lambdaValid2(1) lambdaValid2(2) lambdaValid2(3) lambdaValid2(4)];
                DecMat2 = diag(L2) * D2' * pinv( D2 * diag(L2) * D2'  + Ne * eye(size(D2,1)));
                Xrv2(:,tt)= DecMat2 * Yr([1 3 2 4],tt);
            end
            XrValid1=Xrv1(:);XrValid2=Xrv2(:);
            % step3: 解Invalid
            Yrwj=[real(InValidDec);imag(InValidDec)];
            Hr = [real(H) -imag(H); imag(H) real(H)];
            XrInv = zeros(8,size(InValidDec,2));
            gInValid=zeros(1,8);
            lambdaInValid=zeros(1,8);
            for tt=1:size(InValidDec,2)
                % search g
                for ind4=1:4
                    indg=4*(tt-1)+ind4;
                    for ss=length(TxG)/2:-1:1
                        if(indg<=indGInvalidsum(ss))
                           gInValid(ind4)=TxG(ss);
                           gInValid(ind4+4)=TxG(ss);
                           lambdaInValid(ind4)=TxLambda(ss);
                           lambdaInValid(ind4+4)=TxLambda(ss);
                        else
                           break;
                        end
                    end
                end
                D1 = Hr * diag([gInValid(1) gInValid(7) gInValid(2) gInValid(8)]);
                D2 = Hr * diag([gInValid(3) gInValid(5) gInValid(4) gInValid(6)]);
                D_temp=blkdiag(D1,D2);
                L=[lambdaInValid(1) lambdaInValid(7) lambdaInValid(2) lambdaInValid(8)
                    lambdaInValid(3) lambdaInValid(5) lambdaInValid(4) lambdaInValid(6)];
                DecMat=diag(L) * D_temp' * pinv(D_temp * diag(L) * D_temp' + blkdiag(N0 * eye(size(D_temp,1))));
                XrInv(:,tt)= DecMat * Yrwj([1:2 5:6 3:4 7:8],tt);
            end
            order=[1 3 5 7 6 8 2 4];
            Xtmp=XrInv(order,:);
            Xtmp(3,:)=-Xtmp(3,:);
            Xtmp(6,:)=-Xtmp(6,:);
            XrVideo1=Xtmp(1:4,:);XrVideo2=Xtmp(5:end,:);
            XrInValid1=XrVideo1(:);XrInValid2=XrVideo2(:);
            % step4: 合并解的数据 1024*12288
            
        end
    end


%*******************************Reconstruct*******************************%
    
end






