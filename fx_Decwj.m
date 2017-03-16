function [ Xr ] = fx_Decwj( Yc, H, meta )
    GOP=8;
    blockNum=64;
    Hr = [real(H) -imag(H); imag(H) real(H)];
    Yrwj=[real(Yc);imag(Yc)];
    g=meta.g;
    lambda=meta.lambda;
    N0=meta.N0;
    len=size(Yrwj,2)/blockNum/GOP;% gop & blockNum
    Xr_Dec=[];
    for jj=1:8
        for ii=1:64
            gNum=GOP*blockNum;
            gtmp=ii+(jj-1)*blockNum;
            D1 = Hr * diag([g(gtmp) g(gtmp+gNum) g(gtmp) g(gtmp+gNum)]);
            D2 = Hr * diag([g(gtmp) g(gtmp+gNum) g(gtmp) g(gtmp+gNum)]);
            D_temp=blkdiag(D1,D2);
            Ltmp=[lambda(gtmp) lambda(gtmp+gNum)];
            L=[Ltmp Ltmp Ltmp Ltmp];
            DecMat=diag(L) * D_temp' * pinv(D_temp * diag(L) * D_temp' + blkdiag(N0 * eye(size(D_temp,1))));
            Xr_Dec= [Xr_Dec DecMat * Yrwj([1:2 5:6 3:4 7:8],len*(ii-1)+1:len*ii)];
        end
    end
    order=[1 3 5 7 6 8 2 4];
    Xr=Xr_Dec(order,:);
    Xr(3,:)=-Xr(3,:);
    Xr(6,:)=-Xr(6,:);
end

