function xx1 = fx_Recwj(Xr,bh,bw)
    rx_img{1}=Xr([1 2 3 4],:);
    rx_img{2}=Xr([5 6 7 8],:);
    rx_line_tmp=[];
    rx_line=[];
    seg=size(Xr,2)/64/8;
    for ii=1:1
        rx=rx_img{ii};
        for mm=1:64
            for jj=seg*(mm-1)+1:seg*mm
                for kk=1:4
                    rx_line_tmp=[rx_line_tmp rx(kk,jj)];
                end
            end
%             rx_invzig=invzigzag(rx_line_tmp,cfg.imginfo.cH,cfg.imginfo.cW);
            rx_line = [rx_line;rx_line_tmp];
            rx_line_tmp=[];
        end
     %% reshape 
            z1 = [];
            tt = [];
            for mm = 1:64
                currentBlock=invzigzag(rx_line(mm,:),bh,bw);
                tt = [tt currentBlock];
                if mod(mm,8) == 0 %µ÷Õû¿í¶È
                    z1 = [z1;tt];
                    tt = [];
                end
            end
%             warp_coset=meta.d(1:end/2,:);
%             z1=z1+warp_coset*meta.dsc_mod;
            
            xx=idct2(z1)+128;
            xx1=round(xx);
            xx1(xx1>255)=255;
            xx1(xx1<0)=0;
    end
    
end