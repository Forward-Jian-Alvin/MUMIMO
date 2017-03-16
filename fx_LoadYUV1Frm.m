function y = fx_LoadYUV1Frm (fid, imginfo)
imginfo.cH      = imginfo.H / 2;
imginfo.cW      = imginfo.W / 2;
imginfo.blkH    = 16;
imginfo.blkW    = 16;
imginfo.Ysz     = imginfo.H * imginfo.W;
imginfo.Usz     = imginfo.cH * imginfo.cW;
imginfo.Vsz     = imginfo.cH * imginfo.cW;
[yuv cnt] = fread (fid, imginfo.Ysz + imginfo.Usz + imginfo.Vsz, 'uint8');

if cnt < imginfo.Ysz + imginfo.Usz + imginfo.Vsz
    y = 0;
   
    return;
end

y = yuv(1:imginfo.Ysz);

