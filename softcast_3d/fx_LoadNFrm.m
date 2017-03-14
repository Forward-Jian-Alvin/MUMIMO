% function [y u v] = fx_LoadNFrm (fid, imginfo, N)
function y  = fx_LoadNFrm (fid, imginfo, N)
y = zeros(imginfo.W, imginfo.H, N);
u = zeros(imginfo.cW, imginfo.cH, N);
v = zeros(imginfo.cW, imginfo.cH, N);
for ii = 1:N
%     [r s t] = fx_LoadYUV1Frm (fid, imginfo);
    r  = fx_LoadYUV1Frm (fid, imginfo);
    y(:,:,ii) = reshape(r, imginfo.W, imginfo.H);
%     u(:,:,ii) = reshape(s, imginfo.cW, imginfo.cH);
%     v(:,:,ii) = reshape(t, imginfo.cW, imginfo.cH);
end