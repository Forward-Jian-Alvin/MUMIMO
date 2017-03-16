% function [y u v] = fx_LoadNFrm (fid, imginfo, N)
function y  = fx_LoadNFrm (fid, imginfo, N)
y = zeros(imginfo.H, imginfo.W, N);
u = zeros(imginfo.cH, imginfo.cW, N);
v = zeros(imginfo.cH, imginfo.cW, N);
for ii = 1:N
%     [r s t] = fx_LoadYUV1Frm (fid, imginfo);
    r  = fx_LoadYUV1Frm (fid, imginfo);
    y(:,:,ii) = reshape(r, imginfo.H, imginfo.W);
%     u(:,:,ii) = reshape(s, imginfo.cW, imginfo.cH);
%     v(:,:,ii) = reshape(t, imginfo.cW, imginfo.cH);
end