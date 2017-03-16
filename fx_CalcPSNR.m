function psnr = fx_CalcPSNR(src, rec, varargin)

if ischar(src) && ischar(rec)
    imginfo  = varargin{1};
    FrameNum = varargin{2};
    psnr = PSNR_Video(src, rec, imginfo.H, imginfo.W, FrameNum);
else
    psnr = 0;
    for ii = 1:size(src,3)
        psnr = psnr + 10*log10(255^2/mean(mean((src(:,:,ii) - rec(:,:,ii)).^2)));
    end
    psnr = psnr / size(src,3);
end