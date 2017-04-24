function psnr = Calpsnr(org, rec,NaNnum)
% USAGE:
%    psnr = calcpsnr(org, rec)
% DESCRIPTION:
%    Calculate PSNR between two signals whose dimension
%    can be 1, 2, or 3.
% INPUT:   
%    org  - original signal
%    rec  - reconstructed signal
% OUTPUT:
%    psnr - PSNR
% AUTHORS:
%    Jingyu Yang, Dept. of Automation, Tsinghua University, May 2005
%    yangjy03@mails.tsinghua.edu.cn
% SEE ALSO
%     calsnr
    tmp = org(:,:,1);
    if max(tmp(:))<=1.2
        org = uint8(org*255);
        rec = uint8(rec*255);
    end


    [h, w, k] = size(org);
    mse = sum(sum((double(org) - double(rec)).^2, 2), 1)/(h*w-NaNnum);
    %mse = sum(sum(sum((double(org) - double(rec)).^2,3), 2), 1)/(h*w*k);
    % mse = reshape(mse, 1,[]);
    psnr = mean(10*log10((255^2)./mse));
end

