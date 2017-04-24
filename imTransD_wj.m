% IMTRANSD - Homogeneous transformation of an image.
%
% This is a stripped down version of imTrans which does not apply any origin
% shifting to the transformed image
%
% Applies a geometric transform to an image
%
%  newim = imTransD(im, T, sze, lhrh);
%
%  Arguments: 
%        im     - The image to be transformed.
%        T      - The 3x3 homogeneous transformation matrix.
%        sze    - 2 element vector specifying the size of the image that the
%                 transformed image is placed into.  If you are not sure
%                 where your image is going to 'go' make sze large! (though
%                 this does not help you if the image is placed at a negative
%                 location) 
%       lhrh    - String 'lh' or 'rh' indicating whether the transform was
%                 computed assuming columns represent x and rows represent y
%                 (a left handed coordinate system) or if it was computed
%                 using rows as x and columns as y (a right handed system,
%                 albeit rotated 90 degrees).  The default is assumed 'lh'
%                 though 'rh' is probably more sensible.
%
%
%  Returns:
%        newim  - The transformed image.
%
% See also: IMTRANS
%

% Copyright (c) 2000-2005 Peter Kovesi
% School of Computer Science & Software Engineering
% The University of Western Australia
% http://www.csse.uwa.edu.au/
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.

% April 2000 - Original version.
% April 2010 - Allowance for left hand and right hand coordinate systems.
%              Offset of 1 pixel that was (incorrectly) applied in
%              transformImage removed.

function [newim,invim,x,y,xii,yii] = imTransD_wj(im, T, sze, lhrh);

    if ~exist('lhrh','var'), lhrh = 'l'; end
    
    if isa(im,'uint8')
        im = double(im);  % Make sure image is double     
    end
    
    if lhrh(1) == 'r'     % Transpose the image allowing for colour images
        im = permute(im,[2 1 3]);
    end
    
    threeD = (ndims(im)==3);  % A colour image
    if threeD    % Transform red, green, blue components separately
        im = im/255;  
        r = transformImage_wj(im(:,:,1), T, sze);
        g = transformImage_wj(im(:,:,2), T, sze);
        b = transformImage_wj(im(:,:,3), T, sze);
        
        newim = repmat(uint8(0),[size(r),3]);
        newim(:,:,1) = uint8(round(r*255));
        newim(:,:,2) = uint8(round(g*255));
        newim(:,:,3) = uint8(round(b*255));
        
    else                % Assume the image is greyscale
        [newim,invim,x,y,xii,yii] = transformImage_wj(im, T, sze);
    end
    
    if lhrh(1) == 'r'   % Transpose back again
        newim = permute(newim,[2 1 3]);    
    end
    
%------------------------------------------------------------

% The internal function that does all the work

function [newim,invim,x,y,xii,yii] = transformImage_wj(im, T, sze);
    
    [rows, cols] = size(im);
    
    % Set things up for the image transformation.
    newim = zeros(sze(1),sze(2));
    [xi,yi] = meshgrid(1:sze(2),1:sze(1));    % All possible xy coords in the image.
        
    % Transform these xy coords to determine where to interpolate values
    % from. 
    Tinv = inv(T);
    sxy = homoTrans(Tinv, [xi(:)' ; yi(:)' ; ones(1,sze(1)*sze(2))]);% 画图坐标
    
    xi = reshape(sxy(1,:),sze(1),sze(2));
    yi = reshape(sxy(2,:),sze(1),sze(2));
    
    [x,y] = meshgrid(1:cols,1:rows);% 幕布x，y坐标
 
%    x = x-1; % Offset x and y relative to region origin.
%    y = y-1; 
    newim = interp2(x,y,im,xi,yi,'cubic'); % Interpolate values from source image.
% xi,yi坐标上的像素值已知，插值到x，y上。
        figure(1);subplot(2,2,3);imshow(newim,[]);title('warped image');
    %imwrite(uint8(newim),'.\images\result\warped.png');
    % Place new image into an image of the desired size
    %newim = implace(zeros(sze),newim,0,0);
    %% inverse interplotion
    [xii,yii] = meshgrid(1:sze(2),1:sze(1));
    xy_wj= homoTrans(T, [xii(:)' ; yii(:)' ; ones(1,sze(1)*sze(2))]);
    xii = reshape(xy_wj(1,:),sze(1),sze(2));
    yii = reshape(xy_wj(2,:),sze(1),sze(2));
    invim = interp2(x,y,newim,xii,yii,'cubic');
    
    figure(1);subplot(2,2,4);imshow(invim,[]);title('inverse warped image');
    %imwrite(uint8(invim),'.\images\result\inv_warped.png');
    %% calculate psnr
    NaNNum=0;
    for ii=1:rows
        for jj=1:cols
            if (~(invim(ii,jj)>=0&&invim(ii,jj)<=255))
                invim(ii,jj)=im(ii,jj);
                NaNNum=NaNNum+1;
            end
        end
    end
    invpsnr=InterpCalpsnr(invim,im,NaNNum);
    disp(invpsnr);