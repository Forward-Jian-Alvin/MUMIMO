function [ H1 ] = rmatrix( img1,img2 )
%RMATRIX Summary of this function goes here
%   Detailed explanation goes here
    clear global;
    global fitfn resfn degenfn psize numpar
    fitfn = 'homography_fit';
    resfn = 'homography_res';
    degenfn = 'homography_degen';
    psize   = 4;
    numpar  = 9;
%%
    M = 500;
    thr   = 0.1;  % RANSAC threshold.
    [kp1 ds1] = vl_sift(single(img1),'PeakThresh', 0,'edgethresh',50); %edit by yhj; original is 500
    [kp2 ds2] = vl_sift(single(img2),'PeakThresh', 0,'edgethresh',50); %edit by yhj; original is 500
%     figure(1);subplot(2,2,1);imshow(img{2},[]);title('img_ori image');
%     figure(1);subplot(2,2,2);imshow(img{1},[]);title('img_ref image');
    matches   = vl_ubcmatch(ds1,ds2);
    data_orig = [ kp1(1:2,matches(1,:)) ; ones(1,size(matches,2)) ; kp2(1:2,matches(2,:)) ; ones(1,size(matches,2)) ];

    [ dat_norm_img1 T1 ] = normalise2dpts(data_orig(1:3,:));
    [ dat_norm_img2 T2 ] = normalise2dpts(data_orig(4:6,:));
    data_norm = [ dat_norm_img1 ; dat_norm_img2 ];

    [ par res inx tim ] = multigsSampling(100,data_norm,M,10);
    con = sum(res<=thr);
    [ ~, maxinx ] = max(con);
    inliers = find(res(:,maxinx)<=thr);
    pt1_inlier = data_orig(1:2,inliers);
    pt2_inlier = data_orig(4:5,inliers);

    [h A D1 D2] = feval(fitfn,data_norm(:,inliers));
    H = T2\(reshape(h,3,3)*T1);
    H1 = pinv(H);
end

