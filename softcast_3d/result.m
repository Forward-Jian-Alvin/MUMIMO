close all;
clear
snr = -5:5:25;
%% test sequence :'balloons1'
balloons1=[31.7985   36.1734   40.6844   45.2146   49.6599   53.9401   61.4082];
%% test sequence :'balloons3'
balloons3=[31.8680   36.2772   40.7393   45.2909   49.7190   53.9834   61.6224];
%% test sequence :'balloons5'
balloons5=[31.7876   36.1043   40.5898   45.1451   49.5746   53.8035   61.2268];
plot(snr,balloons1,'-r+',snr,balloons3,'-r*',snr,balloons5,'-rs');
plot(snr,balloons1,'-r+');
xlabel('SNR');ylabel('PSNR')

