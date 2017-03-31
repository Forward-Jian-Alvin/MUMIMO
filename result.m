%% 2017/03/30
close all;
clear;
datestr();
snr = -5:5:20;
softcast_ballons1=[32.1847 36.5923 41.1183 45.6090 50.0263 54.3137];
softcast_ballons3=[32.2950   36.6781   41.1750   45.6886   50.0968   54.3746];
plot(snr,softcast_ballons1,'-r+',snr,softcast_ballons3,'-b*');
%% 