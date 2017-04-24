%% 2017/03/30
close all;
clear;
% datestr(now);
%snr = -5:5:20;
%softcast_ballons1=[32.1847 36.5923 41.1183 45.6090 50.0263 54.3137];
%softcast_ballons3=[32.2950   36.6781   41.1750   45.6886   50.0968   54.3746];
%plot(snr,softcast_ballons1,'-r+',snr,softcast_ballons3,'-b*');
%% tests
Dc=0;
SimSource=0;
Normal=0;
Manned=1;
%% 
snr=-5:5:15;
Proposed0_Video1=[25.9205 29.9364 34.2092 38.5731 43.0958];
Proposed0_Video2=[24.2652 28.1791 32.4688 36.8025 41.3295];
Proposed005_Video1=[25.9622 29.9804 34.3026 38.6052 42.9694];
Proposed005_Video2=[24.3069 28.2796 32.5812 36.8771 41.1892];
Proposed010_Video1=[26.1013 30.0645 34.3018 38.5291 42.6892];
Proposed010_Video2=[24.4635 28.3693 32.5135 36.7606 40.8700];
Proposed015_Video1=[26.2638 30.1100 34.1224 37.9949 41.4218];
Proposed015_Video2=[24.6693 28.4099 32.3297 36.2162 39.4874];
Proposed020_Video1=[26.2171 30.0363 33.9180 37.5399 40.3589];
Proposed020_Video2=[24.5811 28.3142 32.1853 35.6850 38.3825];

Proposed_SameVideo1=[31.5969  35.9370 40.4441 44.9502 49.4214];
Proposed_SameVideo2=[31.5737  36.0538 40.5297 45.0191 49.4789];

Proposed015_Video1_Dc=[26.2533 30.0417 34.1809 38.0899 41.5566];
Proposed015_Video2_Dc=[24.6230 28.3396 32.3892 36.2336 39.5041];

ProposedManned_Video1=[30.4202 33.9246 37.5748 41.4809 45.5444];
ProposedManned_Video2=[29.8860 32.9971 36.2913 39.8033 43.2803];

if Manned
    figure();plot(snr,Proposed0_Video1,'-r+',snr,Proposed005_Video1,'-g+',...
        snr,Proposed010_Video1,'-b+',snr,Proposed015_Video1,'-m+',snr,Proposed020_Video1,'-k+',snr,ProposedManned_Video1,'-k+',snr,Proposed_SameVideo1,'-b*');
    legend('0','0,05','0.10','0.15','0.20','First Gop Same','uplimit');
    figure();plot(snr,Proposed0_Video2,'-r+',snr,Proposed005_Video2,'-g+',...
        snr,Proposed010_Video2,'-b+',snr,Proposed015_Video2,'-m+',snr,Proposed020_Video2,'-k+',snr,ProposedManned_Video2,'-k+',snr,Proposed_SameVideo2,'-b*');
    legend('0','0,05','0.10','0.15','0.20','First Gop Same','upLimit');
end
if Normal
    figure();plot(snr,Proposed0_Video1,'-r+',snr,Proposed005_Video1,'-g+',...
        snr,Proposed010_Video1,'-b+',snr,Proposed015_Video1,'-m+',snr,Proposed020_Video1,'-k+');
    legend('0','0,05','0.10','0.15','0.20');
    figure();plot(snr,Proposed0_Video2,'-r+',snr,Proposed005_Video2,'-g+',...
        snr,Proposed010_Video2,'-b+',snr,Proposed015_Video2,'-m+',snr,Proposed020_Video2,'-k+');
    legend('0','0,05','0.10','0.15','0.20');
end
if SimSource
    figure();plot(snr,Proposed0_Video1,'-r+',snr,Proposed005_Video1,'-g+',...
        snr,Proposed010_Video1,'-b+',snr,Proposed015_Video1,'-m+',snr,Proposed020_Video1,'-k+',snr,Proposed_SameVideo1,'-bs');
    legend('0','0,05','0.10','0.15','0.20','Same Source');
    figure();plot(snr,Proposed0_Video2,'-r+',snr,Proposed005_Video2,'-g+',...
        snr,Proposed010_Video2,'-b+',snr,Proposed015_Video2,'-m+',snr,Proposed020_Video2,'-k+',snr,Proposed_SameVideo2,'-bs');
    legend('0','0,05','0.10','0.15','0.20','Same Source');
end
if Dc
    figure();plot(snr,Proposed015_Video1,'-m+',snr,Proposed015_Video2,'-b+',snr,Proposed015_Video1_Dc,'-ms',snr,Proposed015_Video2_Dc,'-bs');
    legend('Video1','Video2','Video1-dc','Video2-dc');
end
