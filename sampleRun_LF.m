fp=10000:2000:24000;
%fp=26000;
spf=20;
load('fcbw.mat')
bw=bwm_lf11(1:16);
bw(17)=bw(16);
bw(18)=bw(16);
nos_samples=8000;
fs=2e6;
folder=7;
zrec1=0.3;
zrec2=0.5;
% Receivers are 34,38,36,37
rec1no=34;
rec2no=38;
% gains are 25.3, -10.6, 5,
gain_corr=14.7;
% Channels are at 6,7,2,3
ch1=6;
ch2=7;
ws=1;
we=7000;
waveform=1;
block_fill=[1 0 0];
device=3;
