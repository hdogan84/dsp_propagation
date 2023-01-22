fp=10000:2000:24000;
%fp=10000;
spf=20;
load('fcbw.mat')
bw=bwm_lf13(1:16);
bw(17)=bw(16);
bw(18)=bw(16);
nos_samples=8000;
fs=2e6;
folder=8;
zrec1=0.3;
zrec2=0.7;
% Receivers are 34,38,36,37
rec1no=34;
rec2no=36;
% gains are 14.7,10.6,35
gain_corr=25.3;
% Channels are at 6,7,2,3
ch1=6;
ch2=2;
ws=1;
we=6000;
waveform=3;
block_fill=[1 0 0];
device=3;
