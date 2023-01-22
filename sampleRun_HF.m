fp=26000:2000:100000;
%fp=30000;
spf=20;
load('fcbw.mat')
bw=bwm_hf11(1:38);
nos_samples=8000;
fs=2e6;
folder=1;
zrec1=0.3;
zrec2=0.5;
% Receivers are 34,38,36,37
rec1no=38;
rec2no=36;
% gains are 14.7, 0, 5,
gain_corr=0;
% Channels are at 6,7,2,3
ch1=7;
ch2=2;
ws=10;
we=5000;
waveform=1;
block_fill=[1 0 0];
device=1;
