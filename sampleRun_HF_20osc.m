fp=26000:2000:100000;
%fp=30000;
spf=20;
load('fcbw.mat')
bw=bwm_hf13(1:38);
nos_samples=8000;
fs=2e6;
folder=2;
zrec1=0.7;
zrec2=0.9;
% Receivers are 34,38,36,37
rec1no=36;
rec2no=37;
% gains are 6, 3.7, 5,
gain_corr=5;
% Channels are at 6,7,2,3
ch1=2;
ch2=3;
ws=5;
we=5000;
waveform=3;
block_fill=[1 0 0];
device=1;
