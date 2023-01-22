%function [V,ALPHA,nfo]=prop_nov_time_stacked(fp,bw,spf,nos_samples,fs,folder,zrec1,zrec2,rec1no,rec2no,gain_corr,ch1,ch2,ws,we,waveform,block_fill,device);          %
tic

% main function used to process propagation data from the water tank, using
% signals received at two hydrophones
% This uses a time based correlation to detect arrival times (and so group velocity)
% and the amplitudes of the steady state portion of the signals on the two
% hydrophones to commpute attenuation. 
% exactly the same as prop_nov_time excepts uses prop_proc_time_stacked
% which stacks spf shots at each frequency prior to computing arival times
% and amplitudes. This improves SNR but lossess any inflrmation on the shot
% to shot variability (and so the variability errors)
% Need to adjust path to data in prop_proc_time
% uses sub routines get_rec_sens and spreading_noplot

% INPUTS: 
% fp=pump frequencies (Hz)
% bw: bandwidths for tonal pulses in pur water with one value for each
% spf=shots per frequency
% nos_samples=samples acquired per signal 
% fs=sampling frequency (Hz)
% folder: folder number where data stored
% zrec1: S-R separation for first (closer) hydrophone (m)
% zrec2: S-R separation for second (further) hydrophone (m)
% rec1no: last two digits of serial number of first hydrophone
% rec2no: last two digits of serila number of second hydrophone
% gain_corr: Difference in overall (wet and dry end) electronic gain between the first and 
% second hydrophone (dB)
% ch1: row in file in which contains data for first hydrophone
% ch2: row in file in which contains data for second hydrophone
% ws and we:rough start and end of window applied (data points): ensure
% this includes all data points that contain propagated signal on the
% closer and further hydrophone under analysis
% waveform: envelope of tonal wave 1=square, 2 =hanning
% pl: pulse length (seconds)
% block_fill: color for ploting, e.g. [1,0,0] is red
% device : pump source: use 1 for HF device, 2 for MF device and 3 for LF
% device

% OUTPUTS: 
% V group velocity (m/s) including mean value for all shots (row 1)
% and intrinsic error (row 2)
% ALPHA: attenuation (dB/m) including mean value for all shots (row 1)
%  and intrinsic error (row 2)
% nfo: numbersof shots that were successfullt processed at each frequency

% matricies required
nfo=zeros(length(fp),1);
V=zeros(2,length(fp));
ALPHA=zeros(2,length(fp));

% select pump device 
if device==1;radius=0.113/2;   % HF pump
elseif  device==2;radius=0.061/2;   % MF pump
elseif  device==3;radius=0.095/2;   % LF pump
else xxx
end
% pulse length
if waveform==1;PL=1e-3*ones(size(fp));
elseif waveform==2;PL=0.5e-3*ones(size(fp));
elseif waveform==3;PL=(20./fp).*ones(size(fp));
elseif waveform==4;PL=(30./fp).*ones(size(fp));
else xxx
end

% receiver sensitivities of hydrophones
[rs1]=get_rec_sens(fp,rec1no);
[rs2]=get_rec_sens(fp,rec2no);

% stuff
fny=fs/2;
te=1/fs;
t=[0:1/fs:(nos_samples-1)/fs];

% generate blank files
    toa=zeros(length(fp),1);
    tob=zeros(length(fp),1);
    AMP1=ones(length(fp),1);
    AMP2=ones(length(fp),1);
    ae=zeros(length(fp),spf);
    toc
    
for n=1:length(fp);fc=fp(n);BW=bw(n);fc
    % design filters
    %fac=2;[N,Wn]=buttord([fc-fac*BW fc+fac*BW]/fny,[fc-2*fac*BW fc+2*fac*BW]/fny,0.5,5);
    %[b,a]=butter(N,Wn);
    
    % for each frequenccy use sub routine prop_proc_time_stacked to get amplitude and
    % arrival time of signal on each hydrophone
    %[ta,tb,q,amp1,amp2]=prop_proc_time_stack(folder,fc,spf,fs,nos_samples,b,a,ch1,ch2,ws,we,PL(n),waveform,device);
    %   if ta>0; toa(n)=ta;tob(n)=tb;nfo(n)=sum(q);AMP1(n)=amp1;AMP2(n)=amp2;
    %   else
    %   end
end
% correct hydrophone 1 amplitudes for differences in receive amplifier gain
AMP1_corr1=AMP1*10^(gain_corr/20);
AMP1_corr2=zeros(size(AMP1));

% compute group velocities using non-zero arrival times
%  and intrinsic errors
for n=1:length(fp);
    fp(n)/1000
    

% get velocity to compute spreading losses for attenuation calculation
% if mean velocity is sensible use it, if not set to 1450 m/s    
    %if V(1,n)>0;velin=V(1,n);
    %else velloc=find(V(1,:)>0);velin=mean(V(1,velloc));   
    %end
    %if velin<50 || velin>3000; velin=1450;
    %elseif isnan(velin)==1;velin=1450; 
    %end
    
    velin=1400;
% use sub routine to compute dsoreading lossed, assume line
% perpendicular to source and through source centre intersects
% hydrophone centre
    [db_data,x,z]=spreading_noplot(velin,fp(n),radius);
    corr=10^((db_data(round(zrec2/(z(2)-z(1))))-db_data(round(zrec1/(z(2)-z(1)))))/20);
    %corr=1;
    %delta_recgain=rs2(n)-rs1(n);
    %gain_corr1=10^(delta_recgain/20);
    gain_corr1=1;
    AMP1_corr2(n)=gain_corr1*corr*AMP1_corr1(n); 
    alpha=8.686*log(AMP1_corr2(n)./AMP2(n))./(zrec2-zrec1);
    ALPHA(2,n)=alpha.*0.1;
    ALPHA(1,n)=8.686*log(AMP1_corr2(n)./AMP2(n))./(zrec2-zrec1);

end
   
 
figure(100);[out]=blockplot(fp/1000,V(1,:),V(2,:),block_fill);hold on;
set(gca,'fontsize',16);xlabel('frequency (kHz)');ylabel('velocity (m/s)');
title('compressional wave velocity with intrinsic errors');

figure(101);[out]=blockplot(fp/1000,ALPHA(1,:),ALPHA(2,:),block_fill);hold on;
set(gca,'fontsize',16);xlabel('frequency (kHz)');ylabel('attenuation (dB/m)');
title('compressional wave attenuation with intrinsic errors');

ALPHA=ALPHA';
V=V'; 
