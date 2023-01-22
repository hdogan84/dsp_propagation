function [V,ALPHA,nfo]=prop_nov_time(fp,bw,spf,nos_samples,fs,folder,zrec1,zrec2,rec1no,rec2no,gain_corr,ch1,ch2,ws,we,waveform,block_fill,device);          %
addpath 'K:\HARD_DRIVES\2008\110608_prop_mercury'
% main function used to process propagation data from the water tank, using
% signals received at two hydrophones
% This uses a time based correlation to detect arrival times (and so group velocity)
% and the amplitudes of the steady state portion of the signals on the two
% hydrophones to commpute attenuation. 
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
% rec2no: last two digits of serial number of second hydrophone
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
% std (row 2) and intrinsic error (row 3)
% ALPHA: attenuation (dB/m) including mean value for all shots (row 1)
% upward std error (row 2), downward std error (row 3) and intrinsic error (row 4)
% nfo: numbersof shots that were successfullt processed at each frequency

% matricies required
nfo=zeros(length(fp),1);
V=zeros(3,length(fp));
ALPHA=zeros(4,length(fp));

% select pump device 
if device==1;radius=0.113/2;        % HF pump
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

% receive sensitivities of hydrophones
[rs1]=get_rec_sens(fp,rec1no);
[rs2]=get_rec_sens(fp,rec2no);

% stuff
fny=fs/2;
te=1/fs;
t=[0:1/fs:(nos_samples-1)/fs];

% generate blank files
    toa=zeros(length(fp),spf);
    tob=zeros(length(fp),spf);
    AMP1=zeros(length(fp),spf);
    AMP2=zeros(length(fp),spf);
    ae=zeros(length(fp),spf);
    
for n=1:length(fp);fc=fp(n);BW=bw(n);fc
    
    % design filters
    fac=2;[N,Wn]=buttord([fc-fac*BW fc+fac*BW]/fny,[fc-2*fac*BW fc+2*fac*BW]/fny,0.5,5);
    [b,a]=butter(N,Wn);
    
    % for each shot use sub routine prop_proc_time to get amplitude and
    % arrival time of signal on each hydrophone
    for m=1:spf;
        [ta,tb,amp1,amp2]=prop_proc_time(folder,fc,m,fs,nos_samples,b,a,ch1,ch2,ws,we,PL(n),waveform,device);
       if ta>0; toa(n,m)=ta;tob(n,m)=tb;AMP1(n,m)=amp1;AMP2(n,m)=amp2;
       else
       end
    end
end

% correct hydrophone 1 amplitudes for differences in receive amplifier gain
AMP1_corr1=AMP1*10^(gain_corr/20);
AMP1_corr2=zeros(size(AMP1));

% compute group velocities using non-zero arrival times
% also compute std in velocities and intrinsic errors
for n=1:length(fp);
    fp(n)/1000
    ftob=tob(n,:);
    gftob=ftob(find(ftob));
    ftoa=toa(n,:);
    gftoa=ftoa(find(ftob));
    invdt=1./(gftob-gftoa);
    V(1,n)=(zrec2-zrec1).*mean(invdt,2);
    V(2,n)=(zrec2-zrec1).*std(invdt,0,2);
    pex=0.005/(zrec2-zrec1);
    peinvdt=mean(5e-7./(gftob-gftoa),2);
    pev=sqrt(pex^2+peinvdt.^2);
    vout=V(1,n);
    V(3,n)=pev*vout;
    LOC1=find(toa(n,:)>0);
    nfo(n)=length(LOC1);

    % get velocity to compute spreading losses for attenuation calculation
    % if mean velocity is sensible use it, if not set to 1450 m/s
    if V(1,n)>0;velin=V(1,n);
    else velloc=find(V(1,:)>0);velin=mean(V(1,velloc));   
    end
    if velin<0 || velin>10000; velin=1450;
    elseif isnan(velin)==1;velin=1450 
    end

    % use sub routine to compute dsoreading lossed, assume line
    % perpendicular to source and through source centre intersects
    % hydrophone centre
    [db_data,x,z]=spreading_noplot(velin,fp(n),radius);    
    corr=10^((db_data(round(zrec2/(z(2)-z(1))))-db_data(round(zrec1/(z(2)-z(1)))))/20);
    
    % correct hydrophone 1 amplitude for receiver sensitivity and spreading losses
    delta_recgain=rs2(n)-rs1(n);
    gain_corr1=10^(delta_recgain/20);
    AMP1_corr2(n,:)=gain_corr1*corr*AMP1_corr1(n,:); 
    
    % compute mean attenuation, intrinsic errors and asymmetric std errors
    fed=0.005/(zrec2-zrec1);
    feamp1=(1.22e-3)./AMP1(n,LOC1);
    feamp2=(1.22e-3)./AMP2(n,LOC1);
    fear=sqrt(feamp1.^2+feamp2.^2);
    feln=fear./log(AMP1_corr2(n,LOC1)./AMP2(n,LOC1));
    fealpha=8.686*sqrt(fed^2+feln.^2);
    alpha=8.686*log(AMP1_corr2(n,LOC1)./AMP2(n,LOC1))./(zrec2-zrec1);
    ae(n,LOC1)=alpha.*fealpha;
    ALPHA(4,n)=mean(ae(n,LOC1));
    ALPHA(1,n)=8.686*log(mean(AMP1_corr2(n,LOC1))./mean(AMP2(n,LOC1)))./(zrec2-zrec1);
    A1std=std(AMP1_corr2(n,LOC1));
    A2std=std(AMP2(n,LOC1));
    ALPHA(2,n)=8.686*log((mean(AMP1_corr2(n,LOC1))+A1std)./(mean(AMP2(n,LOC1))-A2std))./(zrec2-zrec1)-ALPHA(1,n);
    ALPHA(3,n)=8.686*log((mean(AMP1_corr2(n,LOC1)-A1std))./(mean(AMP2(n,LOC1)+A2std)))./(zrec2-zrec1)-ALPHA(1,n);   
end

    %PLOTTING
figure(102);[out]=blockplot(fp/1000,V(1,:),V(2,:),block_fill);hold on;
set(gca,'fontsize',16);xlabel('frequency (kHz)');ylabel('velocity (m/s)');
title('compressional wave velocity with variability (sd) errors');
figure(103);[out]=blockplot(fp/1000,V(1,:),V(3,:),block_fill);hold on;
set(gca,'fontsize',16);xlabel('frequency (kHz)');ylabel('velocity (m/s)');
title('compressional wave velocity with intrinsic errors');
figure(104);[out]=blockplot_assyerrors(fp/1000,ALPHA(1,:),ALPHA(2,:),ALPHA(3,:),block_fill);
set(gca,'fontsize',16);xlabel('frequency (kHz)');ylabel('attenuation (dB/m)');
title('compressional wave attenuation with variability (sd) errors');
figure(105);[out]=blockplot(fp/1000,ALPHA(1,:),ALPHA(4,:),block_fill);hold on;
set(gca,'fontsize',16);xlabel('frequency (kHz)');ylabel('attenuation (dB/m)');
title('compressional wave attenuation with intrinsic errors');