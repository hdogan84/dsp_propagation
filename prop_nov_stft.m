function [acoustics]=prop_nov_stft(fp,bw,spf,nos_samples,fs,folder,gain_corr,xrec1,xrec2,ch1,ch2,ws,we,waveform,block_fill);          %
% main processing function used to obtain phase velocity and amplitude of 
% signal on two receivers using short-time fourier transfor. 
% Analysis of water based data indocated that the time based correllation produced more reliable results

% !!!! need to adjust rec gains and path to data
 
% INPUTS: 
% fp=pump frequencies (Hz)
% bw as bandwidth in Hz
% spf=shots per frequency
% nos_samples=samples per shot 
% folder
% fs=sampling frequency (Hz)
% block fill

% VARIABLES THAT MAY NEED CHANGING
%gain_corr     % gain correction in dB between 2 receivers 
%xrec1
%xrec2
% channels for rec 1 and rec2 data ch1 and ch2 respectivelly
% start and end of rough window (in datapoints)
% waveform: i.e.  define pulse length 1= 1ms, 2 =0.5 ms, 3=20 osc, and 4=30
% osc. 
folder
% pulse length
if waveform==1;PL=1e-3*ones(size(fp));
elseif waveform==2;PL=0.5e-3*ones(size(fp));
elseif waveform==3;PL=(20./fp).*ones(size(fp));
elseif waveform==4;PL=(30./fp).*ones(size(fp));
else xxx
end

% RECEIVER GAINS
[rs1]=get_rec_sens(fp,34);
[rs2]=get_rec_sens(fp,38);

% STUFF
dx=xrec2-xrec1;
fny=fs/2;
t=[0:1/fs:(nos_samples-1)/fs];

% BLANK FILES
    cp=zeros(length(fp),spf);
    AMP1=zeros(length(fp),spf);
    AMP2=zeros(length(fp),spf);
    acoustics=zeros(4,length(fp));
    nfo=zeros(length(fp),1);
    
for n=1:length(fp);fc=fp(n);BW=bw(n);pl=PL(n);
    if rem(fc,10000)==0;
        fc
    else
    end
    fac=2;[N,Wn]=buttord([fc-fac*BW fc+fac*BW]/fny,[fc-2*fac*BW fc+2*fac*BW]/fny,0.5,5);
    [b,a]=butter(N,Wn);
    for m=1:spf;
        [cpi,amp1,amp2]=prop_process_stft(folder,fc,m,fs,b,a,dx,ch1,ch2,ws,we,pl);
        cp(n,m)=cpi;
        AMP1(n,m)=amp1;
        AMP2(n,m)=amp2;
    end
end
cp
 for n=1:length(fp);fp(n)
 LOC1=find(cp(n,:)>0 & cp(n,:)<1e4);
 nfo(n)=length(LOC1);
 acoustics(1,n)=mean(cp(n,LOC1));
 acoustics(2,n)=std(cp(n,LOC1));
 [db_data,z]=spreading(acoustics(1,n),fp(n));%
 delta_recgain=rs2(n)-rs1(n);
 AMP1_corr=AMP1(n,:)*10^(gain_corr/20)*10^(delta_recgain/20)*10^((db_data(find(z==xrec2))-db_data(find(z==xrec1)))/20); 
 alpha(n,:)=8.686*log(AMP1_corr./AMP2(n,:))./(dx);
 acoustics(3,n)=mean(alpha(n,LOC1));
 acoustics(4,n)=std(alpha(n,LOC1));

 end

 figure(5);[out]=blockplot(fp/1000,acoustics(1,:),acoustics(2,:),block_fill);hold on;
 set(gca,'fontsize',16);xlabel('frequency (kHz)');ylabel('Phase velocity (m/s)');
 figure(6);[out]=blockplot(fp/1000,acoustics(3,:),acoustics(4,:),block_fill);hold on;
 set(gca,'fontsize',16);xlabel('frequency (kHz)');ylabel('Attenuation coeff. (dB/m)');

% %save C:\bubbles\oct_2007_combfreqresults\folder32 f1amp f2amp sumamp
% %diffamp f1doubamp f2doubamp -V6