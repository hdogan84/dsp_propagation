function [ta,tb,q,amp1,amp2]=prop_proc_time_stack(folder,fc,spf,fs,nos_samples,b,a,ch1,ch2,ws,we,pl,waveform,device);
% Subroutine used by prop_nov_time_stacked to get amplitude and arrival time of signal aquired on two hydrophones
% NEED TO DEFINE DATA PATH ON LINE 30

% INPUTS: 
% folder: folder number where data stored
% fc=pump  frequency (Hz)
% spf: shots per frequnecy
% fs=sampling frequency (Hz)
% nos_samples=samples acquired per signal
% b and a: filter coefficients from prop_nov_time
% ch1: row in file in which contains data for first hydrophone
% ch2: row in file in which contains data for second hydrophone 
% ws and we: rough start and end of window applied (data points): use to
% remove possible cross talk signals
% waveform: length of tonal pulse sent to the pump amplifier
% length 1= 1ms, 2 =0.5 ms, 3=20 osc, and 4=30 osc. 
% device : pump source: use 1 for HF device, 2 for MF device and 3 for LF
% device

% OUTPUTS: 
% ta: arrival time of signal on hydrophone 1
% tb: arrival time of signal on hydrophone 2
% q: matrix with 0 for clipped or timed out signals and 1 for good signals
% amp1: amplitude of signal on hydrophone 1
% amp2: amplitude of signal on hydrophone 2
prop1=zeros(spf,we-ws+1);prop2=zeros(spf,we-ws+1);fprop1=zeros(spf,we-ws+1);fprop2=zeros(spf,we-ws+1);
q=ones(1,spf);
for m=1:spf;
    filename=sprintf('C:/Users/hakan.dogan/Desktop/InSitu/2008/150408_prop_mercury/%d/%dkHz%d',folder,fc/1000,m);
    if isfile(filename)==1
      data=load(filename);
      [a1,b1]=size(data);
      data=data';
      % if signal not clipped carry on
      if (abs(max(data(:,ch1)))+abs(min(data(:,ch1))))<20 || (abs(max(data(:,ch2)))+abs(min(data(:,ch2))))<20
        prop1(m,:)=data(ws:we,ch1);prop2(m,:)=data(ws:we,ch2);
        prop1(m,:)-=mean(prop1(m,:)); prop2(m,:)-=mean(prop2(m,:));
        % filter data
        fprop1(m,:)=filtfilt(b,a,prop1(m,:));fprop2(m,:)=filtfilt(b,a,prop2(m,:));
        q(m)=1;
    %clear data
      else
        % If clipped
        disp('clipped');q(m)=0;
      end
    else
       q(m)=0;
    end

    % If no data for that shot (i.e. acquisition timed out)      
  %disp 'acquisition timeout OR error in processing ';      
end


% find shots with non-zero values i.e. no timeout and no clipping
fprop1=fprop1(q>0,:); fprop2=fprop2(q>0,:);
testsum1=sum(abs(fprop1'));testsum2=sum(abs(fprop2'));
loc=find(testsum1>0 & testsum2>0);


if length(loc)>0
    % select non-zero shots
    fprop1keep=fprop1(loc,:);fprop2keep=fprop2(loc,:);
    % stack using median
    stacked1=median(fprop1keep);stacked2=median(fprop2keep);
    % select portion of envelope of hydrophone 1 signal and correllate with
    % hyrdophone 1 and 2 signals
     env1=abs(hilbert(stacked1));env1norm=env1/max(env1);
     env2=abs(hilbert(stacked2));
     loc=find(env1norm>0.25);
     refenv=env1norm(loc);
     tref=[0:1/fs:(length(refenv)-1)/fs];
     cra=xcorr(env1,refenv);
     ta=(find(cra==max(cra))-length(prop1))/fs;
     crb=xcorr(env2,refenv);
     tb=(find(crb==max(crb))-length(prop1))/fs;
     % get steady state regions to find amplitude over
     % for sqy=uare envelope this is all of pulse barr first and last 6
     % oscillations, independent of device in use. 
    if waveform==1;  sp1=round(ta*fs+fs*300e-6);ep1=round(ta*fs+fs*pl-fs*200e-6);
                     sp2=round(tb*fs+fs*300e-6);ep2=round(tb*fs+fs*pl-fs*200e-6);
                     %sp1=round(ta*fs+fs*150e-6);ep1=round(ta*fs-fs*100e-6+fs*pl);
                     %sp2=round(tb*fs+fs*150e-6);ep2=round(tb*fs-fs*100e-6+fs*pl);
     % For hanning envelope this reions depends on the pump source in use 
    elseif waveform==3 && device==1;
         sp1=round(ta*fs+5*fs/fc);ep1=round(ta*fs+fs*pl-5*fs/fc);
         sp2=round(tb*fs+5*fs/fc);ep2=round(tb*fs+fs*pl-5*fs/fc);
    elseif waveform==3 && device==2;
         sp1=round(ta*fs+4*fs/fc);ep1=round(ta*fs+fs*pl-4*fs/fc);
         sp2=round(tb*fs+4*fs/fc);ep2=round(tb*fs+fs*pl-4*fs/fc);            
    elseif waveform==3 && device==3;
         sp1=round(ta*fs+3*fs/fc);ep1=round(ta*fs+fs*pl-3*fs/fc);
         sp2=round(tb*fs+3*fs/fc);ep2=round(tb*fs+fs*pl-3*fs/fc);
    else  xxxx
    end
    % sp1=800;ep1=4200;sp2=500;ep2=5500;
    % get mean amplitude over this setady state region
    if sp1>0 && sp2>0;
     amp1=mean(abs(stacked1(sp1:ep1)));
     amp2=mean(abs(stacked2(sp2:ep2)));
    else ta=0;tb=0;amp1=0;amp2=0;
    end
     
else ta=0;tb=0;amp1=0;amp2=0;
end

% plot filtered stacked data 
tplot=[ws/fs:1/fs:we/fs]';
%figure;subplot(211);plot(stacked1);hold on;plot(stacked2,'r');
%xlabel('time (s)');ylabel('Amplitude (V)');
%title('filtered stacked data for first (blue) and second (red) hydrophone');
%subplot(212);plot(abs(hilbert(stacked1/max(abs(stacked1)))),'k');hold on;
%plot(abs(hilbert(stacked2/max(abs(stacked2)))),'g');
%xlabel('time (s)');ylabel('Amplitude (V)');
%title('normalised envelopes of filtered stacked data for first (blue) and second (red) hydrophone');
