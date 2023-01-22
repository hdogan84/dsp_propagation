function [ta,tb,amp1,amp2]=prop_proc_time(folder,fc,m,fs,nos_samples,b,a,ch1,ch2,ws,we,pl,waveform,device)
% Subroutine used by prop_nov_time to get amplitude and arrival time of signal aquired on two hydrophones
% NEED TO DEFINE DATA PATH ON LINE 28

% INPUTS: 
% folder: folder number where data stored
% fc=pump  frequency (Hz)
% m: shot number
% fs=sampling frequency (Hz)
% nos_samples=samples acquired per signal
% b and a: filter coefficients from prop_nov_time
% ch1: row in file in which contains data for first hydrophone
% ch2: row in file in which contains data for second hydrophone 
% ws and we: rough start and end of window applied (data points): use to
% remove possible cross talk signals
% waveform: envelope of tonal wave 1=square, 2 =hanning
% pl: pulse length (seconds) 
% device : pump source: use 1 for HF device, 2 for MF device and 3 for LF
% device

% OUTPUTS: 
% ta: arrival time of signal on hydrophone 1
% tb: arrival time of signal on hydrophone 2
% amp1: amplitude of signal on hydrophone 1
% amp2: amplitude of signal on hydrophone 2

try
    data=load(sprintf('C:/bubbles/110608_prop_mercury/%d/%dkHz%d',folder,fc/1000,m));
    [a1,b1]=size(data);if b1==nos_samples;data=data';else;end;
    
    % if signal not clipped carry on
    if (abs(max(data(:,ch1)))+abs(min(data(:,ch1))))<20 | (abs(max(data(:,ch2)))+abs(min(data(:,ch2))))<20
        prop1=data(ws:we,ch1);prop2=-1*data(ws:we,ch2);
        
        % filter data
        fprop1=filtfilt(b,a,prop1);fprop2=filtfilt(b,a,prop2);
        
        % get signal envelopes
        env1=abs(hilbert(fprop1));env1=env1/max(env1);
        
        % select portion of hydrophone 1 signal with amplitude >25 %
        % maximum: this is reference pulse
        loc=find(env1>0.25);
        refenv=env1(loc);
        tref=[0:1/fs:(length(refenv)-1)/fs];
        
        % correllate reference pulse with hydrophone 1 and 2 signals to get
        % arrival times
        cra=xcorr(abs(hilbert(fprop1)),refenv);
        ta=(find(cra==max(cra))-length(prop1))/fs;
        crb=xcorr(abs(hilbert(fprop2)),refenv);
        tb=(find(crb==max(crb))-length(prop1))/fs;
        
        env2=abs(hilbert(fprop2(round(tb*fs):round(tb*fs)+length(tref)-1)));
        norm_env2=env2/mean(env2);
        norm_refenv=refenv/mean(refenv);
        % get steady state regions to find amplitude over
        % for any square envelope (i.e. waveform=1) this is all of pulse barr first and last 6
        % oscillations, independent of device in use. 
        if waveform==1;  sp1=round(ta*fs+fs*300e-6);ep1=round(ta*fs+fs*pl-fs*200e-6);
                    sp2=round(tb*fs+fs*300e-6);ep2=round(tb*fs+fs*pl-fs*200e-6);
%         For hanning envelope this reions depends on the pump source in
%         use
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
        
        % get mean amplitude over this setady state region
        amp1=mean(abs(fprop1(sp1:ep1)));
        amp2=mean(abs(fprop2(sp2:ep2)));
        
        % plotting
        tplot=[ws/fs:1/fs:we/fs]';
        if m==1 %& rem(fc,10000)==0;
        figure;plot(tplot,fprop1);hold on;
        plot(tplot,fprop2,'r');
        hold on;
        plot(tplot,abs(hilbert(fprop1)),'k');hold on;
        plot(tplot,abs(hilbert(fprop2)),'r');
        xlabel('time (s)');ylabel('Amplitude (V)');
        title('filtered signals on hydrophone 1 (blue) and hydrophone 2 (red)');
    else
    end
    else
        % If clipped set amplitudes and arrival times to 0
        ta=0;tb=0;amp1=0;amp2=0;
        disp('clipped')
    end
    clear data
catch
        % If no data for that shot (i.e. acquisition timed out)
        % set amplitudes and arrival times to 0
               ta=0;tb=0;amp1=0;amp2=0;
          disp 'acquisition timeout OR error in processing '
 end