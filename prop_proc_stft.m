function [cpi,amp1,amp2]=prop_process_stft(folder,fc,m,fs,b,a,dx,ch1,ch2,ws,we,pl);
% processing function used by prop_nov_stft to obtain phase velocity and amplitude of 
% signal on two receivers using short-time fourier transfor. 
% Analysis of water based data indocated that the time based correllation produced more reliable results
nfft=16384;
try
    load(sprintf('E:/acoustic_data/feb08_water_prop/%d/%dkHz%d',folder,fc/1000,m));
    prop1=data(ws:we,ch1);prop2=data(ws:we,ch2);
    shit
    if (abs(max(prop1))+abs(min(prop1)))<16;
        fprop1=prop1;fprop2=prop2;
    %fprop1=filtfilt(b,a,prop1);fprop2=filtfilt(b,a,prop2);
    win=window(@blackman,fs*1.15*pl);
    [s1,f,t,p1]=spectrogram(fprop1,win,round(0.94*length(win)),nfft,2e6);
    %figure;surf(t,f,10*log10(abs(p1)),'EdgeColor','none');VIEW(180,0)
    [s2,f,t,p2]=spectrogram(fprop2,win,round(0.94*length(win)),nfft,2e6);
    %figure;surf(t,f,10*log10(abs(p2)),'EdgeColor','none');VIEW(180,0)
    floc=min(find(f>fc));
    amp1=max(p1(floc,:));
    amp2=max(p2(floc,:));
    t1=t(find(s1(floc,:)==amp1));
    t2=t(find(s2(floc,:)==amp2));
    cpi=dx/(t1-t2);
    else
        cpi=0;amp1=0;amp2=0;
        disp('clipped')
    end
    clear data
   catch
           cpi=0;amp1=0;amp2=0;
           disp 'timeout occured during acquisition'
   end