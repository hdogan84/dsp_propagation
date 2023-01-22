function [cpi,amp1,amp2]=prop_process_stft(folder,fc,m,fs,b,a,dx,ch1,ch2,ws,we,pl);
nfft=16384;
try
    data=load(sprintf('E:/acoustic_data/feb08_water_prop/lf/%d/%dkHz%d',folder,fc/1000,m));
    prop1=data(ws:we,ch1);prop2=data(ws:we,ch2);
    if (abs(max(prop1))+abs(min(prop1)))<18;
        %fprop1=prop1;fprop2=prop2;
    fprop1=filtfilt(b,a,prop1);fprop2=filtfilt(b,a,prop2);
    win=window(@blackman,round(fs*1.15*pl));
    [s1,f,t,p1]=spectrogram(fprop1,win,round(0.94*length(win)),nfft,2e6);
    %figure;spectrogram(fprop1,win,round(0.94*length(win)),nfft,2e6)
    %figure;surf(t,f,10*log10(abs(s1)),'EdgeColor','none');%VIEW(180,0)
    [s2,f,t,p2]=spectrogram(fprop2,win,round(0.94*length(win)),nfft,2e6);
    %figure;spectrogram(fprop2,win,round(0.94*length(win)),nfft,2e6)
    %figure;surf(t,f,10*log10(abs(s2)),'EdgeColor','none');%VIEW(180,0)
    floc=min(find(f>fc));
    amp1=max(abs(s1(floc,:)));
    amp2=max(abs(s2(floc,:)));
    t1=t(find(abs(s1(floc,:))==amp1));
    t2=t(find(abs(s2(floc,:))==amp2));
    cpi=dx/(t2-t1);
    else
        cpi=0;amp1=0;amp2=0;
        disp('clipped')
    end
    clear data
   catch
           cpi=0;amp1=0;amp2=0;
           disp 'timeout occured during acquisition'
   end