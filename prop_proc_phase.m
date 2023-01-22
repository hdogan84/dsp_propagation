function [cgi,cpi,amp1,amp2]=prop_proc_phase(folder,fc,m,fs,b,a,dx,ch1,ch2,ws,we)
% function used by prop_nov_phase which uses speacral based approach to
% compute group (cgi) and phase (cpi) velocity of propagation data between
% two recivers. Also compute amplitude of steady state signal on receiver 1
% (amp1) and receiver 2 (amp2) which prop_nov_phase converts into
% attenuations
% tests on water based data indicate that time based correllation produces
% more reliable results


nfft=16384;
try
    data=load(sprintf('E:/acoustic_data/feb08_water_prop/lf/%d/%dkHz%d',folder,fc/1000,m));
    prop1=data(ws:we,ch1);prop2=data(ws:we,ch2);
    dt=0;
    %figure;plot(prop1);hold on;plot(prop2,'r')
    if (abs(max(prop1))+abs(min(prop1)))<16;
    %fprop1=prop1;fprop2=prop2;
    fprop1=filtfilt(b,a,prop1);fprop2=filtfilt(b,a,prop2);
    x=fft(fprop1,nfft);x=x(1:round(nfft/10));
    y=fft(fprop2,nfft);y=y(1:round(nfft/10)); 
    s1=unwrap(angle(x));
    s2=unwrap(angle(y));
    s3=unwrap(angle(x.*conj(y)));
    ds=s1-s2;
    fplot=[fs/nfft:fs/nfft:fs/10];
    pt1=min(find(fplot>fc-2000));
    pt2=min(find(fplot>fc+2000));
    ftest=fplot(pt1:pt2)';
    stest1=s3(pt1:pt2);
    %[p]=polyfit(ftest,stest1,1);
    varx1=ftest-mean(ftest);vary1=stest1-mean(stest1);
    ssx1=sum(varx1.^2);ssy1=sum(vary1.^2);
    spxy1=sum(varx1.*vary1);m1=spxy1/ssx1;
    
%      figure;subplot(221);plot(fplot/1000,s1,'r.');
%      subplot(222);plot(fplot/1000,s2,'b.');
%      subplot(223);plot(fplot/1000,ds,'.');
%      subplot(224);plot(fplot/1000,s3,'.');
%      figure;;plot(fplot/1000,abs(x),'r.');hold on;plot(fplot/1000,abs(y),'b.');
    cgi=dx/(m1/(2*pi)+dt);
    %dphi1=ds(max(find(fplot<fc)));c2=dx/((dphi1-pi)/(2*pi*fc)+dt)
    fpos=min(find(fplot>fc));
    dphi2=s3(fpos);
    cpi=dx/((dphi2-pi)/(2*pi*fc)+dt);
    amp1=abs(x(fpos));amp2=abs(y(fpos));
    else
        cpi=0;cgi=0;amp1=0;amp2=0;
        disp('clipped')
    end
    clear data
   catch
           cpi=0;cgi=0;amp1=0;amp2=0;
           disp 'timeout occured during acquisition'
   end