function [nos_file_out,v,ferr_v,att,ferr_att]=data_analysis_aug(fin,bw,spf,nos_samples,startpt,endpt,xbox,dw,vw,ds,fs,gain_corr);          %
% OLD ANALYSIS PROPGRAM USED TO EXAMINE PROPAGATION DATA THROUGH SAMPLE 
% BOX

% !!!! need to adjsut control and samples paths!!!!!!
% !!!! need to adjust 
% process lab data using scenario 2, namely comparison of no box and 
% sample box data, using predicted Transmission losses and time delay incur
% by sample box

% inputs include: 
% fc=analysis central frequencies (Hz)
% bw=bandwidths (Hz)
% spf=shots per frequency
% nos_samples=samples per shot 
% startpt=data point at which to start signals
% endpt=data point at which to terminate signals
% xbox=samples thickness (m)
% dw=density of water (kg/m3)
% vw=velocity water (kg/m3)
% ds=density of sample (kg/m3)
% fs=sampling frequency (Hz)
% gain_corr: gain correction for adjustment of rec gains (multip factor)

% additional fixed parameters
dp=1200;    % density perspex
vp=2700;    % velocity perspex

% matricies required
nos_file_out=zeros(length(fin),1);
v=zeros(1,length(fin));
dt=zeros(1,length(fin));
Ac=zeros(1,length(fin));
As=zeros(1,length(fin));
ferr_v=zeros(1,length(fin));
ferr_att=zeros(1,length(fin));
fny=fs/2;
te=1/fs;
t=[0:1/fs:(nos_samples-1)/fs];


for n=1:length(fin);
    f=fin(n);
    f
     Z=zeros(spf,nos_samples);Y=zeros(size(Z));
    if rem(n,1)==0;
    figure;%set(gcf,'units','normalized');set(gcf,'position',[.1,.1,.8,.8]);
    else 
    end
    for m=1:spf;
       try
    load(sprintf('E:/prop_aug/6/%dkHz%d',f/1000,m))
    Z(m,:)=data(:,2);
    clear data
       catch
        disp 'no such files'
        end
        try
    load(sprintf('E:/prop_aug/1/%dkHz%d',f/1000,m))
    Y(m,:)=data(:,2);
    clear data
    catch
        disp 'no such files'
        end
    end
    
    % remove simultaneous non-existant files
    test1=sum(abs(Y),2); 
    test2=sum(abs(Z),2);
    loc=find(test1>0 & test2>0);
    nos_files_out(n)=length(loc);
    Y1=Y(loc,:);
    Z1=Z(loc,:);
    clear Y Z
    
       % stacking
    
    Y_meanraw=mean(Y1,1);
    Y_mean=Y_meanraw;
    Y_mean(1:startpt-1)=0;        
    Y_mean(endpt+1:nos_samples)=0;
    Y_std=std(Y1,1);
    Z_meanraw=mean(Z1,1);
    Z_mean=Z_meanraw;
    Z_mean(1:startpt-1)=0;        
    Z_mean(endpt+1:nos_samples)=0;    
    Z_std=std(Z1,1);
clear Y Z
    
    % do filter
    BW=bw(n);
    f=fin(n);
    [N,Wn]=buttord([f-0.7*BW f+0.7*BW]/fny,[f-1.4*BW f+1.4*BW]/fny,0.5,5);
    [b,a]=butter(N,Wn);
    Y2=filtfilt(b,a,Y_mean);Yo=Y2/max(abs(Y2));
    Z2=filtfilt(b,a,Z_mean);Zo=Z2/max(abs(Z2));
    if rem(n,5)==0;figure;freqz(b,a,1024);else;end
    
    %correlation to times and errors in velocity  
    envy=abs(hilbert(Yo));
    envz=abs(hilbert(Zo));
    water=envy(find(envy>(max(envy)/2)));
    corr1=xcorr(envy,water);
    corr2=xcorr(envz,water);
    pt1=find(corr1==max(corr1));t1=(pt1*(1/fs));
    pt2=find(corr2==max(corr2));
    t2=pt2*(1/fs);dt1=(t2-t1);
    dt(n)=dt1; 
    ferr_dem=sqrt(te.^2+((1e-3)/vw)^2)/(dt1+((96e-3)/vw));
    ferr_v1=sqrt((1e-3/xbox)^2+ferr_dem^2);
    ferr_v(n)=ferr_v1;
    % amplitudes and errors in attenuation 
    Ac(n)=mean(abs(Y2(min(find(envy>(max(envy)/2))):max(find(envy>(max(envy)/2))))));
    As(n)=mean(abs(Z2([(pt2-length(Z2)):(pt2-length(Z2)+length(water)-1)])));
    ferras=0.005/As(n);
    ferrac=0.005/Ac(n);
    ferr_asac=sqrt(ferras^2+ferrac^2);
    ferr_ln=ferr_asac;
    ferratt=sqrt((1e-3/xbox)^2+ferr_ln^2);
    ferr_att(n)=ferratt;

    %plot time history
    if rem(n,1)==0;
     subplot(121);
     plot(t*1e6,Z_meanraw,'k');hold on;plot(t*1e6,Y_meanraw,'r'); 
     xlabel('time (us)');ylabel('V');
     title([' raw signals at ' int2str(fin(n)/1000) 'kHz']);
     
     subplot(122);
     %plot(t,envz,'k');hold on;plot(t,envy,'r');
     plot(t*1e6,Zo,'k');hold on;plot(t*1e6,Yo,'r');     
     plot(1*[(pt2-length(Zo))/fs:1/fs:(pt2-length(Zo)+length(water)-1)/fs]',water,'m','linewidth',2);
     xlabel('time (us)');ylabel('V');
     title(['Normalised filtered data at ' int2str(fin(n)/1000) 'kHz']);
    else;end
    
%     subplot(224);
%     startpt=ceil(winstart*fs+1);
%     endpt=ceil(winstart*fs+1+winlength*fs/(fin(n)));
%     Z_win=Z_mean(startpt:endpt);
%     [pz,f]=psd(Z_win,4096,fs);
%     pz(1:5)=0;
%     plot(f/1000,pz/max(pz),'k');hold on;
%     Y_win=Y_mean(startpt:endpt);
%     [py,f]=psd(Y_win,4096,fs);
%     py(1:5)=0;
%     plot(f/1000,py/max(py),'r');
%     set(gca,'xlim',[0 150]);
            
end;

% velocities
v=xbox./(dt+ones(size(dt))*xbox/vw);
disp('mean velocity is...');
vs=mean(v)

% transmission coef
zw=dw*vw;
zs=ds*vs;
zp=dp*vp;
%Rwp=(zp-zw)/(zp+zw);Twp=sqrt(1-Rwp^2);
%Rps=(zs-zp)/(zp+zs);Tps=sqrt(1-Rps^2);  
Rws=(zs-zw)/(zw+zs);Tws=sqrt(1-Rws^2)
%factor=1/((Twp^2)*(Tps^2));    % factor if acoustic waves see wall
factor=1/Tws^2
att=(-20/xbox)*log10(factor*(As/gain_corr)./Ac);
    
    %plot final results
    figure;plot(fin/1e3,v,'ro','markersize',5,'linewidth',2);
    hold on;
for n=1:length(fin);
        vup=v(n)*(1+ferr_v(n));
        vdown=v(n)*(1-ferr_v(n));
        plot([fin(n)/1e3,fin(n)/1e3],[vup,vdown],'b','linewidth',2);
end
set(gca,'fontsize',14);title('velocity');
xlabel('frequency (kHz)');
ylabel('velocity (m.s^-^1)')

    figure;
plot(fin/1e3,att,'ro');
hold on;
set(gca,'fontsize', 16);
xlabel('Frequency (kHz)');
ylabel('attenuation (dB.m^-^1)');
for n=1:length(fin);
        aup=att(n)*(1+ferr_att(n)/100);
        adown=att(n)*(1-ferr_att(n)/100);
        plot([fin(n)/1e3,fin(n)/1e3],[aup,adown],'b','linewidth',2);
end

figure;plot(fin/1000,dt*1e6,'ro','markersize',5,'linewidth',2);ylabel('time difference (us)');
figure;plot(fin/1000,As./Ac,'ro');title('amplitude ratio');
