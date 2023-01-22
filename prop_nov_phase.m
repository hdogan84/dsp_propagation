function [acoustics]=prop_nov_phase(fp,bw,spf,nos_samples,fs,folder,gain_corr,xrec1,xrec2,ch1,ch2,ws,we,block_fill);          %
% Main function  which uses speacral based approach to
% compute group velocity, phase velocity and attenuation of  propagation data between
% two recivers. outout acoustics contains phase velocity and std error in rows 1 and 2, 
% group velocity and std error in rows 3 and 4 and attenuation and std error in rows
% 5 and 6
% tests on water based data indicate that time based correllation produces
% more reliable results

% !!!! need to adjust control and samples paths!!!!!!
 
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

% RECEIVER GAINS
[rs1]=get_rec_sens(fp,34);
[rs2]=get_rec_sens(fp,38);

% STUFF
dx=xrec2-xrec1;
fny=fs/2;
t=[0:1/fs:(nos_samples-1)/fs];

% BLANK FILES
    cg=zeros(length(fp),spf);
    cp=zeros(length(fp),spf);
    AMP1=zeros(length(fp),spf);
    AMP2=zeros(length(fp),spf);
    acoustics=zeros(6,length(fp));
    nfo=zeros(length(fp),1);
    %ae=zeros(length(fp),spf);
    
for n=1:length(fp);fc=fp(n);BW=bw(n);fc
    fac=2;[N,Wn]=buttord([fc-fac*BW fc+fac*BW]/fny,[fc-2*fac*BW fc+2*fac*BW]/fny,0.5,5);
    [b,a]=butter(N,Wn);
    for m=1:spf;
        [cgi,cpi,amp1,amp2]=prop_proc_phase(folder,fc,m,fs,b,a,dx,ch1,ch2,ws,we);
        cg(n,m)=cgi;
        cp(n,m)=cpi;
        AMP1(n,m)=amp1;
        AMP2(n,m)=amp2;
    end
end
 for n=1:length(fp);fp(n)
 LOC1=find(cg(n,:)>0);
 nfo(n)=length(LOC1);
 acoustics(1,n)=mean(cp(n,LOC1));
 acoustics(2,n)=std(cp(n,LOC1));
 acoustics(3,n)=mean(cg(n,LOC1));
 acoustics(4,n)=std(cg(n,LOC1));
 [db_data,z]=spreading(mean(cg(n,LOC1)),fp(n));
 delta_recgain=rs2(n)-rs1(n);
 AMP1_corr=AMP1(n,:)*10^(gain_corr/20)*10^(delta_recgain/20)*10^((db_data(find(z==xrec2))-db_data(find(z==xrec1)))/20); alpha(n,:)=8.686*log(AMP1_corr./AMP2(n,:))./(dx);
 acoustics(5,n)=mean(alpha(n,LOC1));
 acoustics(6,n)=std(alpha(n,LOC1));
 
 end

 figure(101);[out]=blockplot(fp/1000,acoustics(1,:),acoustics(2,:),block_fill);hold on;
 set(gca,'fontsize',16);xlabel('frequency (kHz)');ylabel('Phase velocity (m/s)');
 figure(102);[out]=blockplot(fp/1000,acoustics(3,:),acoustics(4,:),block_fill);hold on;
 set(gca,'fontsize',16);xlabel('frequency (kHz)');ylabel('Group velocity (m/s)');
 figure(103);[out]=blockplot(fp/1000,acoustics(5,:),acoustics(6,:),block_fill);hold on;
 set(gca,'fontsize',16);xlabel('frequency (kHz)');ylabel('Attenuation coeff. (dB/m)');