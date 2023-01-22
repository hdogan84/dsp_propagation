function [fc_mean,fc_std,bw_mean,bw_std]=compute_fcbw(fs,fin,spf,folder,ws,ch,wf);
% computes central frequency and bandwidths of a series pulses, 
% which are assumed to consist of
% spf shots at a number of tonal frequencies defines in fin
% NEED to adjust data path on line 40
% 
% INPUTS: 
% fs: sampling frequency (Hz)
% fin: frequencies of tonal pulses (Hz)
% spf: number of shots at each tonal frequency
% folder: folder in data path defined on line to examine
% ws: start of steady state portion of signal to analysis (in data points)
% ch: row of data file to analyse
% wf: defines waveform length used for tonal pulses.
% 
% OUTPUTS:
% fc_mean: mean central frequency across spf shots at each tonal frequency
% fc_std: std in central frequencies across spf shots at each tonal frequency
% bw_mean: mean bandwidth across spf shots at each tonal frequency
% bw_std: std in bandwidths across spf shots at each tonal frequency
% 


close all
fny=fs/2;
% blank matrices
fc=zeros(spf,length(fin));
bw=zeros(spf,length(fin));
if wf==1;len=1.1e-3*fs*ones(1,length(fin));
elseif wf==2;len=0.55e-3*fs*ones(1,length(fin));
elseif wf==3;len=22*round(fs./fin);
elseif wf==4;len=32*round(fs./fin);
elseif wf==5;len=2.2e-3*fs*ones(1,length(fin));
else xxx
end
    
for n=1:length(fin);n
    for m=1:spf;
        try
    data1=load(sprintf('E:/acoustic_data/feb08_water_prop/mf//%d/%dkHz%d',folder,fin(n)/1000,m));
    prop1=data1(ws:ws+len(n),ch);
    [pz,f]=psd(prop1,8192,fs);
    pz(1:20)=0;
    fc(m,n)=f(find(pz==max(pz)));
    bw(m,n)=(max(f(find(pz>(max(pz)/2))))-min(f(find(pz>(max(pz)/2)))));
    if m==1;figure;plot(f/1000,pz/max(pz),'r.');xlabel('frequency (kHz)');ylabel('psd');set(gca,'xlim',[0 200]);else;end
        catch
        end
    end;
end;
fc_mean=zeros(size(fin));
fc_std=zeros(size(fin));
bw_mean=zeros(size(fin));
bw_std=zeros(size(fin));
for n=1:length(fin);
loc=find(fc(:,n)~=0)
fc_mean(n)=mean(fc(loc,n))
fc_std(n)=std(fc(loc,n))
bw_mean(n)=mean(bw(loc,n))
bw_std(n)=std(bw(loc,n))
end
