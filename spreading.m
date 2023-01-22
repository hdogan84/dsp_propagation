function [db_data,x,z]=spreading(c,freq,radius);
% subroutine called by getvpp_cf_focus
% modified version of spreading programs developes by Matt Simpson (2004)
% calclates spreading patern along prependicular 1 m line through centre of a
% piston source
% uses continuous wave: discrepancies between these results and those using
% tonal pulse are negligible

% OUTPUTS: 
% db_data: beam pattern expressed in dB relative to the maximum value along
% the line
% x: x-corodinate of line which equals zero (i.e. centre of source)
% z: points along line at which spreading is computed
% 
% INPUTS: 
% c: sound speed in medium of interest (m/s)
% freq: frequency being emitted by piston source
% radius: active radius of piston face

% model parameters
lambda=c/freq;
x=[0];
y=[0];
z=[0.001:0.001:1.00];       % can modify if larger distances are required

% calculate pressure ate each point along line using subroutine
% ht_circle_prop
countx=1;county=1;
for countz=1:length(z);
    [h_xt, time1(countx,county,countz), time3(countx,county,countz), time_vec]=ht_circle_proc(x(countx),y(county),z(countz),freq,radius,c);
    sig=exp(-j*((2*pi*freq.*time_vec)));
    r(countx,countz,county)=trapz((sig.*h_xt))*(time_vec(end)-time_vec(end-1));
end
x1=[rot90(rot90(-x)) x(2:end)];         % Symmetrical x axis
r2=flipdim(r(2:end, :, :), 1);  
data=[r2;r];                            % Symmetrical data

% convert predicted pressures to dB relative to maximum pressure along the
% line
db_data=-20*log10(max(max(abs(data)))./abs(data));

% plotting
figure;plot(z,db_data);