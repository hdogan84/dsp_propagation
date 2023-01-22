function [h_xt,time1,time3,time_vec]=ht_circle_proc(X,Y,Z,freq,radius,c);
% subroutine used to compute pressure at a certain point in space (X,Y,Z)
% with respect to the centre of a piston source @ (0,0,0)
% only works as subrotine for spreading, i.e. will only compute along
% central line from source. Fpr X>0 and Y>0 use ht_circle
% written by Matt Simspon 2004
rho=sqrt((X^2)+(Y^2));
time1=Z/c;

% perpendicular distance
time2=(1/c)*sqrt((Z^2)+((rho-radius)^2));
%distance to closest edge
time3=(1/c)*sqrt((Z^2)+((rho-radius)^2));
% distance to farthest edge

sample_freq=1/((time3-time1)/60);    % empirical type fix to get roughly right sampling frequency
if sample_freq<10*freq;sample_freq=10*freq;            % to get at least fs=10*f
end
time_vec=[time1:1/sample_freq:time3];       % times needed as impulse response for rest are zero

% next bit is setting correct impulse response functions for points in and outside circle and all times needed.
if rho<=radius;
for count=1:length(time_vec);
if time_vec(count)>=time1 & time_vec(count)<time2;
h_xt(count)=1;                
elseif time_vec(count)>=time2 & time_vec(count)<time3;
h_xt(count)=(1/pi)*acos((((c*time_vec(count))^2)-(Z^2)+(rho^2)-(radius^2))/(2*rho*(sqrt(((c*time_vec(count))^2)-(Z^2)))));
if isnan(h_xt(count))==1; 
h_xt(count)=0;
end
elseif time_vec(count)>=time3;
h_xt(count)=0;
end                                                
end
end