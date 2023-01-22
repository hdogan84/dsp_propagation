function [h_xt, time1, time3, time_vec]=ht_circle(x,y,z,freq,radius,c)
% subroutine called by circle_new_mod_1m and circle_new_mod_asfunction
% to compute pressure at a certain point in space (X,Y,Z)
% with respect to the centre of a piston source @ (0,0,0)

rho=sqrt((x^2)+(y^2));
time1=z/c;% perpendicular distance
time2=(1/c)*sqrt((z^2)+((rho-radius)^2)); %distance to closest edge
time3=(1/c)*sqrt((z^2)+((rho+radius)^2));  % distance to farthest edge

sample_freq=1/((time3-time1)/60);    % empirical type fix to get roughly right sampling frequency
if sample_freq<10*freq;sample_freq=10*freq;            % to get at least fs=10*f
end
sample_freq;
time_vec=[time1:1/sample_freq:time3];       % times needed as impulse response for rest are zero

% next bit is setting correct impulse response functions for points in and outside circle and all times needed.
if rho<=radius;
for count=1:length(time_vec),
if time_vec(count)>=time1 & time_vec(count)<time2,
h_xt(count)=1;                
elseif time_vec(count)>=time2 & time_vec(count)<time3,
h_xt(count)=(1/pi)*acos((((c*time_vec(count))^2)-(z^2)+(rho^2)-(radius^2))/(2*rho*(sqrt(((c*time_vec(count))^2)-(z^2)))));
if isnan(h_xt(count))==1, 
h_xt(count)=0;
end
elseif time_vec(count)>=time3,
h_xt(count)=0;
end                                                
end
end

if rho>radius,
for count=1:length(time_vec),
if time_vec(count)>=time1 & time_vec(count)<=time2,
h_xt(count)=0;                
elseif time_vec(count)>=time2 & time_vec(count)<=time3,
h_xt(count)=(1/pi)*acos((((c*time_vec(count))^2)-(z^2)+(rho^2)-(radius^2))/(2*rho*(sqrt(((c*time_vec(count))^2)-(z^2)))));
elseif time_vec(count)>=time3,
h_xt(count)=0;
end                                                
end
end