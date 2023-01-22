clear all
p0=101000;
ivme=9.81;
cs=1480;
h=0.5;
beta=0.7;
phi=1.3;
sigma=0.072;
mu=1;
rho=1480;
pg0=p0+rho*ivme*h;
G=2.6e6;
R0=330e-6;
f=26e3;
omega=2*pi*f;
m=rho*R0*R0;
first=3*pg0*phi
second=-2*sigma*beta/R0
third=4*G*(1-beta)
acs1=omega*omega*R0*R0*rho;
acs2=1+(omega*R0/cs)^2;
fourth=acs1/acs2
XX=first+second+third+fourth;  %% Dogan
XX=first+second+third;   %% Mantouka

omega0=sqrt(XX/m);
omega0/2/pi/1e3


%RR=(3*pg0*phi+4*G*(1-beta))/rho/omega/omega;
%sqrt(RR)