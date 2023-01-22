function [out]=blockplot(x,y,ye,map);
% program which plots block image of y vs x which incorporates region
% spanned by errors in y (ye)
% map denotes coror of block e.g. [1,0,0] is red and [0,0,0] is black
% output is of no use, simply allows blockplot to operate as a function.
% need to change x and y labels  on last two lines

out=1;
mv=y;ev=ye;xv=x;
ypts=[mv+ev,fliplr(mv-ev),mv(1)+ev(1)];
xpts=[xv,fliplr(xv),xv(1)];
%map=[mapcol+0.1 mapcol mapcol-0.1];
block=fill(xpts,ypts,map);
set(block,'facealpha',0.5);
set(block,'edgecolor',map);
set(gca,'fontsize',16);
xlabel('frequency (kHz)');
ylabel('Attenuation (dB/m)');