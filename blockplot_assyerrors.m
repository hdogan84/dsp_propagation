function [out]=blockplot_assyerrors(x,y,ye_up,ye_down,map);
% program which plots block image of y vs x which incorporates region
% spanned by errors in y (ye_up and ye_down)
% note that ye_down must all be negative values and ye_up positive values
% as block plot except allows for assymetirc errors 
% map denotes coror of block e.g. [1,0,0] is red and [0,0,0] is black
% output is of no use, simply allows blockplot_assyerrors to operate as a function.
% need to change x and y labels  on last two lines

out=1;
mv=y;ev1=ye_up;ev2=ye_down;xv=x;
ypts=[mv+ev1,fliplr(mv+ev2),mv(1)+ev1(1)];
xpts=[xv,fliplr(xv),xv(1)];
block=fill(xpts,ypts,map);
hold on;
set(block,'facealpha',0.5);
set(block,'edgecolor',map);
set(gca,'fontsize',16);
xlabel('frequency (kHz)');
ylabel('attenuation (dB/m)');