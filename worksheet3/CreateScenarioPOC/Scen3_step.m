%depth is y-axes, that means view is on x-z plane (looking from the side)
% the boundaries on top and 
%y-z plane where x = {0,xmax} is not given and set in the simulation
clear all;
close all;

x_len = 30;
z_len = 50;

scen = zeros(x_len, z_len); %ghost layers are not required

blocksize = (x_len+2)/2;
assert(mod(blocksize, 2) == 0);

midp = [x_len-blocksize/2+2, blocksize/2];
scen = drawRect_center(scen, midp, blocksize,blocksize-2,1);
scen = logical(scen);

colormap([0,0,0; 1,1,1])
image(scen);
axis equal;
xlabel('z');
ylabel('x');

pgmwrite(scen, 'Step');