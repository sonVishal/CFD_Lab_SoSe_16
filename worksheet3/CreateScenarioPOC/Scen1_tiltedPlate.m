%depth is y-axes, that means view is on x-z plane (looking from the side)
% the boundaries on top and 
%y-z plane where x = {0,xmax} is not given and set in the simulation
clear all;
close all;

x_len = 212;
z_len = 492;

scen = ones(x_len+2, z_len+2); %+2 is ghoast layers
scen(2:end-1, 2:end-1) = 0; %on all four edges is now boundary

blocksize = 30;
assert(mod(blocksize, 2) == 0);

midp = [1/2*x_len, 1/3*z_len];
left_lower = [midp(1)+blocksize/2, midp(2)-blocksize/2];
right_upper = [midp(1)-blocksize/2, midp(2)+blocksize/2];

scen = drawRect_center(scen, midp, blocksize,blocksize,1);
scen = drawRect_center(scen, left_lower, blocksize,blocksize,1);
scen = drawRect_center(scen, right_upper, blocksize,blocksize,1);

scen = logical(scen);

colormap([0,0,0; 1,1,1])
image(scen); %TODO: make better color-scheme
axis equal;
xlabel('z');
ylabel('x');

%imwrite(scen, 'test.jpg');
pgmwrite(scen, 'test');