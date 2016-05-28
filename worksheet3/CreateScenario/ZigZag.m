%depth is y-axes, that means view is on x-z plane (looking from the side)
% the boundaries on top and 
%y-z plane where x = {0,xmax} is not given and set in the simulation
clear all;
close all;

x_len = 210;
z_len = 490;

scen = zeros(x_len, z_len); %+2 is ghoast layers

p1 = [x_len - 70, 120];
p2 = [71, 250];
p3 = [x_len - 70, 380];


scen = drawRect_center(scen, p1, 40, 140, 1);
scen = drawRect_center(scen, p2, 40, 140, 1);
scen = drawRect_center(scen, p3, 40, 140, 1);
%scen = drawRect_center(scen, right_upper, blocksize,blocksize,1);

scen = logical(scen);

colormap([0,0,0; 1,1,1])
image(scen); %TODO: make better color-scheme
axis equal;
xlabel('z');
ylabel('x');

%imwrite(scen, 'test.jpg');
pgmwrite(scen, 'ZigZag');