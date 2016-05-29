clear all;
close all;

x_len = 210;
z_len = 490;

scen = zeros(x_len, z_len);

bottom_w = 200;
top_w = 30;

start = 100;

for i = 1:(bottom_w-top_w)/2
    for j = 1:bottom_w-2*i
        scen(i, start+i+j) = 1;
        scen(x_len-i+1, start+i+j) = 1;
    end
end


scen = logical(scen);

colormap([0,0,0; 1,1,1])
image(scen); %TODO: make better color-scheme
axis equal;
xlabel('z');
ylabel('x');

%imwrite(scen, 'test.jpg');
pgmwrite(scen, 'Bottleneck');