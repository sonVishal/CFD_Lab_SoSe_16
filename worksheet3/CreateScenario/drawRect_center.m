function [ A ] = drawRect_center(A, center, width, height, val)

assert(mod(width,2) ==0 && mod(height, 2) ==0);

%right upper corner
ruc = [center(1)-height/2, center(2)-width/2];
A(ruc(1):ruc(1)+height, ruc(2):ruc(2)+width) = val;



end

