
% SOURCE: 
% https://youneedtoprogram.wordpress.com/2015/03/12/pgmwrite-m-a-matlab-function-to-write-a-matrix-to-a-pgm-file/
% pgmwrite.m - MATLAB function to write a matrix
% of pixel values to a greyscale PGM image file.
%
% Written by Ted Burke - last updated 12-3-2015
%

function pgmwrite(pixels, name)
    % Check matrix dimensions which will determine
    % the width and height of the image
    s = size(pixels);
    height = s(1);
    width = s(2);
    
    % Open a file for writing
    fid = fopen([name, '.pgm'],'w');
    
    % Write the PGM image header to the file
    fprintf(fid,'P2\n');
    fprintf(fid,'%d %d\n', width, height);
    fprintf(fid,'2\n');
    
    % Write the pixel values from the matrix into
    % the file
    y = 1;
    while y <= height
        x = 1;
        while x <= width
            fprintf(fid, '%i ', pixels(y,x));
            x = x + 1;
        end
        fprintf(fid, '\n');
        y = y + 1;
    end
    
    % Close the file
    fclose(fid)
end
