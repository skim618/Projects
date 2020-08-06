function [ActionImage] = ActionShot(images)
%The ActionShot function creates an action shot image by finding the pixels 
%from a stack of images that are most distant from the median RGB values
% Input: image = A 1xn 1D cell array containing n images, where each element 
%                is an RGB image
% Output: ActionImage = An action shot in the form of an RGB image
% Author: Sooyong Kim

%Obtain the dimensions of the image using one of the images as a sample
[rows,cols,col]=size(images{1});

%Allocate arrays for the action image (rows,cols,col) in uint8 integers; and 
%for the RGB pixels to be used to obtain the most distant pixel.
RGB = zeros(1,length(images),3);
ActionImage = zeros(rows,cols,col, 'uint8');

%This for loop processes each pixel of the images in the same position and 
%stores them in a 1xnx3 RGB array, then calls the MostDistantPixel function 
%to calculate the most distant pixel, then stores the outputs into the ActionImage array.
for i=1:rows
    for j=1:cols
        for k=1:col
            for l=1:length(images)
                RGB(1,l,k)=images{l}(i,j,k);
            end
        end
        [DistRed,DistGreen,DistBlue] = MostDistantPixel(RGB);
        ActionImage(i,j,1) = DistRed;
        ActionImage(i,j,2) = DistGreen;
        ActionImage(i,j,3) = DistBlue;
    end
end

end



