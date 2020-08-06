function [RemovedImage] = RemoveAction(images)
%The RemoveAction function creates an image that has the action removed by
%applying a median filter
% Inputs: images = A 1xn 1D cell array containing n images, where each
%                 element is an RGB image
% Outputs: RemovedImage = An RGB image that was obtained by taking the median
%                         RGB values of the stack of corresponding pixels from
%                         the source images.
% Author: Sooyong Kim

%We must process each pixel individually and find the median

%Obtain the dimensions of the image using one of the images as a sample
[rows,cols,col]=size(images{1});

%Allocate arrays for the action image (rows,cols,col) in uint8 integers;
%and for the RGB pixels to be used to obtain the median pixel.
RGB = zeros(rows,length(images),col);
RemovedImage = zeros(rows,cols,col, 'uint8');

%This for loop processes each pixel of the images in the same position and
%stores them in a nxnx3 RGB array, then calls the ModifiedMedianPixel function to
%calculate the median of the pixels, then stores the outputs into the
%RemovedImage array.
for i=1:cols
    for j=1:rows
        for k=1:col
            for l=1:length(images)
                RGB(j,l,k)=images{l}(j,i,k);
            end
        end
    end
    [Medians] = ModifiedMedianPixel(RGB);
    RemovedImage((1:rows),i,(1:col)) = Medians;
end

end

