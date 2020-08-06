function [MostDistantPoints] = ModifedMostDistantPixel(RGB)
%This modified MostDistantPixel function calculates the pixel from a list that is 
%most distant from the median RGB values of that list. The distance metric 
%to be used is that described in the PixelDistance function
% Inputs: RGB = A nxnx3 3D array of RGB values representing a list of pixels
% Outputs: DistantPoints = Array of the distant points
% Author: Sooyong Kim

%Call the median of each colour using the ModifiedMedianPixel function
[Medians] = ModifiedMedianPixel(RGB);

[rows,cols,col] = size (RGB);

%Allocate an array to store the distances between the median and all points
Distances=zeros(rows,cols);

%This calculates the distance between each point and the median, and stores
%them into the Distanes array.
for j=1:cols
    [DistanceSQArray] = ModifiedPixelDistance(Medians,RGB((1:rows),j,(1:col)));
    Distances((1:rows),j) = DistanceSQArray;
end

MaxDistances=max(Distances,[],2);

%Locate  the position of the maximum distance in the array, where 'Position' 
%is the position of the MaxValue in the array
MostDistantPoints = zeros(rows,1,3);
[row,col]=find(Distances==MaxDistances);
for i=1:length(col)
    MostDistantPoints(i,1,(1:3)) = RGB(row(i),col(i),:) ;
end


%The most distant point is the 'k' element in the array









