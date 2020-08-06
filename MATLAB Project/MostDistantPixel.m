function [DistRed,DistGreen,DistBlue] = MostDistantPixel(RGB)
%The MostDistantPixel function calculates the pixel from a list that is 
%most distant from the median RGB values of that list. The distance metric 
%to be used is that described in the PixelDistance function
% Inputs: RGB = A 1xnx3 3D array of RGB values representing a list of pixels
% Outputs: DistRed = The red value of the most distant pixel
%          DistGreen = The green value of the most distant pixel
%          Dist Blue = The blue value of the most distant pixel
% Author: Sooyong Kim

%Call the median of each colour using the MedianPixel function, and store
%them in an array
[MedRed,MedGreen,MedBlue] = MedianPixel(RGB);
Medians = [MedRed,MedGreen,MedBlue];

[rows,cols,col] = size (RGB);

%Obtain the RGB points of the given RGB array using a for loop and 
%store them in a cell array

for i=1:cols
    for j=1:3
        Point(j)=RGB(1,i,j);
    end
    Points{i}=Point;
    
end

%This calculates the distance between each point and the median, and stores
%them into the Distanes array.
for i=1:cols
    [DistanceSQ] = PixelDistance(Medians,Points{i});
    Distances(i) = DistanceSQ;
end

MaxDistance=max(Distances);

%Locate  the position of the maximum distance in the array, where 'k' is 
%the position of the MaxValue in the array
k=find(Distances==MaxDistance);

%The most distant point is the 'k' element in the array
MostDistantPoint = Points{k};

%Assign them into their respective outputs
DistRed = MostDistantPoint(1);
DistGreen = MostDistantPoint(2);
DistBlue = MostDistantPoint(3);



