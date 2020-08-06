function [MedRed,MedGreen,MedBlue] = MedianPixel(RGB)
%The MedianPixel function calculates the median RGB values from a list of pixels.
% Input: RGB=A 1xnx3 3D array of RGB values representing a list of pixels
% Output: MedRed=The median red value, which will be the median of all the 
%                R values from the list of pixels
%         MedGreen=The median green value, which will be the median of all 
%                  the G values from the list of pixels
%         MedBlue=The median blue value, which will be the median of all 
%                 the B values from the list of pixels
% Author: Sooyong Kim

%Calculate the median for each row in the RGB array, and round it
Medians = round(median(RGB,2));

%Assign the values to its respective outputs
MedRed=Medians(1);
MedGreen=Medians(2);
MedBlue=Medians(3);

end

