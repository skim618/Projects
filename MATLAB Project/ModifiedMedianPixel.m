function [Medians] = ModifiedMedianPixel(RGB)
%The ModifiedMedianPixel function calculates the median RGB values from a list of pixels.
% Input: RGB=A nxnx3 3D array of RGB values representing a list of pixels
% Output: Medians = An array of the medians in each row of the RGB values
% Author: Sooyong Kim

%Calculate the median for each row in the RGB array, and round it
Medians = round(median(RGB,2));
end

