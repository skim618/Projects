function [DistanceSQArray] = ModifiedPixelDistance(ThreeDPoint1,ThreeDPoint2)
%The PixelDistance function calculates the square of the distance between
%two points in colour space.
% Inputs: ThreeDPoint1=An array containing three elements representing a 
%                      point in 3D colour space
%         ThreeDPoint2=An array containing three elements representing a 
%                      second point in 3D colour space
% Outputs: DistanceSQArray=An array of the squares of the distance between 
%                          the two points in 3D colour space.
% Author: Sooyong Kim

%If values are uint8 integers, double them.
ThreeDPoint1=double(ThreeDPoint1);
ThreeDPoint2=double(ThreeDPoint2);

%Element by Element operations:
DistanceSQArray=sum(((ThreeDPoint1-ThreeDPoint2).^2),3);

end

