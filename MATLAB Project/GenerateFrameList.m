function [frames] = GenerateFrameList(F0,SS,n)
%GenerateFrameList generates a list of frames we are interested in which 
%can be used by other functions. In particular it will be useful for 
%determining which frames to extract from a movie file.
% Inputs: F0 = Frames-naut, which is the starting frame number
%         SS = Step size
%         n = The number of frames to generate
% Outputs: frames = a 1 x 'n' 1D array, where n is the desired number of frames (n).
%                   The first element of the array will be the starting frame 
%                   number and each subsequent element will have a frame 
%                   value that is the step size greater than the last.
% Author: Sooyong Kim

%The first value will always be F0
frames(1)=F0;

%Make a for loop to add the subsequent elements into the array
for i=1:(n-1)
    frames(i+1)=frames(i)+SS;
end
end


