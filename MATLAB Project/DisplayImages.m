function DisplayImages(figNum,images,imageNames)
% Displays a sequence of images on a single figure as a row of subplots
% Inputs: figNum, the number of the figure to open
%         image, a 1D cell array of images (each element contains a 3D 
%                array of uint8 values representing an RGB image
%         imageNames, a 1D cell array of image names (each element 
%                contains the name of the corresponding image from the
%                images cell array
% Outputs: none (but a figure is displayed when this function is called)
% Author: Peter Bier

figure(figNum);
numFrames = length(images);
for i=1:numFrames
    subplot(1,numFrames,i);
    image(images{i});
    title(imageNames{i});
    axis off
    axis equal
    axis tight
    
end


