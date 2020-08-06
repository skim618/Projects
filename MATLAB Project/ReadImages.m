function [pictures] = ReadImages(ImageDirectory,FileArray)
%The ReadImages function reads a specified list of images given the 
%filenames and the directory the files are located in.
% Input: ImageDirectory = A string containing the name of the directory the
%                         images are contained in
%        FileArray = A 1xn 1D cell array containing n strings where each 
%                    element is a filename of an image to read
% Output: pictures = 1xn 1D cell array containing n images, where each 
%                    element is an RGB image.
% Author: Sooyong Kim

%Locate and change the directory, while saving the original working space
oldFolder = cd(ImageDirectory);
cd;

%Allocate an array to store the images, where the image is an RGB image
pictures={};

%Create a for loop to process the images in the FileArray, and store them
%in the created pictures array
for i=1:length(FileArray)
    pictures{i}=imread(FileArray{i});
end

%Return to the original working space
cd (oldFolder);
cd;
end


