function [ImageFiles] = GenerateImageList(ImageDirectory,FileExtension)
%The GenerateImageList function retrieves a list of the names of all the 
%images with a particular file extension that arecontained in a specified 
%directory
% Inputs: ImageDirectory = A string containing the name of the directory the
%                         images are contained in
%         FileExtension = A string containing the file extension of the 
%                        images to fetch
% Output: ImageFiles = A 1xn 1D cell array containing n strings where each 
%                     element is the filename of an image from the specified 
%                     directory that has a particular file extension
% Author: Sooyong Kim

%Obtain the structure array using the dir function
FileList = dir (ImageDirectory);

%Extract the file names from the FileList structure array
FileNames = {FileList.name};

%Now we want to compare each element in the cell whether they end in the
%file extension we want. If they do match, then this loop adds that element
%into the allocated ImageFiles cell array.

%Create a counting index to store the matching file names into our 
%ImageFile cell array
i=1;
for j=1:length(FileNames)
    if endsWith(FileNames{j}, FileExtension)==1
        ImageFiles{i} = FileNames{j};
        i=i+1;
    end
end

end

