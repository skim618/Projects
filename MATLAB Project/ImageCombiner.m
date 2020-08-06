% ImageCombiner is the main project script for the enggen131SC 2019
% Matlab project.  Once you have implemented the required functions
% that this script calls, running this script will allow you to create
% action shots or remove action from a set of provided images.
% The images may be in the form of frames extracted from a movie file
% or images that are contained within a particular directory
% Author: Peter Bier
clear
clc
clf

% Determine whether we are processing movies or images from a directory
response = ''; % for testing purpose you may like to hard code the response to m or i here
while ~(strncmp(response,'m',1) || strncmp(response,'i',1))
    disp('What would you like to process?');
    response = lower(input('Enter m for movie frames or i for images from a directory:','s'));
end


if strncmp(response,'m',1)
    % We are reading frames from a video.
    % Note that the video file will need to be in the working directory
    % or the matlab path for Matlab to find it.
    
    % Get the filename from the user, note that Matlab comes with a movie
    % file called xylophone.mp4 and hitting enter will default the name to
    % this
    filename = input('Please enter the name of the video file (or hit enter to default to xylophone.mp4:','s');
    if length(filename) == 0
        filename='xylophone.mp4';
    end
    
    % Get the frame information from the user
    firstFrame = input('Please enter the first frame you want to fetch:');
    stepSize = input('Please enter the step size you want to use:');
    numFrames = input('Please enter the number of frames to fetch:');
    
    % Generate the list of frames to extract from the move
    % You will need to implement GenerateFrameList
    frameNumbers = GenerateFrameList(firstFrame,stepSize,numFrames);
    % Extract the list of frames from the movie using the supplied
    % ReadFrames function
    images = ReadFrames(filename,frameNumbers);
    
    % create names from the frame numbers, useful when displaying the
    % images
    for i=1:length(frameNumbers)
        imageNames{i} = ['Frame ' num2str(frameNumbers(i))];
    end
else
    % We are not reading movies, so must be reading images from a directory    
    
    % Determine the directory to use and the file type to fetch
    dirname = input('Please enter the name of the directory:','s');
    fileType = input('Please enter the three letter file type, e.g. jpg or png:','s');
    
    % Create a list of images to read
    % You will need to implement GenerateImageList
    [imageNames] = GenerateImageList(dirname,fileType);
    
    % Read in the list of specified images
    % You will need to implement ReadImages
    images = ReadImages(dirname,imageNames);
    
    % For convenience you may like to comment out the calls to
    % GenerateImageList and ReadImages
    % function above and temporarily assign some values for images
    
    % ImageNames = {'pesky01.png','pesky02.png','pesky03.png'};
    % images{1} = imread('pesky\pesky01.png')
    % images{2} = imread('pesky\pesky02.png')
    % images{3} = imread('pesky\pesky03.png')
end

% Display the fetched images on figure 1 using the supplied DisplayImages
% function
DisplayImages(1,images, imageNames);

% Remove unwanted data from the sequence of images
% Display the resulting image and write it to a file
staticImage = RemoveAction(images);
figure(2)
image(staticImage);
title('Image with object removed');
imwrite(staticImage,'actionRemoved.png');

% Ceate an action shot from the sequence of images
% Display the resulting composite image and write it to a fiel
actionImage = ActionShot(images);
figure(3)
image(actionImage);

title('Action shot');
imwrite(actionImage,'actionShot.png');
