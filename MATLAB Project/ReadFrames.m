function [frames] = ReadFrames(fileName,frameNumbers)
% Takes a a movie file name, 1D array of frames to fetch, 
% and returns a cell array of RGB images.
% Author: Peter Bier

% Create a video object for the video
vidObj = VideoReader(fileName);
frameCount = 0;

% Iterate through the frame numbers and store the specified frames.
for i=1:length(frameNumbers)
    frameToFetch = frameNumbers(i)
    while frameToFetch > frameCount    
        currentFrame = readFrame(vidObj);
        frameCount = frameCount + 1;
    end
    frames{i}= currentFrame;
end
    
end

