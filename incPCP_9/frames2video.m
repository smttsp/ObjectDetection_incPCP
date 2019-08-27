function [ outputVideo] = frames2video( im_folder)


outputVideo = VideoWriter(fullfile(im_folder,'shuttle_out.avi'));
outputVideo.FrameRate = 30;
open(outputVideo)
% Loop through the image sequence, load each image, and then write it to the video.

imageNames = dir([im_folder '*jpg']);
for ii = 1:length(imageNames)
   img = imread(fullfile(im_folder,imageNames(ii).name));
   writeVideo(outputVideo,img)
end
% Finalize the video file.

close(outputVideo)

end

