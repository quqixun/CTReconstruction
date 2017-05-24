function write_video( views, name, sec )
%%WRITE_VIDEO Write all view images into a video.
%   The 3rd dimension of 'views' indicates the number
%   of views.
%   'name' is the file name of the video.
%   'sec' is the the length of video in seconds.

% Set default length for video file
if nargin < 3 || isempty(sec)
    sec = 10;
end

% Set default name for video file
if nargin < 2 || isempty(name)
    name = 'mip.avi';
end

% Create a video
mip_video = VideoWriter(name);

% Set the frame rate according to the length
mip_video.FrameRate = round(size(views, 3) / sec);

% Write image into video
open(mip_video);
for i = 1 : size(views, 3)
    
    % Normalize the view
    f = image_normalize(views(:, :, i));
    writeVideo(mip_video, f);
    
    % Method to write gif file
%     [A,map] = gray2ind(f,256);
%     delay = 1 / mip_video.FrameRate;
%     if i == 1
%         imwrite(A, map, 'mip_360.gif', 'gif',...
%                 'LoopCount', Inf, 'DelayTime', delay);
%     elseif mod(i, 10) == 0
%         imwrite(A, map, 'mip_360.gif', 'gif',...
%                 'WriteMode', 'append', 'DelayTime', delay);
%     end

end

% Close the video file
close(mip_video);

end