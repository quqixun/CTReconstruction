function img = image_normalize( img )
%%IMAGE_NORMALIZE Convert the range of image's values
% from 0 to 1.

% Convert image from int16 to double
img = double(img);

% Normalize the image
imin = min(img(:));
imax = max(img(:));
img = (img - imin) / (imax - imin);

end