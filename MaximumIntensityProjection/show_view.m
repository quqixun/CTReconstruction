function show_view( img, str )
%%SHOW_VIEW Show image with given title.

if nargin < 2 || isempty(str)
    str = '';
end

figure
imshow(img, [])
title(str)

end
