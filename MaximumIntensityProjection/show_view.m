function show_view( img, str )
%%SHOW_VIEW Show image with given title.

figure
imshow(img, [0, max(img(:))])
title(str)

end
