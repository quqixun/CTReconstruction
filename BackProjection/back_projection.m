function img_bp = back_projection(sg, N, angle)
%%BACK_PROJECTION Implement the back projection approach to reconstruct
%image from given sinogram.
%   Input argument:
%   - sg : the known sinogram
%   - N : the number of projections used in construction,
%         default is the number of columns of sinogram
%   - angle : back projection of one data in sinogram at given angle
%   Output:
%   - img_bp : the reconstruction result

% Set the default value for the number of
% projection that applied to reconstruct
if nargin < 2 || isempty(N)
    N = size(sg, 2);
end

% Set the default value of angle
if nargin < 3 || isempty(angle)
    angle = -1;
else
    N = length(angle);
end

% Obtain the size of sinogram and compute the
% size of reconstructed image
[gl, gt] = size(sg);
hfgl = floor(gl / 2);

iw = 2 * floor(gl / (2 * sqrt(2)));
hfiw = iw / 2;

% Back Projection
% Initialize the reconstructed image
img_bp = zeros(iw);

% Compute some arguments for back projection
% Positions map of reconstructed image
[posX, posY] = meshgrid((1:iw) - hfiw);

if angle == -1
    % If back projection is generated from
    % several data in sinogram
    % Calculate the degree interval
    igt = floor(gt / N);
    angles_array = 1:igt:gt;
else
    % If back projection is created only by
    % one data in sinogram at given angle
    angles_array = angle;
end

% Run N times back projection 
for t = angles_array

    % Calculate the position in sinogram
    pos = posX * cosd(t) + posY * sind(t) + hfgl;
    % Accumulate projection of each degree sinogram
    img_bp = img_bp + interp1(1:gl, sg(:, t), pos);
    
end

% Multiply the factor
img_bp = img_bp * (pi / (2 * N));

% Plot results
% Plot back projection result that
% multiplies the factor
figure
imagesc(img_bp), colormap gray
axis('off')

end