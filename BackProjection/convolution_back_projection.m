function img_cbp = convolution_back_projection(sg, N, filter, angle)
%%CONVOLUTION_BACK_PROJECTION Implement the convolution back projection
%approach to reconstruct image from given sinogram.
%   Input argument:
%   - sg : the known sinogram
%   - filter : the name of the filter used to filter the projection,
%     default is 'hamming'
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

% Set the default value of filter
if nargin < 3 || isempty(filter)
    filter = 'hamming';
end

% Set the default value of angle
if nargin < 4 || isempty(angle)
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

% Convolution Back Projection
% Compute the Ram-Lak filter
gx = [0:hfgl, hfgl - 1:-1:1];
if mod(gl, 2) ~= 0
    gx = [gx, 0];
end

ramlak = 2 * gx / gl;

switch filter

    case 'ramlak' % Use Ramlak filter
        H = ramlak;
    case 'hamming' % Use Hamming filter
        hamming = 0.54 - 0.46 * cos(2 * pi * (0:gl-1) / gl);
        H = [hamming(hfgl:gl), hamming(1:hfgl-1)] .* ramlak;
    otherwise
        fprintf('Wrong filter, it should be ''ramlak'' or ''hamming''.')

end

% Compute inverse Fourier transformation of filter
c = real(ifftshift(ifft(H)));
% Convolve the sinogram with the IFT of filter
sgc = convn(sg, c', 'same');

% Initialize the reconstructed image
img_cbp = zeros(iw);

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
    img_cbp = img_cbp + interp1(1:gl, sgc(:, t), pos);
    
end

% Multiply the factor
img_cbp = img_cbp * (pi / (2 * N));

% Plot results
% Plot back projection result
figure
imagesc(img_cbp), colormap gray
axis('off')

end