%% Implementation Image Reconstruction
%
% In this script, CT image is reconstructed by
% the method of filtered back projection.
%
% Created by Qixun QU
% 2017/04/19

%% Clearn Environment
clc
clear
close all

%% Load Data
load data.mat
load data2.mat

% In this case, g and g2 are sinogram data,
% they are n by 180 matrix, which means in
% g(l, theta), the range of theta is from 1
% to 180, in each degree, l has n values
%sg = g2;
sg = g;

[gl, gt] = size(sg);
hfgl = floor(gl / 2);

iw = 2 * floor(gl / (2 * sqrt(2)));
hfiw = iw / 2;

%% Settings
% method:
% 0 for ramlak, 1 for hamming
method = 1;

% The number of back projection
N = 180;

%% Filtered Back Projection
% Compute the Ramlak filter
gx = [0:hfgl, hfgl - 1:-1:1];
if mod(gl, 2) ~= 0
    gx = [gx, 0];
end

ramlak = 2 * gx / gl;

switch method

    case 0 % Use Ramlak filter
        H = ramlak;
    case 1 % Use Hamming filter
        hamming = 0.54 - 0.46 * cos(2 * pi * (0:gl-1) / gl);
        H = [hamming(hfgl:gl), hamming(1:hfgl-1)] .* ramlak;
    otherwise
        fprintf('Wrong method, it should be 0 or 1.')

end

% Compute Fourier transformation of sinogram
gf = fft(sg, [], 1);
% Multiply the filter in frequency domain
gff = bsxfun(@times, gf, H');
% Do inverse Fourier transformation
gffi = real(ifft(gff, [], 1));

% Initialize the reconstructed image
img_fbp = zeros(iw);

% Compute some arguments for back projection
% Positions map of reconstructed image
[posX, posY] = meshgrid((1:iw) - hfiw);
% The degree interval
igt = floor(gt / N);

% Run N times back projection
for t = 1:igt:gt

    % Calculate the position in sinogram
    pos = posX * cosd(t) + posY * sind(t) + hfgl;
    % Accumulate projection of each degree sinogram
    img_fbp = img_fbp + interp1(1:gl, gffi(:, t), pos);
    
end

% Multiply the factor
img_fbp_f = img_fbp * (pi / (2 * N));

%% Plot results
% Plot back projection result
figure
imagesc(img_fbp), colormap gray
axis('off')

% Plot back projection result that
% multiplies the factor
figure
imagesc(img_fbp_f), colormap gray
axis('off')