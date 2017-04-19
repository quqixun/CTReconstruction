%% Implementation Image Reconstruction
%
% In this script, CT image is reconstructed by
% the method of back projection.
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
sg = g2;
%sg = g;

% Obtain the size of sinogram and compute the
% size of reconstructed image
[gl, gt] = size(sg);
hfgl = floor(gl / 2);

iw = 2 * floor(gl / (2 * sqrt(2)));
hfiw = iw / 2;

%% Setting
% Number of back projection
N = 180;

%% Back Projection
% Initialize the reconstructed image
img_bp = zeros(iw);

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
    img_bp = img_bp + interp1(1:gl, sg(:, t), pos);
    
end

% Multiply the factor
img_bp_f = img_bp * (pi / (2 * N));

%% Plot results
% Plot back projection result
figure
imagesc(img_bp), colormap gray
axis('off')

% Plot back projection result that
% multiplies the factor
figure
imagesc(img_bp_f), colormap gray
axis('off')
