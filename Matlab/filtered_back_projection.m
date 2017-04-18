%% Created by Qixun Qu
% quqixun@gmail.com
% 2017/04/18

%% Implementation Image Reconstruction

%%
clc
clear
close all

%% Load Data
load data.mat
load data2.mat

%sg = g2;
sg = g;

[gl, gt] = size(sg);
hfgl = floor(gl / 2);

iw = 2 * floor(gl / (2 * sqrt(2)));
hfiw = iw / 2;

%% Setting
% method:
% 0 for ramlak, 1 for hamming
method = 1;

% The number of back projection
N = 180;

%% Filtered Back Projection
if mod(gl, 2) == 0
    gx = [0:hfgl, hfgl - 1:-1:1];
else
    gx = [0:hfgl, hfgl - 1:-1:0];
end

ramlak = 2 * gx / gl;

switch method

    case 0
        H = ramlak;
    case 1
        hamming = 0.54 - 0.46 * cos(2 * pi * (0:gl-1) / gl);
        H = [hamming(hfgl:gl), hamming(1:hfgl-1)] .* ramlak;
    otherwise
        fprintf('Wrong method, it should be 0 or 1.')

end

gf = fft(sg, [], 1);
gff = bsxfun(@times, gf, H');
gffi = real(ifft(gff, [], 1));

img_fbp = zeros(iw);
[posX, posY] = meshgrid((1:iw) - hfiw);

igt = floor(gt / N);
for t = 1:igt:gt
   
    pos = posX * cosd(t) + posY * sind(t) + hfgl;
    img_fbp = img_fbp + interp1(1:gl, gffi(:, t), pos);
    
end

img_fbp_f = img_fbp * (pi / (2 * N));

%% Plot results
figure
imagesc(img_fbp), colormap gray
axis('off')

figure
imagesc(img_fbp_f), colormap gray
axis('off')