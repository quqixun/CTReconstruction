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

sg = g2;
%sg = g;

[gl, gt] = size(sg);
hfgl = floor(gl / 2);

iw = 2 * floor(gl / (2 * sqrt(2)));
hfiw = iw / 2;

%% Setting
% Number of back projection
N = 180;

%% Back Projection
img_bp = zeros(iw);
[posX, posY] = meshgrid((1:iw) - hfiw);

igt = floor(gt / N);
for t = 1:igt:gt
   
    pos = posX * cosd(t) + posY * sind(t) + hfgl;
    img_bp = img_bp + interp1(1:gl, sg(:, t), pos);
    
end

img_bp_f = img_bp * (pi / (2 * N));

%% Plot results
figure(3)
imagesc(img_bp), colormap gray
axis('off')

figure
imagesc(img_bp_f), colormap gray
axis('off')
