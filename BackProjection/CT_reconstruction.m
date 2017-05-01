%% Implementation of Image Reconstruction
%
% In this script, CT image is reconstructed by
% the method of back projection, filtered back
% projection and convolution back projection.
% 
% All results are saved in folder 'results'.
%
% Created by Qixun QU
% 2017/05/01

%% Clearn Environment
clc
clear
close all

%% Load Data
load data.mat
load data2.mat

%% Part A. Back Projection
% A2.The back projection of the sinogram data at 30 degree
img_bp_30 = back_projection(g, 1, 30);

% A3.Back projection summation images with different N
N = 180; % Change N
img_bp_N = back_projection(g, N);

%% Part B. Filtered Back Projection
% B2.Filtered back projection with different N
N = 180; % Change N
img_fbp_ramlak = filtered_back_projection(g, N, 'ramlak');

% B3.Apply Hamming window in filtering process
N = 180; % Change N
img_fbp_hamming = filtered_back_projection(g, N, 'hamming');

%% Part C. Convolution Back Projection
% C2.Implement convolution back projection
% with Hamming window with different N
N = 180; % Chang N
img_cbp_hamming = convolution_back_projection(g, N, 'hamming');

%% Part E. Real CT Image
% E1.Reconstruct real CT image by applying three implementations
% Back projection
ri_bp = back_projection(g2);
% Filtered back projection
ri_fbp = filtered_back_projection(g2);
% Convolution back projection
ri_cbp = convolution_back_projection(g2);
