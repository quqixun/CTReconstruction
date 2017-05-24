%% Maximum Intensity Projection (MIP)
% Lab 1 for Diagnostic Imaging
% Group Members: Qixun Qu, Yankun Xu, Zihui Wang
% Script and function are tested in Matlab 2016b.
% 2017/04/22

%% Clean Environment
clc
clear
close all

%% Load Data
% Set the number of input dicom images
% Here, all images are read into the memory
slice_to_load = 1 : 113;

% Load all images, the dimension of V is
% 512 by 512 by 113
V = load_DICOM_volume('spiral_CT_mandible', slice_to_load);

%% Maximum Intensity Projection
% Compute MIP in three different views
mip_axial = MIP(V, 'axial');
mip_coronal = MIP(V, 'coronal');
mip_sagittal = MIP(V, 'sagittal');

% Show results
show_view(mip_axial, 'Axial View')
show_view(mip_coronal, 'Coronal View')
show_view(mip_sagittal, 'Sagittal View')

%% Rotate Volume
% Rotate axial plane, obtain MIP from many angles
% The viewing angle is fixed, rotate the volume
views_360 = rotate_volume(V);

%% Rotate Views
% Rotate axial plane, obtain MIP from many angles
% The volume is fixed, change the perspective
% This method is too slow
% views_360 = rotate_views(V);

%% Write Results into a Video
% Write all views into a video
write_video(views_360, 'mip_360.avi', 10)
