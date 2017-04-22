%% Maximum Intensity Projection (MIP)
% Script and function are tested in Matlab 2016b.
% Coded by Qixun Qu
% quqixun@gmail.com
% 2017/04/22

%% Clean Environment
clc
clear
close all

%% Load Data
% Set the number of input dicom images
slices_num = 113;

% Load all images, the dimension of V is
% 512 by 512 by 113
V = load_volume('dicom_folder', slices_num);

%% Maximum Intensity Projection
% Compute MIP in three different views
mip_axial = MIP(V, 'axial');
mip_coronal = MIP(V, 'coronal');
mip_sagittal = MIP(V, 'sagittal');

% Show results
show_view(mip_axial, 'Axial View')
show_view(mip_coronal, 'Coronal View')
show_view(mip_sagittal, 'Sagittal View')

%% Rotate Axial Plane
% Rotate axial plane, obtain MIP in 360 degrees
views_360 = rotate_axial_plane(V);

%% Write Results into a Video
% Write all views into a video
write_video(views_360, 'mip_360.avi')
