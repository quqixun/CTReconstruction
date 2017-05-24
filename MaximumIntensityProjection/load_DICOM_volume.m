function [v,t] = load_DICOM_volume(directory,slices_to_load)
% LOAD_DICOM_VOLUME Loads in the specified DICOM files from the specified
% directory into a 3D array.
%    V = LOAD_DICOM_VOLUME(D,S) loads in the DICOM files (files with the
%    extension .dcm) from DIRECTORY into the 3D array V with dimensions:
%    row, column, slice. The parameter S is a vector containing a range
%    of slice numbers; e.g. 1:5
 
% keep a record of the current directory
current_directory = cd;
 
% go to the directory containing the DICOM files
cd(directory);
 
% get a list of all the files with the extension .dcm
files = dir('*.dcm');
 
h = waitbar(0,'Loading volume...');
info = dicominfo(files(1).name);
t = info.ContentTime;
 
for i=slices_to_load
    v(:,:,i-slices_to_load(1)+1) = dicomread(files(i).name);
    waitbar(i/length(slices_to_load));
end
close(h)

% return to the original directory
cd(current_directory);
