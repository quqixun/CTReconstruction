function V = load_volume( directory, slices_num )
% LOAD_VOLUME Load volume data from dicom images in given folder.
 
% get a list of all the files with the extension .dcm
files = dir([directory '/*.dcm']);
 
h = waitbar(0, 'Loading ...');
for i = 1 : slices_num
    % Read dicom image and save it into output matrix
    V(:, :, i) = dicomread([directory '/' files(i).name]);
    waitbar(i / slices_num);
    
end
close(h)

end
