function views = rotate_views( V, degrees )
%%ROTATE_VIEWS Rotate the viewing angles on axial plane,
% save all views in given degrees.
%   'V' is the volume data.
%   'degrees' is a vector consists of degrees, default is
%   from 0 to 360.

% Set default value for degrees
if nargin < 2 || isempty(degrees)
    degrees = 0 : 360;
end

% Convert value into double
V = double(V);

% Obtain the size of input volumn data
[w, ~, d] = size(V);

% Obtain the length of 'degrees'
d_num = length(degrees);

% Initialize the output
views = zeros(d, w, d_num);

% Compute the number of rays
p = round(sqrt(2) * w);

% Compute the index of start and stop position
% after the slice being projected
width_diff = round((p - w)/2);
w_start = width_diff + 1;
w_stop = width_diff + w;

% Initialize temporary variables
slice_temp = zeros(d, w);
max_temp = zeros(1, p);

% Conut how many views have been obtained
step = 0;
h = waitbar(0, 'Rotating Axial Plane ...');

for i = degrees % In each viewing angle
    fprintf('Viewing Angle: %d\n', i)
    
    % Compute coefficient matrix, which consists of
    % positions where all rays cross the slice
    A = full(paralleltomo(w, i));
    
    for j = 1 : d % In each slice
        fprintf('\tSlice: %d\n', j)
        
        % Extract a slice from the volume
        Vj = V(:, :, j);
        
        for k = 1 : p % Obtain the maximum in each ray 
            Vj_temp = Vj(A(k, :) > 0);
            if isempty(Vj_temp)
                max_temp(k) = 0;
            else
                max_temp(k) = max(Vj_temp(:));
            end
        end
        
        % slice_temp is the MIP of the volume at one viewing angle
        slice_temp(j, :) = flip(max_temp(w_start : w_stop), 2);
    end
    
    fprintf('\n')
    
    % Save MIP into output
    step = step + 1;
    views(:, :, step) = flip(slice_temp);
    
    waitbar(step / d_num)
    
end
close(h)

end