function mip = MIP(V, view)
%%MIP Compute maximum intensity projection in given view.
%   'V' is the volumn data.
%   Three possible inputs of 'view':
%   'axial'(defaule), 'sagittal', 'coronal'.

% Set default argument
if nargin < 2
    view = 'axial';
end

% Obtain the size of input volumn data
[w, h, d] = size(V);

switch view

    case 'axial'
        % Obtain the MIP in axial view
        mip = max(V, [], 3);
    
    case 'sagittal'
        % Obtain the MIP in sagittal view
        mip = max(V, [], 2);
        mip = flip(reshape(mip, w, d)', 1);
        
    case 'coronal'
        % Obtain the MIP in coronal view
        mip = max(V, [], 1);
        mip = flip(reshape(mip, h, d)', 1);

end

end