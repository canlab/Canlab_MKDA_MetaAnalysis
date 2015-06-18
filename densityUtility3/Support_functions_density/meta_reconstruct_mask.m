% reconstruct mask for a study, given indicator
% v2 = meta_reconstruct_mask(indic, xyzlist, V.dim(1:3), [use values], [V], [imagename]);
%
% This function returns 3D mask values and optionally writes an image file
% if V and imagename are entered as additional arguments
%
% dims are mask dimensions; V is spm_vol structure
% indic: see meta_fast_sphere_conv
% xyzlist: see meta_read_mask
%
% [use values]: optional; use indic values in voxels (for saving count
% values)

function v2 = meta_reconstruct_mask(indic, xyzlist, dims, varargin)

    if length(varargin) > 0, valueflag = varargin{1}; else valueflag = 0; end

    wh = find(indic);   % which indices, in in-mask XYZ list

    % reconstruct mask
    v2 = zeros(dims);
    wh2 = sub2ind(dims, xyzlist(wh,1), xyzlist(wh,2), xyzlist(wh,3));

    if valueflag
        v2(wh2) = indic(wh);
    else
        v2(wh2) = 1;            % faster
    end

    % write image file only if additional arguments are entered
    if length(varargin) > 1
        V = varargin{2};
    else
        return
    end

    if length(varargin) > 2
        V.fname = varargin{3};
    end

    warning off Matlab:DivideByZero   % turn off in case of empty image
    warning off MATLAB:intConvertNaN
    str1 = sprintf('Writing: %s', V.fname);
    fprintf(1, str1);
    spm_write_vol(V, v2);

    erase_string(str1);
end



function erase_string(str1)
    fprintf(1, repmat('\b', 1, length(str1))); % erase string
    %fprintf(1,'\n');
end
