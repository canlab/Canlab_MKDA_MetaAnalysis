function conmask = xyz2density(XYZmm,mask,V,str,radius,studyweight,varargin)
% conmask = xyz2density(XYZmm,mask,V,str,radius,studyweight,[enter vox xyz flag] AND [no write image])
% 
% Take a list of xyz mm coordinates and turn it into a density mask
% with spherical convolution.
% 
% XYZmm: n x 3, mask = zeros of correct dims, V = mapped volume (SPM),
% str = name of output image
% radius: sphere radius in voxels
%
% last args: enter anything are to enter coords in voxels and suppress write output
% for p0_monte_carlo_FWE, to construct a fast density map of random
% coordinates

XYZmm(any(isnan(XYZmm),2),:) = [];

if isempty(XYZmm)
    conmask = mask;
    
    
else
    if length(varargin) > 0
        XYZvox = XYZmm;
    else
        XYZvox = mm2vox(XYZmm,V.mat);   % get voxel coordinates
    end
    
    conmask = xyz2mask(mask,XYZvox);   % put points in mask - could weight by Z-scores here, if desired.
    
    % -----------------------------------------------------
    % * convert to density mask
    % * write density mask before thresholding
    % -----------------------------------------------------
    conmask = mask2density(conmask,radius);

    conmask(conmask > 1) = 1;       % limit to max of 1, in case multiple nearby points in same contrast
                                    % max activation for a single
                                    % contrast is 1.

    conmask = conmask .* studyweight;  % sample size weighting by sqrt relative sample size
end
        
if length(varargin) > 0
    return
else
    V.fname = str;
    warning off, spm_write_vol(V,conmask);, warning on
end


return





% -----------------------------------------------------
% * sub-functions
% -----------------------------------------------------

function XYZout = mm2vox(XYZ,M)
% converts a list of coordinates
% XYZ should be 3-column vector [x y z]
% calls tal2vox

XYZ = XYZ';
for i = 1:size(XYZ,2), 
	XYZout(i,:) = tal2vox(XYZ(:,i),M); 
end

XYZout = round(XYZout);

return



% -----------------------------------------------------
function vox=tal2vox(tal,M)
% converts from talairach coordinate to voxel coordinate
% based on variables from SPM.M (passed here for 
% faster operation)
% e.g., foo=tal2vox([-30 28 -30], VOL)
% from Russ Poldrack's spm ROI utility 

vox=[0 0 0];
vox(1)=(tal(1)-M(1,4))/M(1,1);
vox(2)=(tal(2)-M(2,4))/M(2,2);
vox(3)=(tal(3)-M(3,4))/M(3,3);

return




function mask = xyz2mask(mask,XYZ)

% -----------------------------------------------------
% * make a mask out of XYZ, in space of P
% * deal with repeated coordinates by adding
% -----------------------------------------------------

try
    ind = sub2ind(size(mask),XYZ(:,1),XYZ(:,2),XYZ(:,3));
catch
    warning('FOUND POINTS OUTSIDE OF IMAGE?  RETURNING EMPTY MASK.');
    XYZ
    ind = [];
end

mask(ind) = 1;

ind = sort(ind);
repeats = ind(find(~diff(ind)));            % index values that are repeated
for i = 1:length(repeats)
    mask(repeats(i)) = mask(repeats(i)) + 1;
end
