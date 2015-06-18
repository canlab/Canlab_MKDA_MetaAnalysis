function out = tor_conv_sphere(mask,vox_radius,mmethod,varargin)
%function out = tor_conv_sphere(mask,vox_radius,mmethod,[norm factor])
% norm factor is optional.  scaling factor to divide all values by
% in spherical convolution.  
% used because the convolution does not return a max height of 1 given
% one input point.  
% you'd put in 1 / the max height for a convolution of a single point as
% the norm factor.
%
% returns either a box convolution, faster than smooth3,
% or a truncated gaussian convolution

mask = double(mask);

switch mmethod
case 'box' 
    %x  = floor(s(1)); x = [-x:x]; x(:) = 1;

    x = ones(1,floor(vox_radius) * 2 + 1);
    i  = (length(x) - 1)/2;

    mul = length(x);
    
    y = x; z = x;
    j = i; k = i;

    
    
case 'sphere'    

    % from spm_smooth
    % but assumes isotropic voxels, and input of voxel radius!
    
    % s  = s./VOX;					% voxel anisotropy
    s = [vox_radius vox_radius vox_radius];
    s  = max(s,ones(size(s)));			% lower bound on FWHM
    s  = s/sqrt(8*log(2));				% FWHM -> Gaussian parameter

    x  = round(6*s(1)); x = [-x:x];
    y  = round(6*s(2)); y = [-y:y];
    z  = round(6*s(3)); z = [-z:z];
    x  = exp(-(x).^2/(2*(s(1)).^2)); 
    y  = exp(-(y).^2/(2*(s(2)).^2)); 
    z  = exp(-(z).^2/(2*(s(3)).^2));
    %x  = x/sum(x);
    %y  = y/sum(y);
    %z  = z/sum(z);

    mul = 1; % 
    if length(varargin) > 0, mul = varargin{1};, end
    
    % make max height 1 so that anything less than .5 is less than
    % FWHM, so is outside radius for spherical convolution
    x  = x/max(x);
    y  = y/max(y);
    z  = z/max(z);
    
    % figure out threshold such that 2 * vox_radius voxels are
    % above this threshold on the gaussian convolution
    % this let's us set a threshold for determining the boundaries
    % of the sphere at 2*vox_radius
    % ASSUMING that the max height in the convolved volume for
    % a single point is 1.  (this is why we use "mul").
    
    tmpx = sort(x);
    thresh = tmpx(end-2*vox_radius);
    
    i  = (length(x) - 1)/2;
    j  = (length(y) - 1)/2;
    k  = (length(z) - 1)/2;

otherwise
    disp('Unknown mmethod - choose box or sphere')
end


out = zeros(size(mask));
spm_conv_vol(mask,out,x,y,z,-[i,j,k]);
    
out = out .* mul;

if strcmp(mmethod,'sphere')
    % threshold at thresh to get voxels w/i vox_radius FWHM
    out(out < thresh-.00000001) = 0;
end
    
return


