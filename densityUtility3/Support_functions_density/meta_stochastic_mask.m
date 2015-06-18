function indic = meta_stochastic_mask(xyzlist,n,r)
% indic = meta_stochastic_mask(xyzlist,n,r)
%
% Makes convolved spheres at n randomized locations within analysis mask
% Works in indicator (indic) format for speed
%
% tor wager

v = size(xyzlist,1);

wh = ceil(rand(n,1) .* v);
xyzvox = xyzlist(wh,:);

%indic = meta_fast_sphere_conv(xyzlist,xyzvox,r);
indic = iimg_xyz2spheres(xyzvox,xyzlist,r);
 
return