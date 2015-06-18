function  clusters = smooth_dens_and_mask(P,varargin)
% function  clusters = smooth_dens_and_mask(P,[P2: extract img])
%
% tor wager
% 
% P is an input file, usually a density results mask
% this function smoothes the image with a 3 mm kernel
% and thresholds the results at i1 > 0, so that
% you end up with a mask.  The mask is converted
% into a clusters structure and saved.
% 
% the output files are written with the prefix of sf_*
%
% if a second filename is entered (P2),
% the all_data and Z fields contains the values contained in P2
% clusters = smooth_dens_and_mask('i_density_res_mask.img');
% for i = 1:length(clusters),clusters(i).Z = clusters(i).all_data;,end
% montage_clusters([],clusters,[2 2])
%

V = spm_vol(P);
V.fname = ['sf_' V.fname];

spm_smooth(P,V.fname,3);

V = spm_vol(V.fname); 
v = spm_read_vols(V);

v  = v > 0;

spm_write_vol(V,v);

if length(varargin) > 0
    clusters = mask2clusters(V.fname,varargin{1});
    for i = 1:length(clusters),clusters(i).Z = clusters(i).all_data;,end
else
    clusters = mask2clusters(V.fname);
end

return

