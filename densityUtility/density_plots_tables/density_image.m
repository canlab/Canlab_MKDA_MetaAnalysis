function density_image(varargin)
% function density_image([fname],[radius],[threshold],[plot types])
%
% tor wager, 4/19/03
% 
% This function plots line graphs, a montage, and a surface rendering
% of a density map image, thresholded at a level you choose
%
% fname:        must be an analyze image file containing density values
%               in points / mm3
% radius:       the radius in mm of point smoothing used to create the image
% threshold:    the minimum number of peaks in the radius that will appear
%               on the plots
% plot types:   a vector of integers specifying which plots to produce.
%               if this vector contains an integer from 1 to 5, the corresponding
%               plot will be produced
%               1: line plot of number of significant voxels by threshold
%               2: montage of slices at specified threshold
%               3: surface map of significant densities at threshold
%               4: left medial surface map at threshold
%               5: right medial surface map at threshold
%               default is [1 2 3]
%
% ALSO: surface maps plot vertices within 5 mm of a suprathreshold voxel in color
% img files and clusters structure are saved for future use; clusters in pwd, img in
% original directory it came from.
%
% example:
% fname = 'dens_i_density_rad15_step1.img';
% density_image(fname,15,10,[1 2 3 4])

if length(varargin) == 0, fname = spm_get(1,'Choose density image file');, 
else, fname = varargin{1};, end
if length(varargin) < 2, r = input('Enter radius in mm:');, 
else, r = varargin{2};, end
if length(varargin) < 3, n = input('Enter threshold in # of peaks within radius:');, 
else, n = varargin{3};, end
if length(varargin) < 4, ptype = [1 2 3];, 
else, ptype = varargin{4};, end

r = 15; v = 4*pi*r^3 / 3;

V = spm_vol(fname); vol = spm_read_vols(V);

vol = vol .* v; vol = round(vol);

if any(ptype==1)
    figure('Color','w')
    for i = 0:max(vol(:)), rr(i+1) = sum(vol(:) > i);, end
    figure('Color','w');plot(0:max(vol(:)),rr,'ro-','LineWidth',2),
    ylabel('Number of voxels above threshold','FontSize',18),xlabel('Threshold (number of points < 15 mm)','FontSize',18)
    set(gca,'FontSize',18)
    title(fname)
end

vol(vol < n) = 0;
%sum(vol(:) > 0)
[d,fbase,fext] = fileparts(fname);
V.fname = fullfile(d,[fbase '_r' num2str(r) '_thr' num2str(n) fext]);
spm_write_vol(V,vol);
clusters = mask2clusters(V.fname);
fprintf(1,'Found %6.0f voxels above threshold of %3.0f\n',sum(cat(1,clusters.numVox)),n)

eval(['save clusters_' fbase '_r' num2str(r) '_thr' num2str(n) ' clusters'])

if any(ptype==2)
% montage

    montage_clusters([],clusters,[2 2]);
    set(gca,'Position',[0.1300    0.1100    0.7750    0.8150])
    xlabel(['Number of points < ' num2str(r) ' mm'],'FontSize',18)
end

if any(ptype==3)
% surface

    cluster_surf(clusters,which('surf_single_subj_T1_gray.mat'),5,'heatmap')
    
end

if any(ptype==4)
% surface

    cluster_surf(clusters,which('surf_single_subj_grayL.mat'),5,'heatmap')
    
end

if any(ptype==5)
% surface

    cluster_surf(clusters,which('surf_single_subj_grayR.mat'),5,'heatmap')
    
end

return


