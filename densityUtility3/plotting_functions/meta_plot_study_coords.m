function meta_plot_study_coords(DB, studyname, wh_slice)
% meta_plot_study_coords(DB, studyname, wh_slice)
%
% Plot points on brain slice and then show orthviews for one study
%
% Tor Wager, March 2007
%
% wh_slice can be missing or empty: picks modal slice (most points)

if nargin < 3, wh_slice = []; end

wh = find(strcmp(DB.Study, studyname));
xyz = [DB.x(wh) DB.y(wh) DB.z(wh)];

%%%
[handles, wh_slice, my_z] = plot_points_on_slice(xyz, 'marker', 'o','close_enough',10,'slice',wh_slice);
axis off
set(gcf,'Color','w')
scn_export_papersetup;
saveas(gcf,[studyname '_slice' num2str(wh_slice) '_pts'],'png')

%%% put spheres on brain 
maskname = which('scalped_avg152T1_graymatter_smoothed.img');
[V, maskdata] = iimg_read_img(maskname, 1);
r = 5; % in voxels

xyz = [DB.x(wh) DB.y(wh) DB.z(wh)];
xyzvox = mm2voxel(xyz, V.mat, 1);
ivector = iimg_xyz2spheres(xyzvox, V.xyzlist, r);
modes = mode(xyzvox);
slicedata = iimg_reconstruct_3dvol(ivector, V);
cl = iimg_indx2clusters(ivector, V, 0);
cluster_orthviews(cl,{[0 0 0]});

modes_mm = mode(xyz);
spm_orthviews('Reposition',[modes_mm(1:2) my_z]);

fh = findobj('Tag', 'Graphics');
set(fh, 'InvertHardCopy', 'off');
set(fh, 'Color', 'white');
figure(fh);
scn_export_papersetup;
saveas(gcf,[studyname '_slice' num2str(wh_slice)],'png')

end


% % 
% % %%
% % 
% % 
% % %%P = 'scalped_avg152T1_graymatter_smoothed.img';
% % P = which(P);
% % radius_mm = 10;
% % V = spm_vol(P);
% % mask = zeros(V.dim(1:3));
% % voxsize = diag(V(1).mat)';
% % voxsize = voxsize(1:3);
% % radius = radius_mm ./ mean(abs(voxsize));
% % sphere_vol = 4 * pi * radius_mm ^ 3 / 3;
% % edit Meta_Activation_FWE.m
% % edit meta_prob_activation.m
% % [V, maskdata] = iimg_read_img(DB.maskname, 1);
% % xyzlist = V.xyzlist;
% % maskdata = iimg_reconstruct_3dvol(maskdata, V);
% % [V, maskdata] = iimg_read_img(maskname, 1);
% % mask = Pl
% % mask = P;
% % [V, maskdata] = iimg_read_img(maskname, 1);
% % maskname = P;
% % [V, maskdata] = iimg_read_img(maskname, 1);
% % xyzlist = V.xyzlist;
% % maskdata = iimg_reconstruct_3dvol(maskdata, V);
% % r
% % r = radius
% % 
% % %% put spheres on brain 
% % maskname = which('scalped_avg152T1_graymatter_smoothed.img');
% % [V, maskdata] = iimg_read_img(maskname, 1);
% % r = 5; % in voxels
% % 
% % xyz = [DB.x(wh) DB.y(wh) DB.z(wh)];
% % xyzvox = mm2voxel(xyz, V.mat, 1);
% % ivector = iimg_xyz2spheres(xyzvox, xyzlist, r);
% % modes = mode(xyzvox);
% % slicedata = iimg_reconstruct_3dvol(ivector, V);
% % cl = iimg_indx2clusters(ivector, V, 0);
% % modes_mm = mode(xyz);
% % spm_orthviews('Reposition',modes_mm);
% % 
% % % just the points
% % ivector = iimg_xyz2spheres(xyzvox, xyzlist, 1);
% % cl_pts = iimg_indx2clusters(ivector, V, 0);
% % cluster_orthviews(cl_pts,{[1 0 0]});
% % 
% % %% draw brain slice
% % wh_slice = 29;
% % 
% % ovl = which('scalped_avg152T1.img'); 
% % ovlV = spm_vol(ovl);
% % v = spm_read_vols(ovlV);
% % 
% % mm = voxel2mm([1 1 1; size(v)]', ovlV.mat);
% % xcoords = linspace(mm(1,1), mm(1,2), size(v, 1));
% % ycoords = linspace(mm(2,1), mm(2,2), size(v, 2));
% % zcoords = linspace(mm(3,1), mm(3,2), size(v, 3));
% % [X, Y] = meshgrid(xcoords, ycoords); Z = v(:,:,wh_slice)';
% % figure; contour(X, Y, Z); colormap(gray);
% % axis image
% % 
% % %% put points on slice
% % close_enough = 10;   % in mm
% % 
% % % get z coord in mm of this slice
% % tmp_xyzmm = voxel2mm([1 1 29]', ovlV.mat);
% % my_z = tmp_xyzmm(3);
% % 
% % wh = find(abs(xyz(:,3) - my_z) <= close_enough);
% % hold on; plot(xyz(wh,1), xyz(wh,2), 'ko', 'MarkerFaceColor','k');

