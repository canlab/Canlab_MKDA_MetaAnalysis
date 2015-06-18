function [cl1,cl2,dmt,dmt2,d1,d2,mask,mask2] = density_diff(XYZ,u,XYZ2,u2,varargin)
% function [cl1,cl2,dmt,dmt2,d1,d2,mask,mask2] = density_diff(XYZ,u,XYZ2,u2,[radius_mm],[first outfile name],[2nd outfile name])
% 
% function gets density from a list of XYZ coordinates
% in one condition of a meta-analysis
% and makes and saves clusters
%
% XYZ is 3-column vector of coordinates
%
% Tor Wager 2/18/02

P = ['brain_avg152T1.img'];     % 2 mm voxels
P = 'scalped_avg152T1_graymatter_smoothed.img';
%radius_mm = 10;

if length(varargin) > 0, radius_mm = varargin{1}; 
else, radius_mm = 10;
end
disp(['Radius is ' num2str(radius_mm) ' mm.'])

if length(varargin) > 1, fname = varargin{2}; 
else, fname = input('Enter text string tag for first set of XYZ ','s');
end

if length(varargin) > 2, fname2 = varargin{3}; 
else, fname2 = input('Enter text string tag for second set of XYZ ','s');
end

mask_file = P;
t1 = clock;
%fname = input('Enter text string tag for first set of XYZ ','s');
%fname2 = input('Enter text string tag for second set of XYZ ','s');
f1 = [fname '-' fname2];
f2 = [fname2 '-' fname];
[d f1] = fileparts(f1);, f1 = [f1 '.img'];
[d f2] = fileparts(f2);, f2 = [f2 '.img'];
fprintf(1,'\tdensity_diff.m writing: \n\t%s\n\t%s\n',f1,f2)

% -----------------------------------------------------
% * load standard brain
% -----------------------------------------------------
P = which(P);
V = spm_vol(P);
mask = zeros(V.dim(1:3));
mask2= mask;

voxsize = diag(V(1).mat)';
voxsize = voxsize(1:3);
radius = radius_mm ./ mean(voxsize);
sphere_vol = 4 * pi * radius_mm ^ 3 / 3;
XYZ = mm2vox(XYZ,V.mat);
XYZ2 = mm2vox(XYZ2,V.mat);

% -----------------------------------------------------
% * make a mask out of XYZ, in space of P
% * deal with repeated coordinates by adding
% -----------------------------------------------------

ind = sub2ind(size(mask),XYZ(:,1),XYZ(:,2),XYZ(:,3));
mask(ind) = 1;

ind = sort(ind);
repeats = ind(find(~diff(ind)));            % index values that are repeated
for i = 1:length(repeats)
    mask(repeats(i)) = mask(repeats(i)) + 1;
end

ind = sub2ind(size(mask2),XYZ2(:,1),XYZ2(:,2),XYZ2(:,3));
mask2(ind) = 1;

ind = sort(ind);
repeats = ind(find(~diff(ind)));            % index values that are repeated
for i = 1:length(repeats)
    mask2(repeats(i)) = mask2(repeats(i)) + 1;
end

% -----------------------------------------------------
% * convert to density mask
% * write density mask before thresholding
% -----------------------------------------------------
dm = mask2density(mask,radius,[],sphere_vol);
dm2 = mask2density(mask2,radius,[],sphere_vol);

d1 = dm - dm2;
d2 = dm2 - dm;
clear dm, clear dm2

V.fname = ['dens_' f1];
spm_write_vol(V,d1);

V.fname = ['dens_' f2];
spm_write_vol(V,d2);

%figure; hist(d1(d1>0),50)
%if ~isempty(u),hold on; plot([u u],get(gca,'YLim'),'r','LineWidth',2),end

%figure; hist(d2(d2>0),50)
%if ~isempty(u2),hold on; plot([u2 u2],get(gca,'YLim'),'r','LineWidth',2),end

% -----------------------------------------------------
% * make and write filtered density mask
% -----------------------------------------------------
dmt = maskImg(d1,u,Inf);   
[d,fname,ext] = fileparts(f1);
V.fname = ['dens_' fname '_filtered' ext];
spm_write_vol(V,dmt);

% -----------------------------------------------------
% * get clusters from density mask for imaging
% -----------------------------------------------------
CLU = mask2struct(V.fname,u,0);
if size(CLU.XYZ,2) > 0
    clusters = tor_extract_rois([],CLU,CLU);
else
    disp(['1st XYZ: no significant results.  max is ' num2str(max(max(max(d1))))])
    clusters = [];
end
eval(['save ' fname '_clusters CLU clusters'])
cl1 = clusters;

% -----------------------------------------------------
% * make and write filtered density mask
% -----------------------------------------------------
dmt2 = maskImg(d2,u2,Inf);   
[d,fname,ext] = fileparts(f2);
V.fname = ['dens_' fname '_filtered' ext];
spm_write_vol(V,dmt2);

% -----------------------------------------------------
% * get clusters from density mask for imaging
% -----------------------------------------------------
CLU = mask2struct(V.fname,u2,0);
if size(CLU.XYZ,2) > 0
    clusters = tor_extract_rois([],CLU,CLU);
else
    disp(['2nd XYZ: no significant results.  max is ' num2str(max(max(max(d2))))])
    clusters = [];
end
eval(['save ' fname '_clusters CLU clusters'])
cl2 = clusters;


%try
%    P = spm_get(1,'*img','Select overlay');
%    if ~isempty(cl1)
%        cluster_table(cl1)
%        montage_clusters(P,cl1)
%    end
%    if ~isempty(cl2)
%        cluster_table(cl2)
%        montage_clusters(P,cl2)
%    end
%catch
%    disp('Problem making montages...skipping.')
%end

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



