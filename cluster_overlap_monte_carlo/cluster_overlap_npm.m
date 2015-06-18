function [OUT] = cluster_overlap_npm(cl1,cl2,varargin)
% function [OUT] = cluster_overlap_npm(cl1,cl2,[mask_fname],[iterations])
% tor wager
% takes 2 clusters structures and determines overlap
% - then runs Monte Carlo simulation randomizing cluster
%   centers within a mask
%
% Mask and input files must have the same voxel sizes and dimensions!
% They must be in the same space!
%
% OUT has stuff for stats and overlapping voxels
% see overlap_table.m

% ----------------------------------------------------
% * optional argument setup
% ----------------------------------------------------

if length(varargin) > 0
	P = varargin{1};
else
    P = ['brain_avg152T1.img'];     % 2 mm voxels
end

if length(varargin) > 0, niter = varargin{2}; else niter = 5000; end

% ----------------------------------------------------
% * get list of all coordinates in brain mask
%   vol: row, col, array = x, y, z
% ----------------------------------------------------
P = which(P);
    V = spm_vol(P);
    vol = spm_read_vols(V);
    voxsize = diag(V(1).mat)';
    voxsize = voxsize(1:3);
    
[x,y,z] = ind2sub(size(vol),find(vol));     % find values > 0 in vol
allXYZ = [x y z];                           % XYZ is in 3 columns in this function
%Q       = ones(1,size(allXYZ,1));
%slices  = allXYZ(:,3)';                      % used to pass only certain slices to subfunctions

% ----------------------------------------------------
% * prepare cluster info for faster iteration
%   strip clusters of extraneous fields
%   calculate the amount of overlap
% ----------------------------------------------------
l1 = length(cl1); l2 = length(cl2);
v1 = [];
for i = 1:length(cl1), 
    cla1(i).XYZ = cl1(i).XYZ;
    cla1(i).center = mean(cl1(i).XYZ,2)';
    cla1(i).mm_center = cl1(i).mm_center;
    v1 = [v1 i .* ones(1,size(cla1(i).XYZ,2))];, 
end
cl1 = cla1; clear cla1;

v2 = [];
for i = 1:length(cl2), 
    cla2(i).XYZ = cl2(i).XYZ;
    cla2(i).center = mean(cl2(i).XYZ,2)';
    cla2(i).mm_center = cl2(i).mm_center;
    v2 = [v2 i .* ones(1,size(cla2(i).XYZ,2))];, 
end 
cl2 = cla2; clear cla2;
OUT.cl1 = cl1; OUT.cl2 = cl2;
[OUT.obs_o1,OUT.obs_ovec1,OUT.obs_ovec2,OUT.xyzo,OUT.cent1,OUT.cent2] = overlap(cl1,cl2,v1,v2);

% change voxel centers to mm centers
for i = 1:size(OUT.cent1,1), if ~isnan(OUT.cent1(i,1)),OUT.cent1(i,:)=voxel2mm(OUT.cent1(i,:)',V.mat);,end,end
for i = 1:size(OUT.cent2,1), if ~isnan(OUT.cent2(i,1)),OUT.cent2(i,:)=voxel2mm(OUT.cent2(i,:)',V.mat);,end,end

% get BAs
warning off
for i = 1:size(OUT.cent1,1), if ~isnan(OUT.cent1(i,1)),tmp=mni2BA(OUT.cent1(i,:)');,OUT.ba1{i} = tmp;,end,end
for i = 1:size(OUT.cent2,1), if ~isnan(OUT.cent2(i,1)),tmp=mni2BA(OUT.cent2(i,:)');,OUT.ba2{i} = tmp;,end,end
warning on


t1 = clock;
for i = 1:niter

    % ----------------------------------------------------
    % * pick random voxel coords
    % * adjust each cluster to centers
    % * get overlap (els are overlap for each cluster)
    % ----------------------------------------------------
    
    [c1,c2] = getcenters(l1,l2,allXYZ);
    [cla1,cla2] = movecenters(cl1,cl2,c1,c2);
    [o1(i),ovec1(i,:),ovec2(i,:)] = overlap(cla1,cla2);
    
    if mod(i,500)==0, fprintf(1,'.'), end
end

fprintf(1,' done (in %4.0f s)!\n',etime(clock,t1))
OUT.o1 = o1; OUT.ovec1 = ovec1; OUT.ovec2 = ovec2;

% ----------------------------------------------------
% * print table of results
% ----------------------------------------------------

overlap_table(OUT)


return


% ----------------------------------------------------
% * 
% * Sub-functions
% * 
% ----------------------------------------------------
    

function [c1,c2] = getcenters(l1,l2,XYZ)

tmp = ceil(rand(l1,1) .* size(XYZ,1));
c1 = XYZ(tmp,:);

tmp = ceil(rand(l2,1) .* size(XYZ,1));
c2 = XYZ(tmp,:);

return




function [cla1,cla2] = movecenters(cl1,cl2,c1,c2)

for i = 1:length(cl1)
    tmp = cl1(i).XYZ' - repmat(cl1(i).center,size(cl1(i).XYZ,2),1);
    cla1(i).XYZ = round(tmp + repmat(c1(i,:),size(tmp,1),1))';
end

for i = 1:length(cl2)
    tmp = cl2(i).XYZ' - repmat(cl2(i).center,size(cl2(i).XYZ,2),1);
    cla2(i).XYZ = round(tmp + repmat(c2(i,:),size(tmp,1),1))';
end    

return

    

function [o1,ovec1,ovec2,xyzo,cent1,cent2] = overlap(cla1,cla2,varargin)
% function [o1,o2,ovec1,ovec2,xyzo,centers1,centers2] = overlap(cla1,cla2,[v1,v2])
%
% tor wager
% o1 and o2 are overall measures of overlap
% ovec1 and 2 are for each cluster
% cla1 and 2 are the clusters structures
% v1 and v2 are optional vectors indicating cluster memberships:
% integers code which cluster, for all clusters strung together
% e.g., [1 1 1 1 1 1 2 2 2 2 3 3 3 3 3 ... etc]
% if not included, they will be created
% include them to save time in simulation

if length(varargin) > 0
	v1 = varargin{1}; v2 = varargin{2};
else
    v1 = [];
    for i = 1:length(cla1), v1 = [v1 i .* ones(1,size(cla1(i).XYZ,2))];, end
    v2 = [];
    for i = 1:length(cla2), v2 = [v2 i .* ones(1,size(cla2(i).XYZ,2))];, end 
end
    
xyz1 = cat(2,cla1.XYZ)'; xyz2 = cat(2,cla2.XYZ)';

[xyzo,i1,i2] = intersect(xyz1,xyz2,'rows');
o1 = size(xyzo,1);

% individual clusters

vv1 = v1(i1); vv2 = v2(i2);
if ~isempty(vv1),
    for i = 1:length(cla1), ovec1(i) = sum(vv1 == i);, end                    
else
    ovec1 = zeros(1,length(cla1));
end

if ~isempty(vv2),
    for i = 1:length(cla2), ovec2(i) = sum(vv2 == i);, end
else
    ovec2 = zeros(1,length(cla2));
end



    if nargout > 4  % save overlap centers
        for i = 1:length(cla1),
            [tmp,i1,i2] = intersect(cla1(i).XYZ',xyz2,'rows');
            cent = mean(tmp); if isempty(cent), cent = [NaN NaN NaN];,end
            cent1(i,:) = cent;
        end
    end
    
    if nargout > 4  % save overlap centers
        for i = 1:length(cla2),
            [tmp,i1,i2] = intersect(cla2(i).XYZ',xyz1,'rows');
            cent = mean(tmp); if isempty(cent), cent = [NaN NaN NaN];,end
            cent2(i,:) = cent;
        end
    end
    

return

    

    