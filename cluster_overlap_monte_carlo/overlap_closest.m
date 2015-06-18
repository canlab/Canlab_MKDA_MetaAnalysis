function [clout,o1,ovec1,ovec2,xyzo] = overlap_closest(cla1,cla2,k,varargin)
% function [clout,o1,ovec1,ovec2,xyzo] = overlap_closest(cla1,cla2,k,[v1,v2])
%
% tor wager
% o1 and o2 are overall measures of overlap
% ovec1 and 2 are for each cluster - how many voxels in that 
%   cluster overlap with the other set
% cla1 and 2 are the clusters structures
% v1 and v2 are optional vectors indicating cluster memberships:
% integers code which cluster, for all clusters strung together
% e.g., [1 1 1 1 1 1 2 2 2 2 3 3 3 3 3 ... etc]
% if not included, they will be created
% include them to save time in simulation
%
% this version works on mm coordinates, and
% finds matches within k mm, whereas overlap.m
% works on voxels coords and finds exact matches.
%
% see also overlap.m and find_closest_cluster.m

if length(varargin) > 0
	v1 = varargin{1}; v2 = varargin{2};
else
    v1 = [];
    for i = 1:length(cla1), v1 = [v1; i .* ones(1,size(cla1(i).XYZmm,2))'];, end
    v2 = [];
    for i = 1:length(cla2), v2 = [v2; i .* ones(1,size(cla2(i).XYZmm,2))'];, end 
end
    
xyz1 = round(cat(2,cla1.XYZmm))'; xyz2 = round(cat(2,cla2.XYZmm))';

% exact matches
[xyzoexact,i1,i2] = intersect(xyz1,xyz2,'rows');
xyzo = xyzoexact;
vv1 = v1(i1); vv2 = v2(i2);
disp(['Total exact overlap is ' num2str(size(xyzo,1)) ' voxels'])

% approximate matches - within x mm
for i1 = 1:length(cla1)
    
    clear d, clear wh
    % get distance between closest points
    for i = 1:size(cla1(i1).XYZmm,2)
        x1 = repmat(cla1(i1).XYZmm(:,i)',size(xyz2,1),1);
        tmp = sum((x1 - xyz2) .^2,2) .^ .5; % distances to each in set 2
        d(i) = min(tmp);
        
        t2 = v2(find(tmp==min(tmp)));
        if length(t2) > 1   % break tie - choose largest
            t3 = cat(1,cla2(t2).numVox);
            t3 = find(t3 == max(t3)); t3 = t3(1);
            t2 = t2(t3);
        end
        
        wh(i) = t2;    % contains cluster index of cluster closest to vox in set 1
    end
    t2 = wh(find(d == min(d)));
    if length(t2) > 1   % break tie - choose largest
            t3 = cat(1,cla2(t2).numVox);
            t3 = find(t3 == max(t3)); t3 = t3(1);
            t2 = t2(t3);
    end
        
    cla1(i1).closest =t2;
    cla1(i1).distance = min(d);
    
end
    
omit = cat(1,cla1.distance);
for i=1:length(cla1),cla1(i).index = i;,end
clout = cla1; clout(omit > k) = [];

for i = 1:length(clout)
    tmp = intersect(clout(i).XYZmm',cla2(clout(i).closest).XYZmm','rows');
    clout(i).exact = size(tmp,1);
    
    %if clout(i).distance == 0, keyboard,end
end
    
% print table
[clout] = cluster_ba(clout,1:length(clout));

fprintf(1,'\noverlap_closest.m output table: comparing sets of clusters\nMatch within %3.0f mm\n\n',k)
fprintf(1,'\t\t\tOverlap\t\tSet 1\t\t\t\t\tSet 2\t\n')
fprintf(1,'Index\tIndex1\tIndex2\tx\ty\tz\tx\ty\tz\tVoxels\tMaxZ\tx\ty\tz\tVoxels\tMaxZ\tDistance\tExact Overlap\tBrodmann\n')
for i1 = 1:length(clout)
    
    max1 = clout(i1).index; max2 = clout(i1).closest;
    
    fprintf(1,'%3.0f\t%3.0f\t%3.0f\t%3.0f\t%3.0f\t%3.0f\t%3.0f\t%3.0f\t%3.0f\t%3.0f\t%3.2f\t%3.0f\t%3.0f\t%3.0f\t%3.0f\t%3.2f\t%3.0f\t%3.0f\t%s\n', ...
        i1,max1,max2,clout(i1).mm_center(1),clout(i1).mm_center(2),clout(i1).mm_center(3), ...
        cla1(max1).mm_center(1),cla1(max1).mm_center(2),cla1(max1).mm_center(3),cla1(max1).numVox,max(cla1(max1).Z), ...
        cla2(max2).mm_center(1),cla2(max2).mm_center(2),cla2(max2).mm_center(3),cla2(max2).numVox,max(cla2(max2).Z), ...
        clout(i1).distance,clout(i1).exact,clout(i1).BAstring)
    
end

return