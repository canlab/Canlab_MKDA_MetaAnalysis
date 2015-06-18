function [o1,ovec1,ovec2,xyzo] = overlap(cla1,cla2,varargin)
% function [o1,ovec1,ovec2,xyzo] = overlap(cla1,cla2,[v1,v2])
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
    
xyz1 = round(cat(2,cla1.XYZ))'; xyz2 = round(cat(2,cla2.XYZ))';

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
return