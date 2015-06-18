function dist = distancetosurface(coords,vertices)
% Won't work unless all vertices are on the actual surface!
%
% This function is old and may work for its intended purpose,
% but it is not a general function for calculating distance from a set of
% surface vertices.  therefore, use dist_to_surface.m

for i = 1:size(coords,1)
	
coord = coords(i,:);

% cut down number of vertices by getting into the right quadrant.
if coord(1) >= 0, v2 = vertices(vertices(:,1)>=0,:);,else v2 = vertices(vertices(:,1)<0,:);,end
if coord(2) >= 0, v3 = v2(v2(:,2)>=0,:);,else v3 = v2(v2(:,2)<0,:);,end
if coord(3) >= 0, v4 = v3(v3(:,3)>=0,:);,else v4 = v3(v3(:,3)<0,:);,end


% define center of volume and adjust coordinates
center(1) = 0;
center(2) = (max(vertices(:,2)) + min(vertices(:,2))) /2;
center(3) = (max(vertices(:,3)) + min(vertices(:,3))) /2;
centeredcoord = coord - center;		% now center is the new origin
%adjv4(:,1) = v4(:,1) - center(1);
%adjv4(:,2) = v4(:,2) - center(2);
%adjv4(:,3) = v4(:,3) - center(3);
cdist = (centeredcoord(1) - v4(:,1)).^2 + (centeredcoord(2) - v4(:,2)).^2 + (centeredcoord(3) - v4(:,3)).^2;
cdist = min(cdist);
dist(i,1) = sqrt(cdist(1));

end
return
