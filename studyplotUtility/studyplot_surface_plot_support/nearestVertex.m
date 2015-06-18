function [surfacevert,dist] = nearestVertex(coord,vertices,varargin)

if nargin > 3, plotit = 1;,else plotit = 0;,end

% cut down number of vertices by getting into the right quadrant.
if coord(1) >= 0, v2 = vertices(vertices(:,1)>=0,:);,else v2 = vertices(vertices(:,1)<0,:);,end
if coord(2) >= 0, v3 = v2(v2(:,2)>=0,:);,else v3 = v2(v2(:,2)<0,:);,end
if coord(3) >= 0, v4 = v3(v3(:,3)>=0,:);,else v4 = v3(v3(:,3)<0,:);,end

if size(v4,1) > 10000, status = 1;,else status = 0;,end

% define center of volume and adjust coordinates
center(1) = 0;
center(2) = (max(vertices(:,2)) + min(vertices(:,2))) /2;
center(3) = (max(vertices(:,3)) + min(vertices(:,3))) /2;
centeredcoord = coord - center;		% now center is the new origin
%adjv4(:,1) = v4(:,1) - center(1);
%adjv4(:,2) = v4(:,2) - center(2);
%adjv4(:,3) = v4(:,3) - center(3);


% extrapolate point to somewhere outside the brain, so we can find the CLOSEST point on the surface
scale = 200 / sqrt((coord(1)-center(1))^2 + ((coord(2)-center(2))^2) + ((coord(3)-center(3))^2));	% the distance from the center should be 200...
ocoord = scale * centeredcoord;
adjcoord = centeredcoord - ocoord;		% now ocoord, outside the brain, is the new origin
xadj = center(1) + ocoord(1);
yadj = center(2) + ocoord(2);
zadj = center(3) + ocoord(3);
adjv4(:,1) = v4(:,1) - xadj;			% adjust vertices to be centered on ocoord
adjv4(:,2) = v4(:,2) - yadj;
adjv4(:,3) = v4(:,3) - zadj;


unit = adjcoord ./ norm(adjcoord);

% choose a column of points close to the line
% ------------------------------------------- %
for i = 1:size(adjv4,1)
   % choose small orthod = sum((vertex - coord of prjn)^2) - which is the sum of squares of the orthogonal projection onto the line.
   orthod(i) = sum( (adjv4(i,:) - (dot(unit,adjv4(i,:))*unit)) .^2) ;
   if status & mod(i,5000) == 0, disp(['done ' num2str(i)]),end   
end

% restrict points to those than lie w/i 1 mm of line
cylinderrad = 16;	% radius of acceptable squared dist to line
if nargin>2,cylinderrad = varargin{1};,end

if sqrt(min(orthod)) > 4,	disp(['Warning: for point [' num2str(coord) '] closest point to prjn is ' num2str(sqrt(min(orthod))) ' mm']),end
adjv5 = adjv4(orthod < cylinderrad,:);
v5 = v4(orthod < cylinderrad,:);

if size(adjv5,1) == 0,
    disp(['WARNING! No points within cylinder radius for ' num2str(coord)])
    surfacevert = [0 0 0]; dist = Inf;
    return
end

% for points close to line, minimize projection of vertex on line
% ------------------------------------------- %
for i = 1:size(adjv5,1)
   % compute dot product, which is proportional to projection but faster
	dotp(i) = dot(unit,adjv5(i,:));
end

surfacevert = v5(dotp == min(dotp),:);		% minimize the dot product - the closest vertex
sufacevert = surfacevert(1,:);
dist = sqrt((coord(1)-surfacevert(1))^2 + ((coord(2)-surfacevert(2))^2) + ((coord(3)-surfacevert(3))^2));

if plotit & dist < 20,
	plot3([center(1) coord(1)],[center(2) coord(2)],[center(3) coord(3)],'-go','MarkerSize',8,'LineWidth',2)
	plot3([ocoord(1)],[ocoord(2)],[ocoord(3)],'-m+','MarkerSize',8,'LineWidth',2)
	plot3([ocoord(1) surfacevert(1)],[ocoord(2) surfacevert(2)],[ocoord(3) surfacevert(3)],'-b','MarkerSize',8,'LineWidth',2)
	%text(ocoord(1), ocoord(2), ocoord(3),num2str(coord))
	text(ocoord(1), ocoord(2), ocoord(3),[' ' num2str(dist)])
end

return
