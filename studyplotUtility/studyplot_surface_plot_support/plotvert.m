function h = plotvert(vert,ptcolor,varargin)
% function h = plotvert(vert,pointcolor,vertices[opt],lncolor [opt])
%
% draws lines optionally, if you specify verties list and linecolor

h(1) = plot3(vert(1),vert(2),vert(3),ptcolor,'MarkerSize',4,'MarkerFaceColor',ptcolor(1));

if nargin > 3

vertices = varargin{1};
lncolor = varargin{2};

minx = min(vertices(:,1)) - 75;
miny = min(vertices(:,2)) - 75;
minz = min(vertices(:,3)) - 75;
maxx = max(vertices(:,1));
maxy = max(vertices(:,2));
maxz = max(vertices(:,3));

% lines
h = plot3([minx vert(1)],[vert(2) vert(2)],[vert(3) vert(3)],lncolor,'LineWidth',2);
h = plot3(minx,vert(2),vert(3),ptcolor,'MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor',ptcolor(1));

h = plot3([vert(1) vert(1)],[miny vert(2)],[vert(3) vert(3)],lncolor,'LineWidth',2);
h = plot3(vert(1),miny,vert(3),ptcolor,'MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor',ptcolor(1));

h = plot3([vert(1) vert(1)],[vert(2) vert(2)],[minz vert(3)],lncolor,'LineWidth',2);
h = plot3(vert(1),vert(2),minz,ptcolor,'MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor',ptcolor(1));

end

return
