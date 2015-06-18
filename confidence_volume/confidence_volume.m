function results = confidence_volume(x,varargin)
% function results = confidence_volume(x,color [opt])
% x is 3-column [x y z] matrix of coordinates
% plots if color argument is specified
%
% Tor Wager and Robert Welsh, 10/15/01
% 
% Modified by Tor Wager, 11/16/02 
% to use the method of Johnson & Wichern, 
% Applied Multivariate Statistical Analysis, p. 237, ex. 5.3
%
% Minor modification/documentation by Tor Wager, Aug 2010
% Now correctly indicates that permutation test is used for p-values.
% See help conf_region for details of the test.

if nargin > 1, mycolor = varargin{1};  doplot = 1; end
if nargin > 2, dopts = 0; else dopts = 1; end

[ax,ci,cen,ub,lb,S,e,lam,Fm,Fc,F,pval, msb] = conf_region(x);
ci = diag(ci);


% rows of ub are x,y,z coordinates
% cols are axis 3, axis 2, axis 1 along axes of greatest variation


% ----------------------------------------------------------------------
% now move the boundaries from center 0,0,0 back to original center.
% ----------------------------------------------------------------------

for i = 1:size(ub,2)
    ub(:,i) = ub(:,i) + cen;
    lb(:,i) = lb(:,i) + cen;
end


% ----------------------------------------------------------------------
% get pointlist of coordinates in conf vol
% ----------------------------------------------------------------------
vox = ellipsoidMask(ci(1),ci(2),ci(3));
respts = rotateByEigen(vox,[],[],e);
XYZ = respts.XYZ;
for i = 1:size(XYZ,1)
    XYZ(i,:) = XYZ(i,:) + cen';
end

% ----------------------------------------------------------------------
% image the ellipsoid and rotate to proper position
% need to rotate ellipse xyz into frame of reference of new axes
% then shift ellipse to actual xyz center
% ----------------------------------------------------------------------
[X,Y,Z] = ellipsoid(0,0,0,ci(1),ci(2),ci(3),10);
results = rotateByEigen(X,Y,Z,e);

results.xP = results.xP + cen(1);
results.yP = results.yP + cen(2);
results.zP = results.zP + cen(3);

results.F = F; 
results.F_descrip = 'F-value for test against origin; needs df check. Not used for P-values';

results.msb = msb;
results.msb_descrip = 'Euclidean distance of average [x,y,z] vector from origin; used for P-values';
results.pval = pval;

if doplot
    hold on
    results.p = surf(results.xP,results.yP,results.zP,'EdgeColor','none','FaceColor',mycolor(1),'FaceAlpha',.6);
    % next line plots all points in x
    if dopts
        hold on; for i = 1:size(x,1),plot3(x(i,1),x(i,2),x(i,3),mycolor,'MarkerFaceColor',mycolor(1)),end
    end    
 end

results.XYZci = XYZ;

results.ax = ax;
results.fcrit = Fc;
results.fcrit_descrip = 'Critical F-value for test against origin; needs df check. Not used for P-values';

results.ci = ci;
results.ub = ub;
results.lb = lb;
results.eig = e;

return


