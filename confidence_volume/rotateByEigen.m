

% 
% function to rotate a 3d object to a new orientation defined
% by an eigenvalue problem.
%
% inputs are either:
%   X, Y, Z matrices defining points to rotate
%   n x 3 pointlist of [X Y Z] values, with 2nd and 3rd args empty
%
% returns: structure with X Y Z rotated, and XYZ rotated pointlist
%
% results = rotateByEigen(X,Y,Z,V)
% results = rotateByEigen(X,[],[],V)
%
% Robert Welsh, 11/11/01
% Modified by Tor Wager to take point list as well.

function results = rotateByEigen(X,Y,Z,V)

xP = 0*X;
yP = 0*Y;
zP = 0*Z;

if isempty(Y) & isempty(Z)
    xP = [];
    nX = 1;
    nY = size(X,1);
else
    [nX nY] = size(X);
end

for iX = 1:nX
    for iY = 1:nY
        if isempty(Y) & isempty(Z)
            xOld = X(iY,:)';
        else
            xOld = [X(iX,iY);Y(iX,iY);Z(iX,iY)];
        end
        
        xNew = V*xOld;  
        
        if isempty(Y) & isempty(Z)
            XYZ(iY,:) = xNew';
        else
            xP(iX,iY) = xNew(1);
            yP(iX,iY) = xNew(2);
            zP(iX,iY) = xNew(3);
        end
        
        
    end
end

if isempty(Y) & isempty(Z)
    results.XYZ = XYZ;
else
    results.xP = xP;
    results.yP = yP;
    results.zP = zP;
end

%
% all done
%

return
