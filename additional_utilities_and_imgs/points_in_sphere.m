function [xyz2,wh,d] = points_in_sphere(XYZ,r,varargin)
% function [xyz2,wh,d] = points_in_sphere(XYZ,r,[xyz2], ['any'])
%
% returns coordinates within a radius of r units from center [x y z]
% if optional [xyz2] coords are listed, returns only those coords 
% also within that set.
%
% wh = output containing indices of which coords are within sphere.
% d = Euclidean distances
% see sphere_mask for function with mm to voxel coordinates
%
% Tor Wager
% Updated Oct 2008: Find coordinates within r mm of ANY input coords
%
% Optional input: 'any': find coords within r mm of ANY input coords
% Enter coordinates in n x 3 list.
%
% Example:
% [xyz3, wh2, d2] = points_in_sphere(rdacc(1).XYZmm', 10, [], 'any');

%XYZ = mm2voxel(XYZmm,V);
%rv = abs(r ./ diag(V.M(1:3,1:3))');   % voxel coords

% ANY mode: Recursively run this function for each coord
% Take intersection
% ------------------------------------------------------
if length(varargin) > 1 && strcmp(varargin{2}, 'any')

    varargin(2) = [];  % get rid of 'any' to avoid infinite loop!
    verbose = 0;

    ncoords = size(XYZ, 1);
    xyz2 = cell(ncoords, 1);
    wh = cell(ncoords, 1);
    d = cell(ncoords, 1);

    if verbose
        fprintf('ANY mode: finding coordinates near any input coordinate\n');
        fprintf('Running spheres for %3.0f coordinates\n', ncoords)
    end

    for i = 1:ncoords

        if verbose, fprintf('%3.0f ', i); end

        [xyz2{i}, wh{i}, d{i}] = points_in_sphere(XYZ(i, :), r, varargin{:});

        xyz2{i}(any(isnan(xyz2{i}), 2), :) = [];
        
    end

    if verbose, fprintf('\n'), end

    xyz2 = unique(cat(1, xyz2{:}), 'rows');
    
    [wh, indx] = unique(cat(1, wh{:}));

    d = cat(1, d{:});
    d = d(indx);
    return
end


% ------------------------------------------------------

build = 0;

if length(varargin) > 0 && ~isempty(varargin{1})
    xyz2 = varargin{1};
else
    build = 1;
end


if build
    % Build list of all coordinates

    % build box (list of XYZ coords) to improve speed

    lim = round([XYZ - r; XYZ + r]);
    diffs = diff(lim);

    xtmp = prod([diffs(2)+1 diffs(3)+1]);
    ztmp = prod([diffs(1)+1 diffs(2)+1]);

    x = repmat((lim(1,1):lim(2,1))',xtmp,1);

    y = []; for i=1:diffs(2)+1,
        ytmp = repmat(lim(1,2)+i-1,diffs(1)+1,1); y = [y;ytmp];,
    end
    y = repmat(y,diffs(3)+1,1);

    ztmp = repmat(lim(1,3),ztmp,1);
    z = [];
    for i = 1:diffs(3)+1
        z = [z; ztmp+i-1];
    end

    xyz2 = [x y z];

else

    % do not build = use existing coords

    %xyz2mm = voxel2mm(xyz2',V.mat)';
end



% build list of repeated sphere center
xyzc = repmat(XYZ,size(xyz2,1),1);

if isempty(xyzc), error('No coordinates in sphere.'), end

% distances from center
d = sum((xyz2 - xyzc) .^ 2,2) .^ .5;

wh = find(d <= r);
xyz2(d > r,:) = [];
d = d(wh);

return

