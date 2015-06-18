function [u,maxd,usum,maxsum] = fast_max_density(n,iterations,radius_mm,varargin)
% [u,maxd,usum,maxsum] = fast_max_density(n,iterations,radius_mm,[brain_mask_name])
%
% Improved faster algorithm for null-hypothesis monte carlo simulation
% Same results as density_pdf.m
%
% u = critical threshold at alpha = .05
% maxd = distribution of maximum density values under Ho
%
% tor wager
% Feb 4,2003
%
% t1 = clock; [u,maxd]=fast_max_density(10,100,20); t2 = etime(clock,t1) \
% takes about 80 s on my 1.3 PIII, or 7:20 for n = 100
% estimated 6 hrs 20 mins for 5000 iterations at n = 100

if length(varargin) > 0
    P = varargin{1};
else
    P = ['brain_avg152T1.img'];     % 2 mm voxels
    P = 'scalped_avg152T1_graymatter_smoothed.img';
end

if length(n) > 1
    % we have z-scores, so use them
    zscores = n;
    n = length(zscores);
end
    
% ----------------------------------------------------
% * read the input brain mask image
% ----------------------------------------------------
if isempty(P)
    P = spm_get(1,'*.img','Select brain or gray matter mask');
end
if ischar(P)
    P = which(P);
    V = spm_vol(P);
    vol = spm_read_vols(V);
    voxsize = diag(V(1).mat)';
    voxsize = voxsize(1:3);
else
    vol = P;
    voxsize = NaN;
end

vox_radius = radius_mm ./ mean(voxsize);
sphere_vol = 4 * pi * radius_mm ^ 3 / 3;
diam = 2 * vox_radius;   % diameter of sphere in voxels

% ----------------------------------------------------
% * get list of all coordinates in brain mask
%   set up stuff for computing density 
%   as in mask2density.m, but faster
%   vol: row, col, array = x, y, z
% ----------------------------------------------------
[x,y,z] = ind2sub(size(vol),find(vol));     % find values > 0 in vol
allXYZ = [x y z];                           % XYZ is in 3 columns in this function
Q       = ones(1,size(allXYZ,1));
slices  = allXYZ(:,3)';                      % used to pass only certain slices to subfunctions


% ----------------------------------------------------
% * set up the gui figure
% ----------------------------------------------------
f = figure('Color','w');
tmp = get(gcf,'Position') .* [1 1 .5 .1];
set(gcf,'Position',tmp)
set(gcf,'MenuBar','none','NumberTitle','off')
figure(f), set(gca,'Xlim',[0 100])



% ----------------------------------------------------
% * do the pdf
% ----------------------------------------------------

for i = 1:iterations

    % select n coordinates at random
    whichv = ceil(rand(1,n) * size(allXYZ,1));
    XYZ = allXYZ(whichv,:);
    
    % get distances to all units (neural network toolbox function)
    d = dist(XYZ');

    % threshold the distances at 2 * radius (within-sphere diameter)
    d2 = d; d2(d2 <= diam) = 1; d2(d2 > diam) = 0;

    if exist('zscores') == 1
        % weight each row by z-scores for each point, if entered
        zmat = repmat(zscores',size(d2,1),1);
        d2 = d2 .* zmat;
    end
    
    % find the point that has the maximal number of other within-sphere points
    s = sum(d2,2);  % sum across columns = z-score weighted sum if zscores entered
    wh = find(s == max(s));
    wh = wh(1);

    % restrict set of XYZ coordinates
    vr = XYZ(find(d2(wh,:)),:)';
    [vr,si] = sortrows(vr',3);
    vr = vr';
    
    if exist('zscores') == 1
        z2 = zscores(find(d2(wh,:)),:);
        z2 = z2(si);    % sort z-scores to match vr (restricted voxels)
    end
    
    % compute density
    upd = [1 diff(vr(3,:)) > 0];    % points where slice changes
    dm      = zeros(size(vol));
    
    for j = 1:size(vr,2)                                                % for each reported point...
    
        if upd(j)
            vin = allXYZ(slices <= vr(3,j) + vox_radius & slices > vr(3,j) - vox_radius,:)';    % these are coordinates that may be in-sphere
        end
    
        if exist('zscores') ~= 1
            % count of number of peaks at this exact voxel center
            count   =  sum(vr(1,:) - vr(1,j) == 0 & vr(2,:) - vr(2,j) == 0 & vr(3,:) - vr(3,j) == 0);
        else
            % we have z-scores
            % sum z-scores for each point at this exact voxel center
            yesses = (vr(1,:) - vr(1,j) == 0 & vr(2,:) - vr(2,j) == 0 & vr(3,:) - vr(3,j) == 0);
            count = sum(yesses .* z2'); 
        end
                                       
        dm      = add_density_sphere(vr(:,j),vin,Q,count,dm,vox_radius);        % add density to sphere around this point
    end
    
    mymax = max(dm(:));
    maxd(i) = mymax(1);

    
    if rem(i,10) == 0
        figure(f), try,barh(100*i / iterations),catch,end
        set(gca,'Xlim',[0 100]),set(gca,'YTickLabel',i),drawnow
    end
    
end
    
maxn = maxd;
maxsum = maxd;
maxd = maxd ./ sphere_vol;
u = prctile(maxd,95);
usum = prctile(maxsum,95);
    
close    
    
    
    
return




function [dm,j] = add_density_sphere(coord,vindex,Q,count,dm,radius)
    % j is index numbers of in-sphere coordinates
    % coord is coordinate of sphere center
    % vindex is XYZ index of all slices in the vicinity
    % Q is ones of (3,length vindex)
    % count is the number of counts at the sphere center
    
    j       = find(sum((vindex - coord*Q(:,1:size(vindex,2))).^2) <= radius^2);

    % turn index of restricted vindex space back to whole-mask coordinates
    XYZ = vindex(:,j);
    ind = sub2ind(size(dm),XYZ(1,:),XYZ(2,:),XYZ(3,:));
    
    % add the mask value at the point to voxels in sphere; mask may be > 1
    % do not divide by sphere_vol - instead, do this at the end.
    dm(ind) = dm(ind) + count;
    
    
return
