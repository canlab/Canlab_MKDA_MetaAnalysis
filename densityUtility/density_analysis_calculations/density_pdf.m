function [u,OUT] = density_pdf(n,iterations,radius_mm,varargin)
%function [u,OUT] = density_pdf(n,iterations,radius_mm,[output filename],[mask_file],[opt] RESUME)
%
% given a brain mask and n coordinate points, this function computes
% a null hypothesis probability density function based on the number 
% of iterations specified, and returns significant density thresholds
% at p < crit_p (default = 0.05).
%
% height thresholds
%
% n is the number of random points
% iteration is number of points in monte carlo distribution to sample
% radius_mm is the radius of the smoothing kernel in mm
% output filename is self-explanatory	[optional]
% mask_file is the name of an Analyze .img file containing 1's and 0's
%	1's denote the allowable region for point placement in the Monte Carlo
%	default is 'brain_avg152T1.img'
% RESUME: optional argument is file containing saved variables
% from an analysis that quit.
%
% considerations
% - edge effects in determining density
% - speed
% - multiple coordinates at the same point
%
% 
%
% Tor Wager, 2/6/2002

P = ['brain_mask_1mm.img'];
P = ['brain_avg152T1.img'];     % 2 mm voxels
P = 'scalped_avg152T1_graymatter_smoothed.img';

% 	radius_mm = 10;	% now an input argument
th = [.001 .005 .01];
cls = [10 100 1000];
crit_p = .05;
starti = 1;

if length(varargin) > 0
	fname = varargin{1};
else
	fname = input('Enter file name for the output ','s');
end

if length(varargin) > 1
	mask_file = varargin{2};
else
	mask_file = P;
end

t1 = clock;


% -----------------------------------------------------
% * process resume, if entered
% -----------------------------------------------------
if length(varargin) > 2
    fprintf(1,'Resuming...')
    if ~exist(varargin{1}) == 2, error('Resume file does not exist.'),end
    eval(['load ' varargin{1}])
    starti = length(maxd) + 1;
    fprintf(1,'Starting at iteration %3.0f of %3.0f\n ',starti,iterations)
end

for i = starti:iterations
    
    % -----------------------------------------------------
    % * get density mask of random values
    % -----------------------------------------------------
   
    if i == 1,
        [mask,searchmask,voxsize] = get_random_mask(n,P); 
        t3 = clock;, 
        oldrad = radius_mm;
        radius = radius_mm ./ mean(voxsize);
        sphere_vol = 4 * pi * radius_mm ^ 3 / 3;
        P = searchmask;     % stores the actual array in P  
    else
        [mask] = get_random_mask(n,P);
    end
                 
    dm = mask2density(mask,radius,searchmask,sphere_vol);
    
    
    % -----------------------------------------------------
    % * checking and setup stuff for 1st iteration
    % -----------------------------------------------------
    
    if i == 4   
        t4 = clock;
        if any(size(mask) - size(P)), error('Mask sizes do not match!'), end
        fprintf(1,'density.pdf | %3.0f voxels, radius = %3.0f mm = %3.2f voxels',sum(searchmask(:)>0),oldrad,radius);
        fprintf(1,'\t%3.0f iterations x %3.0f s per iteration\n',iterations,etime(t4,t3)./4)  
        clear oldrad, clear t3, clear t4
        figure;imagesc(dm(:,:,round(size(dm,3)./2))); pause(5); close
    end   
    
    % -----------------------------------------------------
    % * max height
    % -----------------------------------------------------
    maxd(i) = max(max(max(dm)));
    
    % -----------------------------------------------------
    % * cluster size distribution
    % -----------------------------------------------------
    [numc(:,:,i), cl_size(i,:)] = density2clusters(dm,th,cls);
    
    % -----------------------------------------------------
    % * save if we need to
    % -----------------------------------------------------
    if mod(i,20) == 0,
        eval(['save ' fname ' mask_file radius th cls n i dm maxd numc cl_size'])
        fprintf(1,'%d . ',i)
    end
    
    if mod(i,500) == 0
        t2 = clock;
        fprintf(1,'%3.0f s\n',etime(t2,t1))
    end
    
end

% -----------------------------------------------------
% * get the 95% value for all variables
% -----------------------------------------------------
[u,xpdf,ypdf,xcdf,ycdf] = d_pdf(maxd,crit_p);

index = 0;
for j = 1:size(numc,2)
    for i = 1:size(numc,1)

        index = index + 1;
        vec = str2num(num2str(numc(i,j,:)));
        cl_u(index) = d_pdf(vec,crit_p);

    end
end
cl_u = reshape(cl_u,size(numc,1),size(numc,2));

% -----------------------------------------------------
% * cluster size thresholds
%   for every specified nominal height threshold,
%   get the maximum cluster size that exceeds the 95% chance value
%   this gives us the size at which a cluster may be
%   considered significant; 95% of the time, no cluster
%   in the map will exceed this size.
% -----------------------------------------------------
for i = 1:size(cl_size,2)
    sz_u(i) = d_pdf(cl_size(:,i),crit_p);
end


OUT.fname = fname;
OUT.voxsize = voxsize;

OUT.sz_u = sz_u;
OUT.maxd = maxd;

OUT.mask_file = mask_file;
OUT.u    = u;
OUT.cl_u = cl_u;
OUT.crit_p = crit_p;
OUT.radius = radius;
OUT.th  = th;
OUT.cls = cls;
OUT.n   = n;
OUT.xpdf = xpdf;
OUT.ypdf = ypdf;
OUT.xcdf = xcdf;
OUT.ycdf = ycdf;
OUT.maxd = maxd;
OUT.numc = numc;
OUT.cl_size = cl_size;
%OUT.legend = mylegend;

eval(['save ' fname ' u OUT'])

[f1,f2,OUT] =  plot_dens_results(OUT);

t2 = clock;
fprintf(1,'\ndone!\t%3.0f s, average %3.0f s per iteration\n',etime(t2,t1),etime(t2,t1)/iterations) 

return


