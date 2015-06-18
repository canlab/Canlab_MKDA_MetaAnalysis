function [u,OUT] = density_pdf(n,n2,radius_mm,iterations,varargin)
%function [u,OUT] = density_pdf(n,n2,radius_mm,iterations,[opt] RESUME)
%
% given a brain mask and n coordinate points, this function computes
% a null hypothesis probability density function based on the number 
% of iterations specified, and returns significant density thresholds
% at p < crit_p (default = 0.05).
%
%
% it also computes pairs of density maps, and takes the difference
% between them, making pdfs for the difference between two maps
% thresholds are returned for each subtraction.
%
% height thresholds
%
% n     number of random points for first mask
% n2    number of random points for second mask
%
% critical thresholds at alpha = crit_p are denoted with variable u
% u     max density height threshold for map 1 (using n)
% u2    max density height threshold for map 2 (using n2)
% du    max density height threshold for map 1 - map 2
% du2   max density height threshold for map 2 - map 1
%
% sz_u  max size thresholds at alpha = th (also sz_u2, sz_du, etc.)
% cl_u  max number thresholds at alpha = th and extent = cls
%       (also cl_u2, cl_du, etc.)
%
% numc  : matrix of number of clusters
%         rows index density thresholds
%         columns index cluster size thresholds
%         with (th x cls) elements
%
% RESUME: optional argument is file containing saved variables
% from an analysis that quit.
% 
%
% Tor Wager, 2/6/2002

P = ['brain_mask_1mm.img'];
P = ['brain_avg152T1.img'];     % 2 mm voxels
P = 'scalped_avg152T1_graymatter_smoothed.img';

% radius_mm = 15;  this is now an input argument
th = [.0005 .001 .005 .01];
cls = [10 100 1000];
crit_p = .05;
starti = 1;

mask_file = P;
t1 = clock;
%fname = input('Enter file name for the output ','s');
fname = varargin{1};
fname2 = varargin{2};

% -----------------------------------------------------
% * process resume, if entered
% -----------------------------------------------------
if length(varargin) > 3
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
   
    if i == starti,
        [mask,searchmask,voxsize] = get_random_mask(n,P); 
        mask2 = get_random_mask(n2,P);
	disp(['mask1: ' num2str(n) ' points entered, ' num2str(sum(sum(sum(mask)))) ' points in mask.'])
 
	disp(['mask2: ' num2str(n2) ' points entered, ' num2str(sum(sum(sum(mask2)))) ' points in mask.'])
        t3 = clock;, 
        oldrad = radius_mm;
        radius = radius_mm ./ mean(voxsize);
        sphere_vol = 4 * pi * radius_mm ^ 3 / 3;
        P = searchmask;     % stores the actual array in P  
    else
        mask = get_random_mask(n,P);
        mask2 = get_random_mask(n2,P); 
    end
                 
    dm = mask2density(mask,radius,searchmask,sphere_vol);
    clear mask
    dm2 = mask2density(mask2,radius,searchmask,sphere_vol);
    clear mask2


    


    % -----------------------------------------------------
    % * checking and setup stuff for 1st iteration
    % -----------------------------------------------------
    
    if i == starti + 3   
        t4 = clock;
        % if any(size(mask) - size(P)), error('Mask sizes do not match!'), end
        fprintf(1,'density.pdf | %3.0f voxels, radius = %3.0f mm = %3.2f voxels',sum(searchmask(:)>0),oldrad,radius);
        fprintf(1,'\t%3.0f iterations x %3.0f s per iteration\n',iterations,etime(t4,t3)./4)  
        clear oldrad, clear t3, clear t4
        figure;imagesc(dm(:,:,round(size(dm,3)./2))); pause(5); close
    end   
    
    % -----------------------------------------------------
    % * max height
    % -----------------------------------------------------
    maxd(i) = max(max(max(dm)));
    maxd2(i) = max(max(max(dm2)));
    
    if i == starti,
	figure; subplot(1,2,1), 
	hist(dm(:),25)
	subplot(1,2,2)
	hist(dm2(:),25)
	drawnow
	disp(['maxima are ' num2str(maxd(i)) ' and ' num2str(maxd2(i))])
    end

    % -----------------------------------------------------
    % * cluster size distribution
    % -----------------------------------------------------
    [numc(:,:,i), cl_size(i,:)] = density2clusters(dm,th,cls);
    [numc2(:,:,i), cl_size2(i,:)] = density2clusters(dm2,th,cls);
    
    % -----------------------------------------------------
    % * differences: max height and cluster size
    % -----------------------------------------------------  
    d1 = dm - dm2;
    d2 = dm2 - dm;
    clear dm, clear dm2
    
    maxdiff(i) = max(max(max(d1)));
    maxdiff2(i) = max(max(max(d2)));
    
    [numcdiff(:,:,i), cl_size_diff(i,:)] = density2clusters(d1,th,cls);
    [numcdiff2(:,:,i), cl_size_diff2(i,:)] = density2clusters(d2,th,cls);
    
    
    % -----------------------------------------------------
    % * save if we need to
    % -----------------------------------------------------
    if mod(i,20) == 0,
        eval(['save ' fname ' mask_file radius th cls n i maxd numc cl_size maxd2 numc2 cl_size2 maxdiff maxdiff2 numcdiff numcdiff2 cl_size_diff cl_size_diff2'])
        fprintf(1,'%d . ',i)
    end
    
    if mod(i,500) == 0
        t2 = clock;
        fprintf(1,'%3.0f s\n',etime(t2,t1))
    end
    
end

fprintf(1,' | compiling stats\n')

% -----------------------------------------------------
% * max height thresholds: get the 95% value for all v
% -----------------------------------------------------
u = d_pdf(maxd,crit_p);
u2 = d_pdf(maxd2,crit_p);
du = d_pdf(maxdiff,crit_p);
du2 = d_pdf(maxdiff2,crit_p);

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
    sz_u2(i) = d_pdf(cl_size2(:,i),crit_p);
    sz_du(i) = d_pdf(cl_size_diff(:,i),crit_p);
    sz_du2(i) = d_pdf(cl_size_diff2(:,i),crit_p);
end
    
% -----------------------------------------------------
% * number of clusters thresholds
%   for every combination of height threshold and extent threshold,
%   get the number of clusters that exceeds the 95% chance value
%   this gives us the expected number of clusters for omnibus test
%   at some nominal threshold th
% -----------------------------------------------------
index = 0;
for j = 1:size(numc,2)
    for i = 1:size(numc,1)

        index = index + 1;
        vec = str2num(num2str(numc(i,j,:)));
        cl_u(index) = d_pdf(vec,crit_p);
        
        vec = str2num(num2str(numc2(i,j,:)));
        cl_u2(index) = d_pdf(vec,crit_p);
        
        vec = str2num(num2str(numcdiff(i,j,:)));
        cl_du(index) = d_pdf(vec,crit_p);
        
        vec = str2num(num2str(numcdiff2(i,j,:)));
        cl_du2(index) = d_pdf(vec,crit_p);

    end
end
cl_u = reshape(cl_u,size(numc,1),size(numc,2));
cl_u2 = reshape(cl_u2,size(numc2,1),size(numc2,2));
cl_du = reshape(cl_du,size(numcdiff,1),size(numcdiff,2));
cl_du2 = reshape(cl_du2,size(numcdiff2,1),size(numcdiff2,2));

OUT.fname = fname;
OUT.mask_file = mask_file;
OUT.voxsize = voxsize;

OUT.u    = u;
OUT.sz_u = sz_u;
OUT.cl_u = cl_u;
OUT.maxd = maxd;
OUT.cl_size = cl_size;
OUT.numc = numc;

OUT.u2    = u2;
OUT.sz_u2 = sz_u2;
OUT.cl_u2 = cl_u2;
OUT.maxd2 = maxd2;
OUT.cl_size2 = cl_size2;
OUT.numc2 = numc2;

OUT.du    = du;
OUT.sz_du = sz_du;
OUT.cl_du = cl_du;
OUT.maxdiff = maxdiff;
OUT.cl_size_diff = cl_size_diff;
OUT.numcdiff = numcdiff;

OUT.du2    = du2;
OUT.sz_du2 = sz_du2;
OUT.cl_du2 = cl_du2;
OUT.maxdiff2 = maxdiff2;
OUT.cl_size_diff2 = cl_size_diff2;
OUT.numcdiff2 = numcdiff2;

OUT.crit_p = crit_p;
OUT.radius = radius;
OUT.th  = th;
OUT.cls = cls;
OUT.n   = n;
OUT.n2  = n2;


eval(['save ' fname ' u OUT'])

density_results_table(OUT);
%[f1,f2,OUT] =  plot_dens_results(OUT);

t2 = clock;
fprintf(1,'\ndone!\t%3.0f s, average %3.0f s per iteration\n',etime(t2,t1),etime(t2,t1)/iterations) 

return


