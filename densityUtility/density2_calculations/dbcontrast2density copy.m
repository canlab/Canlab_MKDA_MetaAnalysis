function [dmt,clusters,dm,OUT] = dbcontrast2density(DB,u,varargin)
% function [dmt,clusters,dm,OUT] = dbcontrast2density(DB,u,[radius_mm],[testfieldname],[contrast over levels],[con_dens_images])
% 
% XYZ is 3-column vector of coordinates
%
% Tor Wager 11/20/04
%
% step 1: read database
% step 2: database2clusters -- get database structure (don't use clusters; can enter any clusters here)
% 
%
% examples:
% [dmt,clusters,dm,OUT] = dbcontrast2density(EMDB,.05,10,'valence');
%
% load meta_analysis_master_file
% dbcontrast2density(PAINDB3,.05,15,'Right_vs_Left');
% dbcontrast2density(PAINDB3,.05,15,'Right_vs_Left',[1 1 0]);  % for r vs l
%

P = ['brain_avg152T1.img'];     % 2 mm voxels
P = 'scalped_avg152T1_graymatter_smoothed.img';

if length(varargin) > 0, radius_mm = varargin{1}; 
else, radius_mm = 10;
end
disp(['Radius is ' num2str(radius_mm) ' mm.'])


%if(length(fname) <= 4), error('Enter longer file name or type filename.img'),end
%if ~strcmp(fname(end-4:end), '.img'), [d fname] = fileparts(fname);, fname = [fname '.img'], end

mask_file = P;
t1 = clock;
fprintf(1,'Setup. ')

% -----------------------------------------------------
% * load standard brain
% -----------------------------------------------------
P = which(P);
V = spm_vol(P);
mask = zeros(V.dim(1:3));

voxsize = diag(V(1).mat)';
voxsize = voxsize(1:3);
radius = radius_mm ./ mean(voxsize);
sphere_vol = 4 * pi * radius_mm ^ 3 / 3;



% -----------------------------------------------------
% identify independent contrasts
% -----------------------------------------------------
testfield = 'valence';
if length(varargin) > 1, testfield = varargin{2};,end

% weight by sample size; field is Subjects
[allconditions,pointcounts,studynames,allcondindic,allcondnames] = dbcluster2indic(DB,DB,{testfield},1);

xyz = [DB.x DB.y DB.z];
eval(['allcond = DB.' testfield ';']); 
try, study = DB.Study;, catch, study = DB.study;, end

fprintf(1,'%3.0f Independent contrasts.',length(studynames))

OUT.allconditions = allconditions; 
OUT.pointcounts = pointcounts;
OUT.studynames = studynames;
OUT.allcondindic = allcondindic;
OUT.allcondnames = allcondnames;

% -----------------------------------------------------
% for each contrast, make a density image
% -----------------------------------------------------


PP = [];     % image names
if length(varargin) > 3, PP = varargin{4}; ,end            % load input filenames here!!
            
if isempty(PP)
    
    % normalize density maps so that 1 activation = value of 1 in map
    tmp = mask; tmp2 = round(size(mask)./2);
    tmp(tmp2(1),tmp2(2),tmp2(3)) = 1;
    spm_smooth(tmp,tmp,radius); % in vox if not a mapped vol! 
    normby = max(tmp(:));


    % Try to find images first, and if they don't exist in current dir,
    % create.
    fprintf(1,'Creating images.')
    
    for i = 1:length(studynames)
    
        wh = find(strcmp(study,studynames{i}) & strcmp(allcond,allconditions{i}));
        XYZmm = xyz(wh,:);
        XYZmm(any(isnan(XYZmm),2),:) = [];
    
        str = [studynames{i} '_' testfield '_' allconditions{i} '.img'];
        if exist(str) == 2, disp(['Found existing: ' str])
            
        else
            % doesn't exist yet
        if isempty(XYZmm)
            conmask = mask;
        else
            XYZvox = mm2vox(XYZmm,V.mat);   % get voxel coordinates
    
            conmask = xyz2mask(mask,XYZvox);   % put points in mask - could weight by Z-scores here, if desired.
            spm_smooth(conmask,conmask,radius)   %_mm); % smooth it! radius in voxels if no header info is given.
            
            conmask = conmask ./ normby;    % normalize so center of study activation = 1
            conmask(conmask > 1) = 1;       % limit to max of 1, in case multiple nearby points in same contrast
                                            % max activation for a single
                                            % contrast is 1.
    
        end
        
        
        V.fname = str;
        warning off, spm_write_vol(V,conmask);, warning on
        
        end % if find existing
        PP = str2mat(PP,str);
    end % loop thru contrasts
    PP = PP(2:end,:);
else
    % we have filenames already created
end

fprintf(1,'Done %3.0f images in %3.0f s',length(studynames),etime(clock,t1))
    
OUT.PP = PP;

% -----------------------------------------------------
% across all density images, compute prob(activation), p(a)
% -----------------------------------------------------
tor_spm_mean_ui(PP,'pa_overall.img');

% -----------------------------------------------------
% compute prob(task) for each allcondition
% -----------------------------------------------------



%%% THRESHOLD THE OVERALL P(ACTIVATION) IMAGE
% ----------------------------------------------------
% * get list of all coordinates; 
% ----------------------------------------------------
    for i = 1:length(studynames)
    
        wh = find(strcmp(study,studynames{i}) & strcmp(allcond,allconditions{i}));
        XYZmm = xyz(wh,:);
        XYZmm(any(isnan(XYZmm),2),:) = [];
        
        XYZmm_all{i} = XYZmm;
    end
        
    
    % ----------------------------------------------------
% * get list of all coordinates; select n at random
%   vol: row, col, array = x, y, z
% ----------------------------------------------------
vol = spm_read_vols(V);
[x,y,z] = ind2sub(size(vol),find(vol));     % find values > 0 in vol
allXYZ = [x y z];                           % XYZ is in 3 columns in this function

iterations = 5000;
testp = [32 24 32]; % ignores edge effects!

% make convolution kernel -- weighting function for distance
s = radius/sqrt(8*log(2));    % kernel st. dev in VOXEL distance
m = ceil(4*s);                      % limit of kernel -- 4 std. deviations
                              % kernel for positive distances
y = normpdf([0:m],0,s); y = y ./ max(y);

for i = 1:length(XYZmm_all), if isempty(XYZmm_all{i}), whom(i)=1;,else, whom(i)=0;,end,end
XYZmm_all(find(whom)) = [];

disp('Iterating null hypothesis for overall density.')
for it = 1:iterations
    
    for i = 1:length(XYZmm_all)
        
        % select random coordinates for each study
        % n is coords for this study, whichv indexes n random in-brain-mask
        % coords
        n = size(XYZmm_all{i},1);
        whichv = ceil(rand(1,n) * size(allXYZ,1));
        XYZ = allXYZ(whichv,:);

        % euclidean distance from test point, in voxel space (coords) 
        d = distance(testp,XYZ);
        
        % weight distances by gaussian function
        % to get sum of weighted distances from fixed test point
        d(d>m-1) = [];    % more than 4 stds
        
        dstudy(i) = sum(y(round(d)+1));   % this is the weighted sum, the statistic value of interest
        
    end
    
    dstudy(dstudy > 1) == 1;        % cap at 1, because no study can exceed 1 
                                    % for z-score weighting, cap would be
                                    % z(study)
    pa_h0(it) = mean(dstudy);       % weighted p(a) overall, by chance.  save this distribution.
                                    
                                    
end
       


% threshold the overall density image
cl = dens_fdr_thresh('pa_overall.img',pa_h0,1);

save pa_overall_cl cl




% NOW WORK ON THE CHI2 VALUES FOR DIFFS AMONG CONDITIONS
% COMPUTE CHI2 VALUE AT EACH VOXEL
% WE NEED ROWS = VOXELS, COLUMNS = P(ACT|TASK X), OR % STUDIES IN TASK X

% -----------------------------------------------------
% across images for each cond., compute prob(activation) given condition,
% p(a)|task
% -----------------------------------------------------

% select conditions to test across
if length(varargin) > 2, contrast = varargin{3}; ,
else
    contrast = ones(length(OUT.allcondindic));
end
    

allcondindic = OUT.allcondindic(:,find(contrast));
allcondnames = OUT.allcondnames(find(contrast));


% number of contrasts in each map
numcons = sum(allcondindic);
    
for i = 1:size(allcondindic,2)  % for each level
    names = OUT.PP(find(allcondindic(:,i)),:);
    nstr = ['pa_given_' allcondnames{i} '.img'];
    tor_spm_mean_ui(names,nstr);
    
    % convert probability scores to count for chi2 test
    V = spm_vol(nstr); v = spm_read_vols(V);
    v = v .* numcons(i);    % counts!
    ncstr = ['counts_given_' allcondnames{i} '.img'];
    V.fname = ncstr;
    spm_write_vol(V,v);
    
    if i == 1,
        OUT.cntimgnames = ncstr;
    else
        OUT.cntimgnames = str2mat(OUT.cntimgnames,ncstr);
    end
end

% load images and make voxels x conditions matrix countdat
V = spm_vol(OUT.cntimgnames); v = spm_read_vols(V);

Vmask = spm_vol(which(mask_file)); vm = spm_read_vols(Vmask);
wh = find(vm>0);    % wh contains indices in mask, for use later to convert back to voxels!
    
clear countdat
for i = 1:size(v,4) % for each condition
    vtmp = v(:,:,:,i);
    countdat(:,i) = vtmp(wh);
end

% eliminate voxels w/ zeros in all conds
whzero = find(all(countdat<eps,2));
countdat(whzero,:) = [];
wh(whzero) = [];

% eliminate voxels w/ no more than 5 counts in any condition
mymax = max(countdat')';
whlow = find(mymax < 5);
countdat(whlow,:) = [];
wh(whlow) = [];

% compute chi2 values for all eligible voxels
disp('Computing chi2 values')
chi2 = computeynchi2(countdat,numcons);

% write p-value mask
pmap = ones(size(v(:,:,:,1)));
pmap(wh) = chi2(:,3);   % p-values
pmap(vm==0 | isnan(vm)) = NaN;
pmap = double(pmap);
nstr = []; for i =1:length(allcondnames),nstr = [nstr '_' allcondnames{i}];,end
pname = ['chi2_pmap' nstr];
V = spm_vol(which(mask_file));
V.fname = pname;
spm_write_vol(V,pmap);

% write chi2 mask
chimap = pmap;
chimap(wh) = chi2(:,1);   % chi2-values
pname = ['chi2_map_' nstr];
V.fname = pname;
spm_write_vol(V,chimap);

% FDR threshold -- use all in-mask voxels as search area
tmp = pmap; tmp(vm==0 | isnan(vm))=[];   % eliminate out-of-analysis voxels
pt = FDR(tmp,.05);  clear tmp;
if ~isempty(pt)
    whfdr = logical(pmap <= pt);
else
    disp('No FDR corrected chi2 results.')
    whfdr = [];
end
warning off
fdrmap = zeros(V.dim(1:3));
fdrmap(whfdr) = 1;
pname = ['chi2_sig_fdr_' nstr];
V.fname = pname;
spm_write_vol(V,fdrmap);
warning on

% .05 threshold
sigmap = zeros(V.dim(1:3));
sigmap(wh) = chi2(:,5); % has 1's and 0's for sig already
%whsig = wh(find(chi2(:,4)));
%sigmap(whsig) = 1;
pname = ['chi2_sig_05_' nstr];
V.fname = pname;
warning off
spm_write_vol(V,sigmap);
warning on
if sum(sigmap) < eps, 
    disp('No chi2 results at p < .05');,
else
    disp([num2str(sum(sigmap(:))) ' Voxels at p < .05.'])
    cl = mask2clusters([pname '.img']);
    save chi2_p05_cl cl
end

save OUT OUT

return






% -----------------------------------------------------
% * sub-functions
% -----------------------------------------------------

function XYZout = mm2vox(XYZ,M)
% converts a list of coordinates
% XYZ should be 3-column vector [x y z]
% calls tal2vox

XYZ = XYZ';
for i = 1:size(XYZ,2), 
	XYZout(i,:) = tal2vox(XYZ(:,i),M); 
end

XYZout = round(XYZout);

return



% -----------------------------------------------------
function vox=tal2vox(tal,M)
% converts from talairach coordinate to voxel coordinate
% based on variables from SPM.M (passed here for 
% faster operation)
% e.g., foo=tal2vox([-30 28 -30], VOL)
% from Russ Poldrack's spm ROI utility 

vox=[0 0 0];
vox(1)=(tal(1)-M(1,4))/M(1,1);
vox(2)=(tal(2)-M(2,4))/M(2,2);
vox(3)=(tal(3)-M(3,4))/M(3,3);

return




function mask = xyz2mask(mask,XYZ)

% -----------------------------------------------------
% * make a mask out of XYZ, in space of P
% * deal with repeated coordinates by adding
% -----------------------------------------------------

ind = sub2ind(size(mask),XYZ(:,1),XYZ(:,2),XYZ(:,3));
mask(ind) = 1;

ind = sort(ind);
repeats = ind(find(~diff(ind)));            % index values that are repeated
for i = 1:length(repeats)
    mask(repeats(i)) = mask(repeats(i)) + 1;
end
