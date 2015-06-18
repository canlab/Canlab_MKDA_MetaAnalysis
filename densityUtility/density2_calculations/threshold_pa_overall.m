function cl = threshold_pa_overall(OUT,varargin)
% cl = threshold_pa_overall(OUT,[task level name])
%
% if no opt. argument, works with overall pa and pa_overall.img
% Otherwise, restricts analysis to specified task level only
%
% saves output struture OUT in STATS_VARS.mat
% saves cl structure to file
%
% called by dbcontrast2density
% tor wager
%
% example:
% to add a task-specific Ho distribution to OUT,
% load OUT
% cl = threshold_pa_overall(OUT,'pos')
%
% to get a new variable to test, and add an Ho distribution
% OUT = get_contrast_indicator(DB,'Method');     % method is the new var
% cl = threshold_pa_overall(OUT,'pos')

% variables we need

studynames = OUT.studynames;    % study names
allconditions = OUT.allconditions;  % conditions for each study(contrast)
allcond = OUT.allcond;          % task conditions for all peaks
study = OUT.study;              % study names for all peaks
xyz = OUT.xyz;                  % all coordinates
V = OUT.maskV;                  % memory-mapped mask volume
studyweight = OUT.studyweight;  % sqrt of relative sample size w/i task condition
radius = OUT.radius;            % FWHM for smoothing in voxels

% specify overall or specific level
filename = 'pa_overall.img'; level = [];
if length(varargin) > 0, level = varargin{1}; filename = ['pa|' levelname];, end

if ~(exist(filename) == 2)
    % we don't have overall image; let's create it
    % ***
end


if isfield(OUT,'pa_h0')
    
    % we already have a null H distribution -- skip the Monte Carlo
    disp('Ho distribution already saved ... skipping Monte Carlo')
    pa_h0 = OUT.pa_h0;
    
else
    
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
 
    % if specific level is entered, restrict to those studies w/correct
    % level value
    if ~isempty(level), 
        wh = find(strcmp(level,allconditions));
        XYZmm_all = XYZmm_all(wh);
    end
    
% ----------------------------------------------------
% * get list of all coordinates; select n at random
%   vol: row, col, array = x, y, z
% ----------------------------------------------------
vol = spm_read_vols(V);                     % read the mask image
[x,y,z] = ind2sub(size(vol),find(vol));     % find values > 0 in vol
allXYZ = [x y z];                           % XYZ is in 3 columns in this function

iterations = 5000;
testp = [32 24 32]; % ignores edge effects! pick a voxel in the middle of the brain to test

% make convolution kernel -- weighting function for distance
s = radius/sqrt(8*log(2));    % kernel st. dev in VOXEL distance
m = ceil(4*s);                      % limit of kernel -- 4 std. deviations
                              % kernel for positive distances
y = normpdf([0:m],0,s); y = y ./ max(y);

for i = 1:length(XYZmm_all), if isempty(XYZmm_all{i}), whom(i)=1;,else, whom(i)=0;,end,end
XYZmm_all(find(whom)) = [];

fprintf(1,'Iterating null hypothesis for overall density: %3.0f iterations\n')

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
                                    
    dstudy = dstudy * studyweight' ./ length(dstudy);    % weight by sqrt relative sample size 
    
    pa_h0(it) = mean(dstudy);       % weighted p(a) overall, by chance.  save this distribution.
                        
    if mod(it,round(iterations/10)) == 0, fprintf(1,'%3.0f . ',i),end
                                    
end
     
str = ['OUT.pa_h0' level ' = pa_h0;'];  % save ho dist with correct name appended
save STATS_VARS OUT

end % if we already have Ho distribution


% threshold the overall density image
cl = dens_fdr_thresh(filename,pa_h0,1);

str = ['save ' filename '_cl cl'];
eval(str)

if ~isempty(cl), 
    spm_image('init',filename);
    cluster_orthviews(cl,{[1 0 0]},'add');
end

return

