function [p0,pa_h0] = p0_monte_carlo(DB)
% [p0,pa_h0] = p0_monte_carlo(DB)
%
% needs fields:
% DB.XYZmm    % mm coords for each study
% DB.maskV    % mask mapped volume
% DB.radius   % in voxels


% ----------------------------------------------------
% * get list of all coordinates; select n at random
%   vol: row, col, array = x, y, z
% ----------------------------------------------------
vol = spm_read_vols(DB.maskV);              % read the mask image
[x,y,z] = ind2sub(size(vol),find(vol));     % find values > 0 in vol
allXYZ = [x y z];                           % XYZ is in 3 columns in this function

iterations = 15000;
testp = [32 24 32]; % ignores edge effects! pick a voxel in the middle of the brain to test

XYZmm_all = DB.XYZmm;

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
        
        dstudy(i) = any(d <= DB.radius);   % this is the weighted sum, the statistic value of interest
        
    end
    

    dstudy = dstudy * DB.studyweight;    % weight by weights: e.g., sqrt relative sample size 
    
    pa_h0(it) = mean(dstudy);       % weighted p(a) overall, by chance.  save this distribution.
                        
    if mod(it,round(iterations/10)) == 0, fprintf(1,'%3.0f . ',it),end
                                    
end
     
p0 = prctile(pa_h0,95);


return

