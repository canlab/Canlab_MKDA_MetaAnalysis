function [p0,pa_h0] = p0_monte_carlo_FWE(DB,iterations)
% [p0,pa_h0] = p0_monte_carlo_FWE(DB,iterations)
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

XYZmm_all = DB.XYZmm;

for i = 1:length(XYZmm_all), if isempty(XYZmm_all{i}), whom(i)=1;,else, whom(i)=0;,end,end
XYZmm_all(find(whom)) = [];

fprintf(1,'Iterating null hypothesis for overall density: %3.0f iterations\n')

fprintf(1,'Iteration: %04d',0)
tic

conmask = [];
    
for it = 1:iterations
    
    
    for i = 1:length(XYZmm_all)
        
        % select random-location coordinates for each study
        % n is coords for this study, whichv indexes n random in-brain-mask
        % coords
        n = size(XYZmm_all{i},1);
        whichv = ceil(rand(1,n) * size(allXYZ,1));
        XYZ = allXYZ(whichv,:);

        conmask(:,:,:,i) = xyz2density(XYZ,DB.maskzeros,[],[],DB.radius,DB.studyweight(i),1);
        
        dstudy = sum(conmask,4);   % this is the weighted sum, the statistic value of interest
        
    end
    
    pa_h0(it) = max(dstudy(:));  % weighted max p(a) overall, by chance.  save this distribution.
                  
    fprintf(1,'\b\b\b\b%04d',it)
                                    
end
     
fprintf(1,' Done.\n')
toc

p0 = prctile(pa_h0,95);


return

