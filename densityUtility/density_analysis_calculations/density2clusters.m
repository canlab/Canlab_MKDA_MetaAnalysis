function [numc, max_clsize] = density2clusters(dm,th,cls)
% 
% dm    : density mask
% th    : height threshold, or vector of thresholds
% cls   : vector of cluster sizes to test for
%
% numc  : matrix of number of clusters
%         rows index density thresholds
%         columns index cluster size thresholds

for i = 1:length(th)
    
    [maskedImage maskingImage] = maskImg(dm,th(i),Inf);
    [x,y,z] = ind2sub(size(maskingImage),find(maskingImage));     % find coords of values > 0 in mask
    XYZ = [x y z]';
    
    % oh, so slow
    %XYZ = (mask2voxel(maskingImage))';
    
    if size(XYZ,2) > 0
        
        % -----------------------------------------------------
        % * get list of all clusters and their sizes
        % -----------------------------------------------------
        cl = spm_clusters(XYZ);
    
        num(i) = max(cl);
    
        clear clsize
        for j = 1:num(i)
            clsize(j) = sum(cl == j);
        end
    
        max_clsize(i) = max(clsize);
        
        % -----------------------------------------------------
        % * get number of clusters at each threshold level
        % -----------------------------------------------------
        for j = 1:length(cls)
            numc(i,j) = sum(clsize > cls(j));
        end
        
    else
        numc(i,1:length(cls)) = 0;
        max_clsize(i) = 0;
    end
    
end