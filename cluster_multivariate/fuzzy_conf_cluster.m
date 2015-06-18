function OUT = fuzzy_conf_cluster(im,xyz)
% OUT = fuzzy_conf_cluster(im,xyz)
%
% tor wager
% 
% input:
% this function takes a set of indicator vectors (im)
% coded as 1's and 0's, and a list of varables (xyz)
% 
% output: 
% a permutation test for whether there are separate regions 
% of the space defined by the columns of xyz
% occupied by objects in the classes specified in the columns 
% of im
% 
% this can be thought of as a clustering algorithm,
% based on something similar to the silhouette index of 
% Kaufmann and Rousseeuw (1990)
% with one important difference:
% the algorithm confirms whether classes occupy different regions
% of space, rather than trying to isolate the space they occupy
% (as do data-driven clustering algorithms)
%
% think of it as a confirmatory version of the partitioning
% around medoids algorithm
%
% objects can have membership in more than one class simultaneously.
%
% can replace and probably be better than linear discriminant analysis
% because it makes fewer assumptions.

niter = 5000;

% ---------------------------------------------
% correct permutation
% ---------------------------------------------

[OUT.avgq, OUT.q,center] = doquality(im,xyz);

% ---------------------------------------------
% bootstrap
% ---------------------------------------------

t1 = clock;
for i = 1:niter

    % permute columns of xyz
    for j = 1:size(xyz,2)
        xyz2(:,j) = getRandom(xyz(:,j));
    end

    [OUT.ho_avgq(i), OUT.ho_q(i,:)] = doquality(im,xyz2);
    
    if mod(i,500)==0, fprintf(1,'.'), end
end

fprintf(1,' done (in %4.0f s)!\n',etime(clock,t1))

% ---------------------------------------------
% table of results
% ---------------------------------------------

fprintf(1,'\nPermutation test on distribution of classes over xyz space\n')
fprintf(1,'Omnibus: \tq = \t%3.3f\t, expected q = \t%3.3f\t, p = \t%3.4f\t\n', ...
    OUT.avgq, mean(OUT.ho_avgq), sum(OUT.avgq <= OUT.ho_avgq) ./ length(OUT.ho_avgq))

fprintf(1,'\t95%% sig. level is at q = %3.3f\n',prctile(OUT.ho_avgq,95))

fprintf('\nIndividual classes\n')
fprintf(1,'x\ty\tz\tq\texp. q\tthreshold\tp\t\n')

for i = 1:size(OUT.ho_q,2)
    tmp = OUT.ho_q(:,i);
    
    if ~any(isnan(tmp))
        fprintf(1,'%3.0f\t%3.0f\t%3.0f\t%3.3f\t%3.3f\t%3.3f\t%3.4f\t\n', ...
        center(i,1), center(i,2), center(i,3), ...
        OUT.q(i), ...
        mean(tmp), ...
        prctile(tmp,95), ...
        sum(OUT.q(i) <= tmp) ./ length(tmp))
    else
        fprintf(1,'Insufficient objects for this class\n')
    end
    
end

return








function [avgq,q,center] = doquality(im,xyz)

for i = 1:size(im,2)    % for each class
    tmp = mean(xyz(find(im(:,i)),:),1);     % get center
    center(i,:) = tmp;
    
    % dist of all points from each center
    d(:,i) = sum((repmat(tmp,size(xyz,1),1) - xyz).^2,2).^0.5;
    
end

for i = 1:size(im,2)    % for each class
    
    tmp = mean(d(find(im(:,i)),:),1);   % get mean dist to all centers
    a = tmp(i);                         % distance to same class center
    tmp(i) = [];
    b = min(tmp);                       % distance to nearest neighboring class
    q(i) = (b - a) ./ max(a,b);            % quality index for this class
    
end
    
avgq = nanmean(q);

return





