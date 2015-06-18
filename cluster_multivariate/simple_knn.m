function [clout,testclass] = simple_knn(data,clas,k,varargin)
% function [clout,testclass] = simple_knn(data,clas,k,[test set])
% used in tor_knn
% tor wager

% prior class probabilities - percentage of points in each class
for i = 1:max(clas), pcl(i) = sum(clas == i);, end
pcl = pcl ./ sum(pcl);

% weights for each class
pcl = 1 - pcl;


  
  if length(varargin) > 0
      origsize = size(data,1);
      data = [data; varargin{1}];
  end
  
  
[Ntst,nent]=size(data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of the K-nearest neighbours in the training set

  dst=[];ist=[];
  
  for i = 1:Ntst,
      
	dist=sum( ((ones(Ntst,1)*data(i,:)) - data)' .^2 )';
    if length(varargin) > 0, dist = dist(1:origsize);,end
	[dss,iss]=sortrows([clas dist],2);
   
    
    % trouble is that with categorical data, many values may be equally
    % distant.  to choose among them, we don't want to take just the first
    % k; we want to assign a class when distances are equivalent based on
    % the overall proportions of class memberships.  thus, if 5 points are
    % 3 units away from the target point, and 2 are class 1 and 3 are class
    % 2, then we want to choose based on the mean of those: mean[1 1 2 2
    % 2], weighted by the overall probability, say 60% of points in class 1
    
    % this part is new, invented by me...so decide for yourself...
   
   tmp = dss(1:k,:);
   % at the last point, we may have arbitrarily selected among many points
   % of same distance, so check:
   wh = tmp(:,2) == tmp(end,2);     % which neighbors share same dist as end
   w = clas(dist == tmp(end,2));    % which points are largest included distance
   
   if length(w) > sum(wh)   % if more at this distance than included
       
        clear we
        for j = 1:length(w), we(j) = pcl(w(j));, end
        we = we ./ sum(we);     % we is weights; normalize
        tmp(wh,1) = we * w;     % fuzzy class; weighted avg.
    end
   
    clout(i,1) = mean(tmp(:,1));  
    
  end
  
  clout = round(clout); 
  
  if length(varargin) > 0
      testclass = clout(origsize+1:end);
      clout = clout(1:origsize);
  end

  
  return
  