function clnew = reclassify_clusters(cl,varargin);
% clnew = reclassify_clusters(cl,varargin);
%
% Takes a set of database peaks (XYZall) within a clusters structure
% and performs a clustering algorithm on those peaks
% Then creates a new clusters structure based on the classifications of
% peaks to clusters.
% See also xyz2cluster_class.m
% 
% requires XYZall field for all [x y z] points
% meant to use with database2clusters


% input args and defaults

k = 8;  % test case

if length(varargin) > 0, k = varargin{1};,end



% setup

N = fieldnames(cl(1));
XYZall = cat(1,cl(:).XYZall);

NN = {'title' 'threshold' 'voxSize' 'M' 'name' 'numVox' 'numpeaks' 'center' 'mm_center' 'P' 'imP' 'Z' 'XYZmm' 'XYZ'};
%[N] = reorder('delete',NN,N,N);

for i = 1:length(NN), wh = find(strcmp(N,NN{i}));, N(wh) = [];, end

c = clusterdata(XYZall,'maxclust',k,'linkage','average');

% create confidence volume for each, and use that as the new clusters
    
    
            
for m = 1:length(cl)
        
   c2 = c(1:length(cl(m).Study)); 
   c(1:length(cl(m).Study)) = [];
   
   
   
   for i = 1:max(c2)
        
       clnew(i).title = ['class' num2str(i)];

        wh = find(c2 == i);

        for j = 1:length(N)
        
            eval(['tmp = ''' N{j} ''';'])
            if ~isfield(clnew(i),tmp)
                eval(['clnew(i).' N{j} ' = [];'])
            end       
                
            str = ['clnew(i).' N{j} ' = [clnew(i).' N{j} '; cl(m).' N{j} '(wh,:)];'];
            eval(str);
            
        end
        
    end
    
end
        