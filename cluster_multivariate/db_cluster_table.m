function db_cluster_table(clusters,fnames)
% function db_cluster_table(clusters,fnames)
%
% tor wager
% counts studies and contrasts in each cluster
%
% clusters is a struct in which XYZmm field has list of points
% --should be output of database2clusters, which contains all fields
%
% Prints one row per coordinate in any cluster, along with study names,
% cluster ID, and values in fields you specify in fnames
%
% fnames is a cell array that has names of all fields in clusters to use
% example:
% db_cluster_table(PAIN,{'Right_vs_Left'})
%
% see also dbcluster_contrast_table.m - a very useful function!

fprintf(1,'Study\tCluster\tx\ty\tz\t')
for i = 1:length(fnames)
    fprintf(1,'%s\t',fnames{i})
end
fprintf(1,'\n')

for i = 1:length(clusters)
    
    cl = clusters(i);
    
    for k = 1:length(clusters(i).x)
    
        fprintf(1,['%s\t%3.0f\t%3.0f\t%3.0f\t%3.0f\t'],cl.Study{k},i,cl.x(k),cl.y(k),cl.z(k))
        
        for j = 1:length(fnames)

            eval(['tmp1 = cl.' fnames{j} '(k);'])
            if iscell(tmp1)
                fprintf(1,'%s\t',tmp1{1})
            else
                fprintf(1,'%3.0f\t',tmp1)
            end
            
        end
        
        fprintf(1,'\n')
        
    end
    
end
         