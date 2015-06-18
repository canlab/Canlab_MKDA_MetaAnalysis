function overlap_table(OUT)
% function overlap_table(OUT)
% 
% prints table from output of cluster_overlap_npm.m
% tor wager


fprintf(1,'\nMonte Carlo test on spatial overlap of clusters\n')
fprintf(1,'Omnibus: \tobs = %3.0f voxels\t, exp = %3.2f voxels\t, p = %3.4f\t\n', ...
    OUT.obs_o1, mean(OUT.o1), sum(OUT.o1 >= OUT.obs_o1) ./ length(OUT.o1))
fprintf(1,'\t95%% sig. level is at %3.0f voxels of overlap\n',prctile(OUT.o1,95))

fprintf('\nIndividual clusters: first set\n')
fprintf(1,'Index\tx\ty\tz\tvoxels\tobs. overlap\tx\ty\tz\texp. overlap\tthreshold\tp\t')
    if isfield(OUT,'ba1'), fprintf('BA\t'),end
    fprintf(1,'\n')

for i = 1:size(OUT.ovec1,2)
    tmp = OUT.ovec1(:,i);
    fprintf(1,'%3.0f\t%3.0f\t%3.0f\t%3.0f\t%3.0f\t%3.0f\t%3.0f\t%3.0f\t%3.0f\t%3.2f\t%3.0f\t%3.4f\t', ...
    i, OUT.cl1(i).mm_center(1), OUT.cl1(i).mm_center(2), OUT.cl1(i).mm_center(3), ...
    size(OUT.cl1(i).XYZ,2), ...
    OUT.obs_ovec1(i), ...
    OUT.cent1(i,1),OUT.cent1(i,2),OUT.cent1(i,3), ...
    mean(tmp), ...
    prctile(tmp,95), ...
    sum(tmp >= OUT.obs_ovec1(i)) ./ length(tmp))

    if isfield(OUT,'ba1'), try,fprintf('%s\t',OUT.ba1{i}),catch,end,end
    fprintf(1,'\n')
end

fprintf('\nIndividual clusters: second set\n')
fprintf(1,'Index\tx\ty\tz\tvoxels\tobs. overlap\tx\ty\tz\texp. overlap\tthreshold\tp\t')
    if isfield(OUT,'ba1'), fprintf('BA\t'),end
    fprintf(1,'\n')
    
for i = 1:size(OUT.ovec2,2)
    tmp = OUT.ovec2(:,i);
    fprintf(1,'%3.0f\t%3.0f\t%3.0f\t%3.0f\t%3.0f\t%3.0f\t%3.0f\t%3.0f\t%3.0f\t%3.2f\t%3.0f\t%3.4f\t', ...
    i, OUT.cl2(i).mm_center(1), OUT.cl2(i).mm_center(2), OUT.cl2(i).mm_center(3), ...
    size(OUT.cl2(i).XYZ,2), ...
    OUT.obs_ovec2(i), ...
    OUT.cent2(i,1),OUT.cent2(i,2),OUT.cent2(i,3), ...
    mean(tmp), ...
    prctile(tmp,95), ...
    sum(tmp >= OUT.obs_ovec2(i)) ./ length(tmp))

    if isfield(OUT,'ba2'), try,fprintf('%s\t',OUT.ba2{i}),catch,end,end
    fprintf(1,'\n')
end

return
