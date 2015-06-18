function [nms, contrastcounts] = meta_explore_field(DB, fieldname)
% [nms, contrastcounts] = meta_explore_field(DB, fieldname)
%
% Print a table of the peak and contrast counts for unique values of a
% field.
%
% e.g., [nms, contrastcounts] = meta_explore_field(DB, 'Valence')

labels = DB.(fieldname);

nms = unique(DB.(fieldname));
labels2 = labels(DB.pointind);

fprintf('Label\tPeaks\tContrasts\n')

for i = 1:length(nms)
    contrastcounts(i, 1) = sum(strcmp(labels2, nms{i}));
    
    fprintf('%s\t%3.0f\t%3.0f\n', nms{i}, sum(strcmp(labels, nms{i})), contrastcounts(i));
end

end
