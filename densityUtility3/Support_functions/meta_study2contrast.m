function DB = meta_study2contrast(DB)
% DB = meta_study2contrast(DB)
% helper function to turn unique studies stored in DB.Study
% into unique contrast codes stored in DB.Contrast
%
% in case contrasts are not manually entered.

% contrast numbers in order -- unique with no sorting
[u,ii]=unique(DB.Study); ii = sort(ii); u = DB.Study(ii);
for i=1:length(u)
wh = find(strcmp(DB.Study,u{i})); DB.Contrast(wh) = i;, end
DB.Contrast = DB.Contrast';
DB.Subjects = ones(size(DB.Contrast));

return
