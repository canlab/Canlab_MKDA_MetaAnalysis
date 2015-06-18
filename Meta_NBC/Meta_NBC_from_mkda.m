function [nbc, data] = Meta_NBC_from_mkda(DB, MC_Setup, fieldname, names_cell)
% Take output from Meta_Setup and Meta_Activation_FWE (DB and MC_Setup) and
% run Naive Bayes classifier, on labels of your choice.
% 
% [nbc, data] = meta_nbc_from_mkda(MC_Setup, DB, fieldname, names_cell)
%
% If MC_Setup is empty, will run Meta_Activation_FWE to get it.
% 
% Examples:
% -----------------------------------------------------------------
% [nbc, data] = Meta_NBC_from_mkda(DB, MC_Setup, 'Valence', {'positive' 'negative'})
%
% fieldname = 'Stimuli';
% [nms, contrastcounts] = meta_explore_field(DB, fieldname);
% labels = nms(contrastcounts > 20)';
% [nbc, data] = Meta_NBC_from_mkda(DB, MC_Setup, fieldname, labels);
%
% Tor Wager, June 2012


if isempty(MC_Setup)
    MC_Setup = Meta_Activation_FWE('setup', DB); 
end

data = meta_dataset('MC_Setup', MC_Setup);
data.params.min_vox = 50; % do not exclude sparse maps, just those < 50 vox

wh_names = DB.(fieldname);
wh_names = wh_names(DB.pointind);
data.connames = wh_names;

% the below is not needed because setup_data does it
% 
% wh_include = false(size(wh_names));
% 
% for i = 1:length(names_cell)
%     wh_include(strcmp(wh_names, names_cell{i})) = true;
% end
% 
% if sum(wh_include) == 0
%     error('No studies meet criteria. Perhaps you entered the wrong labels?');
% else
%     fprintf('Selected %3.0f total studies.\n', sum(wh_include));
% end
%
% the below is not needed because setup_data does it 
% data.dat = data.dat(:, wh_include);
% data.connames = data.connames(wh_include);

data = setup_data(data, names_cell);
nbc = train(meta_nbc, data); % need this to define terms field for report

nbc = cv(nbc, data);
nbc = report(nbc);

end
