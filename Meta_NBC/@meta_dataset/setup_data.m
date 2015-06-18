function obj = setup_data(obj, terms)
% obj = setup_data(obj, terms)
%
% Given a meta_dataaset object and terms, do the following:
% 1) Select contrasts and voxels that meet cutoffs
% 2) Load term text files
% 3) Select studies that involve only one term
%
%
% Parameters used, with default values from meta_dataset constructor:
%
% min_vox: 5000
% term_freq_cutoff: 1.0000e-03
% min_prop_active: 0.0300


% assign all parameters in object to their respective variable names
N = fieldnames(obj.params);

for i = 1:length(N)
    eval([N{i} ' = obj.params.' N{i} ';'])
end

% ----------------------------------------------------------------------
% Database setup
% ----------------------------------------------------------------------

t1 = tic;
fprintf('Database setup\n');

fprintf('Original training set: %3.0f voxels x %3.0f studies\n', size(obj.dat));
fprintf('Selecting training set (studies and voxels) ')

% filter out studies with poor quality -- e.g., too few voxels
n_vox = sum(obj.dat, 1);
obj.include_cases = n_vox >= min_vox;

n_con = size(obj.dat, 2); % # of studies/contrasts

% filter out voxels with too few active studies
n_active = sum(obj.dat, 2)/n_con; % proportion of active studies
obj.include_vox = n_active > min_prop_active;

fprintf('%3.0f total contrasts.\n%3.0f meet min voxel criteria. %3.0f with no activations.\n', full(length(n_vox)), full(sum(obj.include_cases)), full(sum(n_vox == 0)))

if exist('test_data', 'var')
    test_data = test_data(obj.include_vox, :);
end

% later, we could set up with indices and exclude at train-time
obj.dat = obj.dat(obj.include_vox, obj.include_cases);
obj.connames = obj.connames(obj.include_cases);

%fprintf('Trimmed training set: %3.0f voxels x %3.0f studies\n', size(obj.dat));

% set classes
n_class = max(size(terms));
classes = zeros(n_con, 1);

fprintf(' %3.0f sec\n', toc(t1)); t1 = tic;

% Below was for Neurosynth. This is for Meta_Activation_FWE
ismetaactfwe = 1;

if ismetaactfwe
    
    [~, terms, classes] = string2indicator(obj.connames, terms);
    
    if sum(classes) == 0, error('No valid studies. Perhaps terms are entered incorrectly?'); end
    
    fprintf('%3.0f studies with excluded labels.\n', full(sum(classes == 0)))

else
    % ----------------------------------------------------------------------
    % Load terms and further select dataset
    % ----------------------------------------------------------------------
    
    % load terms
    fprintf('Loading %3.0f term lists', n_class);
    
    studynames = {};
    study_exclude = cell(1, n_class);
    [study_exclude{:}] = deal(false(n_con, 1)); % which to exclude
    
    for i = 1:n_class
        term = terms{i};
        termdata = importdata(sprintf('Input/%s.txt', term), '\t', 1);
        studynames{i} = termdata.textdata(2:end, 1);
        
        % use frequency cut-off
        freq = termdata.data(:,2)./termdata.data(:,1);
        
        % later, we could set up with indices and exclude at train-time
        study_exclude{i}(freq < term_freq_cutoff) = 1;
        studynames{i}(study_exclude{i}) = [];
        
    end
    
    fprintf(' %3.0f sec\n', toc(t1)); t1 = tic;
    fprintf('Selecting studies with only one term\n');
    
    % remove all studies that contain <> 1 term. note that this isn't
    % perfect because you could have a study that mentions one term above
    % frequency threshold and another just below threshold (e.g., 6
    % mentions versus 4 where the cut-off is 5). A better way would be to
    % use a more liberal exclusion than inclusion criterion. But that would
    % leave too few studies.
    for i = 1:n_class
        names = studynames{i};
        
        for j=1:n_class
            if i==j, continue; end
            names = setdiff(names, studynames{j});
        end
        
        [tf, ind] = ismember(names, obj.connames);
        valid = sum(ind>0);
        fprintf('%d studies contain term %s.\n', valid, terms{i});
        
        classes(ind(ind>0)) = i;
    end
    
    % end switch db type
end

% trim universe to studies with terms
obj.dat = obj.dat(:,classes > 0);
obj.connames = obj.connames(classes > 0);
classes(classes==0) = [];

obj.classes = classes;
obj.terms = terms;

fprintf(' Terms selected in %3.0f sec\n', toc(t1)); t1 = tic;
fprintf('Final training set (with exactly 1 term): %3.0f voxels x %3.0f studies\n', size(obj.dat));

end