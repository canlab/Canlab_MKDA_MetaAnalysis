function [nbc, data, group_data] = Meta_NBC(MC_Setup, terms, test_images)
% [nbc, data, group_data] = Meta_NBC(MC_Setup, terms, test_images)
%
% You will need:
% Inputs: a folder of text files with term labels
%
% MC_Setup = a database saved from an MKDA meta-analysis
% terms = {'pain', 'working.memory', 'emotion'};
% test_images = 'pain_test_data.img'
%

% Load group data up front, so we know if OK...
% --------------------------------------------------
if nargin > 2
    group_data = load_test_data(MC_Setup, test_images);
else
    group_data = [];
end

% Set up data object with default parameters
% --------------------------------------------------
data = meta_dataset('MC_Setup', MC_Setup);
clear MC_Setup

% Choose terms and restrict dataset
% --------------------------------------------------
data = setup_data(data, terms);

% Train classifier on all data
% --------------------------------------------------
nbc = train(meta_nbc, data);

% Test on same data to get apparent misclass rate
% --------------------------------------------------
nbc = test(nbc, data.dat);
nbc = get_error(nbc, data.classes);

% Display apparent correct classification rate
fprintf('\nApparent classification accuracy:\n---------------------------------------------\n')
nbc = report(nbc);

% Cross-validate and get cv error rate
% --------------------------------------------------
nbc = cv(nbc, data);

fprintf('\nCross-validated accuracy for meta-analysis dataset:\n---------------------------------------------\n')
nbc = report(nbc);

% Test on new image data -- if entered
% --------------------------------------------------

if nargin > 2
    
    nbc = test(nbc, group_data.dat(:, data.include_vox)', 1);
    
    test_class = ones(size(group_data.dat, 1), 1);
    
    nbc = get_error(nbc, test_class);
    
    fprintf('\nAccuracy for test data:\n---------------------------------------------\n')
    nbc = report(nbc);
end

end % main function




function group_data = load_test_data(MC_Setup, test_images)


t1 = tic;
disp('Loading mask from MC_Setup')
[discard, maskname, ee] = fileparts(MC_Setup.volInfo.fname);
maskname = [maskname ee];

if exist(maskname, 'file')
    mask = fmri_mask_image(which(maskname));
else
    error('Mask %s is not on path.', maskname);
end


disp('Reading data in mask space')
group_data = fmri_data(test_images, mask, 'sample2mask');


fprintf('Loaded in %3.0f sec\n', toc(t1)); 

end

% methods(group_data)
% plot(group_data);

% % RESAMPLE TEST DATA TO NEW SPACE
% 
% group_data = resample_space(group_data, volInfo);
% 
% fprintf('Resampled in %3.0f sec\n', toc(t1)); t1 = tic;
