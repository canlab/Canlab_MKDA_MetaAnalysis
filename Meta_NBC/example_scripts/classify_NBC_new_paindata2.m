%% LOAD DATABASE

clear all
t1 = tic;
load('/Users/tor/Documents/Tor_Documents/Presentations/2010/Columbia_fMRI_Prediction_Workshop_2010/NBC/full_data.mat', 'MC_Setup')

volInfo = MC_Setup.volInfo;
MC_Setup.cl = [];
MC_Setup.unweighted_study_data = logical(MC_Setup.unweighted_study_data);

fprintf('Loaded in %3.0f sec\n', toc(t1)); t1 = tic;

%% LOAD TEST DATA (IN DIFFERENT SPACE!)

clear a group_data obj
clear functions

load('/Users/tor/Documents/Tor_Documents/Presentations/2010/Columbia_fMRI_Prediction_Workshop_2010/group_scaling=none_amplitude_cond002_Pain/group_data.mat')
group_data.dat = group_data.X'; group_data.volInfo = group_data.mask.volInfo;
group_data.dat(isnan(group_data.dat)) = 0;

fprintf('Loaded in %3.0f sec\n', toc(t1)); t1 = tic;

methods(group_data)
plot(group_data);

%% RESAMPLE TEST DATA TO NEW SPACE

group_data = resample_space(group_data, volInfo);

fprintf('Resampled in %3.0f sec\n', toc(t1)); t1 = tic;

% could possibly rescale here ....

%% Classify the new test data

cd('/Users/tor/Documents/Tor_Documents/Presentations/2010/Columbia_fMRI_Prediction_Workshop_2010/NBC')

% Set up data object with default parameters
data = meta_dataset('MC_Setup', MC_Setup);

% Choose terms and restrict dataset
data = setup_data(data, {'pain', 'working.memory', 'emotion'});

% Train classifier on all data
nbc = train(meta_nbc, data);

% Test on same data to get apparent misclass rate
nbc = test(nbc, data.dat);
nbc = get_error(nbc, data.classes);

% Display apparent correct classification rate
nbc = report(nbc);

% Cross-validate and get cv error rate

nbc = cv(nbc, data);
nbc = report(nbc);

% Test on new image data
nbc = test(nbc, group_data.dat(data.include_vox, :));

test_class = ones(size(group_data.dat, 2), 1);
nbc = get_error(nbc, test_class);
nbc = report(nbc);

%% AN ALTERNATE WAY TO RUN

[meta_classes, test_probs] = multi_class_nbc_tor({'pain', 'working.memory', 'emotion'}, MC_Setup, group_data.dat, test_class);

