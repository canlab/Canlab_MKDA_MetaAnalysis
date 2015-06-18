%% LOAD DATABASE

clear all
t1 = tic;
% load('/Users/tor/Documents/Tor_Documents/Presentations/2010/Columbia_fMRI_Prediction_Workshop_2010/NBC/full_data.mat', 'MC_Setup')
% 
% volInfo = MC_Setup.volInfo;
% MC_Setup.cl = [];
% MC_Setup.unweighted_study_data = logical(MC_Setup.unweighted_study_data);

load('/Users/tor/Documents/Tor_Documents/Presentations/2010/Columbia_fMRI_Prediction_Workshop_2010/NBC/MC_Setup.mat')

fprintf('Loaded in %3.0f sec\n', toc(t1)); t1 = tic;

%% RUN IT
% This function runs several simple commands (you can look at the function
% to see how to use them flexibly) to run a Naive Bayes Classifier analysis
% on the meta-analytic dataset, including 10-fold cross-validation and
% prediction of new test data from a set of images.

% You will also need:
% Inputs: a folder of text files with term labels

terms = {'pain', 'working.memory', 'emotion'};
test_images = fullfile(pwd, 'example_data', 'pain_test_data.img');  % a set of 17 pain images from the Wager lab NSF project

[nbc, data, group_data] = Meta_NBC(MC_Setup, terms, test_images);


%% NOW LOAD NEW TEST DATA 

% Load new group data in the same space as the mask

% If these were image files, we would do this:
% test_images = files_here....
% 
% [dummy, maskname, ee] = fileparts(MC_Setup.volInfo.fname);
% maskname = [maskname ee];
% mask = fmri_mask_image(which(maskname));
% 
% group_data = fmri_data(test_images, mask, 'sample2mask');

load('/Users/tor/Documents/Tor_Documents/Presentations/2010/Columbia_fMRI_Prediction_Workshop_2010/NBC/example_data/latest_imgs.mat', 'testdata', 'classes')

%% APPLY THE ALREADY-TRAINED CLASSIFIER

% Test on new image data
nbc = test(nbc, testdata(data.include_vox, :));

test_class = classes;
nbc = get_error(nbc, test_class);
nbc = report(nbc);

%% IT is onvenient to save the test data in an object...
% then we can do things with it more easily

test_data = fmri_data; test_data = create(test_data, 'dat', testdata);
test_data = create(test_data, 'Y', classes);
test_data.source_notes = 'Data from Tal, pain wm emo, 10/26/10';

clear testdata

% cd example_data
% save latest_test_data test_data;
% cd ..

plot(test_data);

%% Now center the images and try classification again

test_data = rescale(test_data, 'centerimages');

nbc = test(nbc, test_data.dat(data.include_vox, :));
nbc = get_error(nbc, test_data.Y);
nbc = report(nbc);


%% Now z-score the images and try classification again

test_data = rescale(test_data, 'zscoreimages');

nbc = test(nbc, test_data.dat(data.include_vox, :));
nbc = get_error(nbc, test_data.Y);
nbc = report(nbc);

