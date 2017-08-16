%% Set up the analysis
% ------------------------------------------------------------------------

% Add the MKDA Toolbox with subfolders to your Matlab path:
% (Customize folder name for your computer)
g = genpath('/Users/tor/Documents/Code_Repositories/Canlab_MKDA_MetaAnalysis');
addpath(g)
savepath

% Create a directory for the meta-analysis results and go to it
mkdir('tmp_results'); 
cd('tmp_results')

% Load the sample database. 
% When prompted, type Agency_meta_analysis_database to save loaded
% variables in a file with this name.
dbname = 'Agency_meta_analysis_database.txt';
read_database

% Run Meta_Setup, which checks and sets up important variables
% Save setup database structure DB
DB = Meta_Setup(DB);

% Set up the MKDA analysis and data matrix, save setup files
Meta_Activation_FWE('setup', DB)

% Add 100 iterations to the permutation test:
num_iterations = 100;
Meta_Activation_FWE('mc', num_iterations);

% Run it again to add 100 more
Meta_Activation_FWE('mc', num_iterations);

% At least 5,000, or preferably 10,000, iterations is recommended.

% Get results
Meta_Activation_FWE('results', 1);

% You can also do this to run all steps at once: 
% Meta_Activation_FWE('all', DB, 5000);

% Make images of the results: Montages and surface rendering
% There are many options using the CANlab "Core Tools" and display
% functions in the CANlab mediation toolbox.  Here are just a couple of
% quick examples:

slices_fig_h = cluster_orthviews_montage(6, 'axial');
slices_fig_h = cluster_orthviews_montage(6, 'coronal');
hh = findobj(gcf, 'Type', 'text', 'Color', 'g')
set(hh, 'FontSize', 16)
?
load('Activation_clusters.mat')
all_surf_handles = mediation_brain_surface_figs(cl(1:3), []);

