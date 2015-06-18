function [MC_Setup,activation_proportions,icon] = meta_prob_activation(DB,Xi,contrasts,connames,Xinms)
% [MC_Setup,activation_proportions,icon] = meta_prob_activation(DB,[X indicator mtx],[contrasts],[con. names],[Xinms])
%
% This function writes 'Activation_proportion.img'
% --a weighted map of the proportion of independent contrasts activating
% each voxel in the brain --
%
% ...And it returns MC_Setup setup information for the whole-brain FWE
% corrected Monte Carlo
%
% Additional input arguments Xi...contrasts...etc will create maps by
% condition, and contrasts across conditions
% This is a SECOND KIND of contrast: contrasts across study types
% Use Meta_Logistic_Design to set up these inputs
%
% This function is part of Meta_Activation_Fwe
%
% uses:
% DB.maskname
% DB.radius
% DB.Contrast      (study ID number; contrast within study)
% DB.x, DB.y DB.z
% DB.studyweight
%
% Tor Wager, May 2006

% Programmers' notes:
% 4.7.2013  There was a bug in Meta_Activation_FWE that made it not robust to using 
% non-ascending numerical order of contrasts in DB.Contrast.  
% The function meta_prob_activation was reconstructing the maps in ascending order of 
% contrasts, returning maps in ascending sorted contrast order in 
% MC_Setup.unweighted_study_data.  MC_Setup.Xi, the indicators for task type, are 
% sorted in the order of contrasts entered in the database, which respects the order 
% in DB.pointind and DB.connumbers from Meta_Setup.  Tor changed the sort order in 
% meta_prob_activation to use DB.connumbers.  
% 
% This affects MKDA difference analyses (contrasts) when contrast numbers are not 
% entered in ascending numerical order in your database. It also affects classification 
% with Meta_NBC_from_mkda, but not Meta_SVM_from_mkda.


if ~isfield(DB, 'maskname'), error('Missing DB.maskname.  Try running Meta_Setup.'); end
if ~isfield(DB, 'radius'), error('Missing DB.radius.  Try running Meta_Setup.'); end
if ~isfield(DB, 'Contrast'), error('Missing DB.Contrast.  Try running Meta_Setup.'); end
if ~isfield(DB, 'x'), error('Missing DB.x,y,z.  Try running Meta_Setup.'); end
if ~isfield(DB, 'studyweight'), error('Missing DB.studyweight.  Try running Meta_Setup.'); end

% set up mask
str1 = sprintf('Reading mask data'); fprintf(1,str1);

%[maskdata,xyzlist,V] = meta_read_mask(DB.maskname);

% the line above works, but this has the extra fields we need.
[V, maskdata] = iimg_read_img(DB.maskname, 1);
xyzlist = V.xyzlist;
maskdata = iimg_reconstruct_3dvol(maskdata, V);

erase_string(str1);

str1 = sprintf('Setting up contrasts and initializing output'); fprintf(1,str1);

% ---------------------------------------
% set up contrasts and sizes
% ---------------------------------------
% cons = unique(DB.Contrast, 'stable');     % contrasts (studies)
%  tor: 4/7/2013: added 'stable' : to match order returned in Meta_Logistic_Design
%  but this does not always return same order as DB.connumbers
cons = DB.connumbers;

nc = length(cons);              % number of contrasts
v = size(xyzlist, 1);            % voxels in mask
r = DB.radius;                  % radius in voxels

%ivectors = zeros(v,nc,'uint8');       % initialize output
ivectors = sparse(v, nc);
n = zeros(1, nc);

% set up weights (returns empty if no weights, or wt values for each con)
wts = setup_weights(DB, nc);

erase_string(str1);

% ---------------------------------------
% Convolve activation images with kernel
% ---------------------------------------
str1 = sprintf('Convolving %03d contrast maps: %03d',nc,0);
fprintf(1,str1);

for c = 1:nc

    fprintf(1,'\b\b\b%03d', c);
    % get actual voxel coords for study (in-mask only)
    xyzvox = meta_get_voxel_coords(DB, cons(c), V, maskdata);

    % save number for Monte Carlo
    n(c) = size(xyzvox, 1);

    % get convolved indicator function
    %ivectors(:,c) = meta_fast_sphere_conv(xyzlist,xyzvox,r);  OLD, not
    %used.
    ivectors(:,c) = iimg_xyz2spheres(xyzvox, xyzlist, r);
end

erase_string(str1);

MC_Setup.unweighted_study_data = ivectors;

% ---------------------------------------
% apply weights (fastest method)
% ---------------------------------------

str1 = sprintf('Applying weights'); fprintf(1, str1);

for i=1:nc, ivectors(:,i) = ivectors(:,i) .* wts(i); end

erase_string(str1);

% ---------------------------------------
% get summary map and reconstruct into 3D
% ---------------------------------------
str1 = sprintf('Constructing average map.'); fprintf(1,str1);
activation_proportions = sum(ivectors, 2);
erase_string(str1);
maskdata = meta_reconstruct_mask(activation_proportions, xyzlist, V.dim(1:3), 1, V, 'Activation_proportion.img');

% ---------------------------------------
% Contrasts across conditions
% (Apply weights in process of contrast computation
% for compatibility with MC)
% ---------------------------------------
if nargin > 1

    if nargin < 5, Xinms = define_Xinms(Xi); end
        
    if nargin < 4, [contrasts,connames] = meta_enter_contrasts(Xi); end
    str1 = sprintf('Computing contrast maps.'); fprintf(1,str1);
    
    [icon,ctxtxi,prop_by_condition] = meta_apply_contrast(MC_Setup.unweighted_study_data, Xi, wts, contrasts);

    for i = 1:size(prop_by_condition,2)  % for each condition, write map
        meta_reconstruct_mask(prop_by_condition(:,i), xyzlist, V.dim(1:3), 1, V, [Xinms{i} '.img']);
    end
    
    for i = 1:size(icon,2)  % for each contrast, write map
        meta_reconstruct_mask(icon(:,i), xyzlist, V.dim(1:3), 1, V, [connames{i} '.img']);
    end
    
    erase_string(str1);
end

% ---------------------------------------
% save output for MC simulation
% ---------------------------------------
%MC_Setup.xyzlist = xyzlist;
MC_Setup.volInfo = V;
MC_Setup.n = n;
MC_Setup.wts = wts;
MC_Setup.r = r;

if nargin > 1
    MC_Setup.contrasts = contrasts;
    MC_Setup.connames = connames;
    MC_Setup.Xi = Xi;
    MC_Setup.Xinms = Xinms;
    MC_Setup.ctxtxi = ctxtxi;
    MC_Setup.con_data = icon;
end




return





function wts = setup_weights(DB,nc)
if isfield(DB,'studyweight') && ~isempty(DB.studyweight);
    wts = DB.studyweight;
    if length(wts) ~= nc,
        error('Length of study weight vector DB.studyweight does not match number of contrasts.');
    end
end
return


function erase_string(str1)
fprintf(1,repmat('\b',1,length(str1))); % erase string
%fprintf(1,'\n');
return


function define_Xinms(Xi)
n = size(Xi,2);
Xinms = cell(1,n);
for i = 1:n
    Xinms{i} = input(['Enter name for condition ' num2str(i)],'s');
end
return
