function DB = Meta_Chisq(DB,varargin)
% DB = Meta_Chisq(DB,[verbose],[mask image name],[control strings])
% 
% NEEDS:
%
% to set up design:
%DB.(fields)   % lists of fields containing task conditions for each coordinate point
%DB.pointind   % indices of which coord points are in which unique
%                contrast
%DB.X
%
% to run analysis:
% DB.PP         % list of contrast image names
% DB.studyweight  % weights
% 'Activation.img' or a mask input argument (2nd arg).
%
% 3rd to nth arguments are optional 'control strings'
% 'load'    : load from SETUP.mat file; do not recreate
% 'empty'   : to specify empty contrasts for each condition
%
% Examples:
% run chi-sq with mask image Activation.img thresholded at 0.05
% DB = Meta_Chisq(DB,2,.05) 
%
% To load already-created SETUP.mat and run:
% DB = Meta_Chisq(DB,2,[],'load')
% DB = Meta_Chisq(DB,2,'../Disgust/Activation_thresholded.img','load');
%
% To run, specifying empty contrasts
% DB = Meta_Chisq(DB,1,DB.maskname,'empty');

global dims
global cols

spm_defaults; defaults.analyze.flip = 0;
doload = 0; doempty = 0;

controlflags = varargin;    % to pass into subfunctions
if length(varargin) > 0, vb = varargin{1};, else, vb = 1;, end
% mask is varargin{2}
if length(varargin) > 2,
    for i = 3:length(varargin)
        if strcmp(varargin{i},'load'), doload = 1;, end
        if strcmp(varargin{i},'empty'), doempty = 1;, end
    end
end

% ---------------------------------------
% Set up design matrix
% names = contrast names, also output image names
% prunes DB
% ---------------------------------------
if doload
    load SETUP
else
    [X,names,DB,Xi,Xinms,condf,testfield] = Meta_Logistic_Design(DB,controlflags);
end

cols = size(X,2);

w = DB.studyweight;
w = w ./ mean(w);   
% weights should be mean 1 for logistic regression
% also preserve relative weighting

DB.X = X; DB.Xnames = names;
                 
% Save names of contrast images to use
DB.PP = check_valid_imagename(DB.PP);

contrast_images_used = DB.PP;

if ~doload
    save SETUP contrast_images_used X names DB w Xi Xinms condf testfield
end

% ---------------------------------------
% Set up images
% ---------------------------------------

% image dimensions, etc.
V = spm_vol(DB.PP(1,:));
Vall = spm_vol(DB.PP);

dims = V.dim(1:3);

if vb > 0, fprintf(1,'\nNew image dims: %3.0f %3.0f %3.0f %3.0f ',dims(1), dims(2), dims(3), cols),end


% ---------------------------------------
% brain masking
% e.g., could mask with overall activation image
% ---------------------------------------

mask = []; 
if length(varargin) > 1, 
    mask = varargin{2};, 
    if isstr(mask), Vm = spm_vol(mask);, 
        if vb > 0, fprintf(1,'\nMask is: %s\n');, end
        mask = spm_read_vols(Vm);,
            % if a mask file, then load it
    elseif length(mask) == 1
        % if a number, treat this as threshold value for pa_overall
        % I recommend .001, which is a 0.1% chance of activating across all
        % tasks

        if ~(exist('Activation.img') == 2),
            try
                Vm = spm_vol('../Activation.img'); 
            catch
                P = spm_get(1,'Select Activation.img');
                Vm = spm_vol(P);
            end
        else
            Vm = spm_vol('Activation.img');
        end
        maskimg = spm_read_vols(Vm);,
        mask = real(maskimg > mask);  % threshold 
    end
end
if isempty(mask), if vb > 0, fprintf(1,'\nNo mask. Running all voxels!\n');, end, mask = ones(dims);,end

nvox = mask(:); nvox = sum(nvox>0);
fprintf(1,'\nVoxels in mask: %3.2f\n',nvox);





% -------------------------------------------------------------------
% * for each slice...
% -------------------------------------------------------------------

fprintf(1,'\nRunning Chi-square Slice %03d',0)

for slicei = 1:dims(3)
    
    if vb > 1, t1 = clock;, fprintf(1,'\b\b\b%03d',slicei),end
    
    %if vb > 0, fprintf(1,'      '); , end   % to offset Load prompt
    
    % returns image names
    process_slice(slicei,Vall,X,w,names,vb,mask(:,:,slicei),[],Xi,Xinms,condf);


    %if vb > 1, fprintf(1,'\t%6.0f s for slice',etime(clock,t1)),end
end





return







% -------------------------------------------------------------------
%
%
%
%
% * SUB-FUNCTIONS
%
%
%
%
% -------------------------------------------------------------------    
    
    
function process_slice(slicei,Vall,X,w,names,vb,varargin)

mask = []; if length(varargin) > 0, mask = varargin{1};, end
sl = []; if length(varargin) > 1, sl = varargin{2};, end

% set up weighted average calc
Xi = []; 
if length(varargin) > 2, 
    Xi = varargin{3};, Xinms = varargin{4};, 
    W = diag(w); % ./ sum(w));
    Xi = Xi' * W;           % multiply this by data to get weighted avgs
    sumxi = sum(Xi,2);
end

condf = []; if length(varargin) > 4, condf = varargin{5};, end

global dims
global cols

nimgs = length(Vall);
ncontrasts = size(Xi,2);

% if empty images specified, ncontrasts will be greater than nimgs
doempty = ncontrasts - nimgs;
if doempty > 0
    addtoy = zeros(doempty,1);
end

% output images -- initialize
% --------------------------------------------

omnchi2 = NaN * zeros([dims(1:2)]);
chi2pmap = ones([dims(1:2)]);
chi2warn = omnchi2;
chi2nonp = omnchi2;
emptyimg = NaN .* zeros(dims);     % in case we need to create a new volume, used later as well

if ~isempty(Xi), avgs = NaN * zeros([dims(1:2) size(Xi,1)]);, end



% -------------------------------------------------------------------
% * load the slice
% -------------------------------------------------------------------
et = clock;

if isempty(sl)

    if vb > 0, fprintf(1,'Load >'), end

    if ~isempty(mask) && ~(sum(sum(sum(mask))) > 0)
        % skip it; write only
        % -------------------------------------------------------------------
        % * write output images
        % -------------------------------------------------------------------

        warning off
        write_beta_slice(slicei,Vall(1),omnchi2,emptyimg,'chi2_',{'map'});
        write_beta_slice(slicei,Vall(1),chi2pmap,emptyimg,'chi2p_',{'map'});
        write_beta_slice(slicei,Vall(1),chi2warn,emptyimg,'chi2_warning_',{'map'});
        write_beta_slice(slicei,Vall(1),chi2nonp,emptyimg,'chi2_nonpar_',{'map'});

        if ~isempty(Xi)
            class_avg_images = write_beta_slice(slicei,Vall(1),avgs,emptyimg,'avgs_',Xinms);
            save SETUP -append class_avg_images
        end
        
        return
    else
        sl = timeseries_extract_slice(Vall,slicei);
    end

    if vb > 0, fprintf(1,'\b\b\b\b\b\b');,end

else
    % we've already loaded the slice!

end



if ~isempty(mask), sl(:,:,1) = sl(:,:,1) .* mask;, end


% -------------------------------------------------------------------
% * find in-mask voxels
% -------------------------------------------------------------------

sumsl = sum(sl,3);  % sum of all contrasts

% find indices i and j for all in-mask voxels
wvox = find(sumsl ~= 0 & ~isnan(sumsl));
[i,j] = ind2sub(size(sl(:,:,1)),wvox);

%fprintf(1,[repmat('\b',1,19) '%3.0f voxels > '],length(i))
if vb > 0, 
    if slicei > 1, fprintf(1,'\b\b\b\b\b\b\b\b\b\b\b\b\b\b');, end
    str = sprintf(' elapsed: %3.0f s vox:%04d fitting',etime(clock,et),length(i));
    fprintf(1,str);
end

% -------------------------------------------------------------------
% * compute regression for each in-mask voxel
% -------------------------------------------------------------------



% loop through voxels
% --------------------------------------------
et = clock;
warning off     % because of iteration limit
    
fprintf(1,'%04d',0);

for k = 1:length(i)

    y = squeeze(sl(i(k),j(k),:)); 
    
    if doempty > 0      % empty contrasts
        y = [y; addtoy];
    end
    
    % convert from weighted to indicator (on/off)
    % only works for spherical convolution!!
    % * no longer necessary, because images are no longer weighted.
    y = double(y>0);

    % chi-square test, weighted
    % last argument is nonparametric flag to test close results more
    % carefully
    % --------------------------------------------
    [chi2,df,chi2p,sig,warn,tab,e,isnonpar] = chi2test([y condf],'obs',w,1);
 
    omnchi2(i(k),j(k)) = chi2;
    chi2pmap(i(k),j(k)) = chi2p; 
    chi2warn(i(k),j(k)) = warn;
    chi2nonp(i(k),j(k)) = isnonpar;
    
    % --------------------------------------------
    
    % weighted average probability images
    % --------------------------------------------
    if ~isempty(Xi)
        avg = Xi * y; avg = avg ./ sumxi;
        avgs(i(k),j(k),:) = avg;
    end
    
        
    
    if mod(k,10)==0, fprintf(1,'\b\b\b\b%04d',k);,end
    %if k == 1000, fprintf(1,'%3.0f s per 1000 vox.',etime(clock,et)), end
end
fprintf(1,'\b\b\b\b');

warning on

if vb > 1, fprintf(1,[repmat('\b',1,14+9) ' elapsed: %3.0f s'],etime(clock,et));, end

clear sl



% -------------------------------------------------------------------
% * write output images
% ------------------------------------------------------------------- 

warning off

write_beta_slice(slicei,Vall(1),omnchi2,emptyimg,'chi2_',{'map'});
write_beta_slice(slicei,Vall(1),chi2pmap,emptyimg,'chi2p_',{'map'});
write_beta_slice(slicei,Vall(1),chi2warn,emptyimg,'chi2_warning_',{'map'});
write_beta_slice(slicei,Vall(1),chi2nonp,emptyimg,'chi2_nonpar_',{'map'});

if ~isempty(Xi)
    class_avg_images = write_beta_slice(slicei,Vall(1),avgs,emptyimg,'avgs_',Xinms);
    save SETUP -append class_avg_images
end

if length(i) > 0, save SETUP -append df, end

warning on

if vb > 0, fprintf(1,[repmat('\b',1,14+9+1)]);, end
    
    
return







function Pw = write_beta_slice(slicei,V,betas,emptyimg,varargin)
% Pw = write_beta_slice(sliceindex,V,data,empty,prefix)
% Slice-a-metric version
warning off % due to NaN to int16 zero conversions
V.dim(4) = 16; % set to float to preserve decimals

prefix = 'check_program_';
names = {''};               % for field only
if length(varargin) > 0, prefix = varargin{1};,end      % field name
if length(varargin) > 1, names = varargin{2};,end       % level names

for voli = 1:size(betas,3) % for each image/beta series point
    %if voli < 10, myz = '000';, elseif voli < 100, myz = '00';, else myz = '000';,end
    V.fname = [prefix names{voli} '.img'];
    V.descrip = ['Created by Meta_Logistic'];

    % create volume, if necessary
    if ~(exist(V.fname) == 2), spm_write_vol(V,emptyimg);,end
        
    spm_write_plane(V,betas(:,:,voli),slicei);
    
    if ~(exist('Pw')==1), Pw = which(V.fname);, 
    else Pw = str2mat(Pw,which(V.fname));
    end
    
end



warning on
return