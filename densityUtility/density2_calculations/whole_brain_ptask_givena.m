function OUT = whole_brain_ptask_givena(OUT,varargin)
% OUT = whole_brain_ptask_givena(OUT,verbose level, mask image or threshold value(for pa_overall) )
% OUT = whole_brain_ptask_givena(OUT,2,.001)
% verbose flag
if length(varargin) > 0, vb = varargin{1};, else, vb = 1;, end

P = OUT.PP;
task_indicator = OUT.allcondindic;
field = OUT.testfield;
names = OUT.allcondnames;

% image dimensions, etc.
V = spm_vol(P(1,:));
Vall = spm_vol(P);

global dims
global cols

dims = V.dim(1:3);

if vb > 0, fprintf(1,'\n\t\tNew image dims: %3.0f %3.0f %3.0f %3.0f ',dims(1), dims(2), dims(3), cols),end

% ---------------------------------------
% brain masking
% ---------------------------------------

mask = []; 
if length(varargin) > 1, 
    mask = varargin{2};, 
    if isstr(mask), Vm = spm_vol(mask);, mask = spm_read_vols(Vm);,
            % if a mask file, then load it
    elseif length(mask) == 1
        % if a number, treat this as threshold value for pa_overall
        % I recommend .001, which is a 0.1% chance of activating across all
        % tasks
        Vm = spm_vol('pa_overall.img');  maskimg = spm_read_vols(Vm);,
        mask = real(maskimg > mask);  % threshold 
    end
end
if isempty(mask), mask = ones(dims);,end


% -------------------------------------------------------------------
% * define output image names
% -------------------------------------------------------------------
for i = 1:length(names)
    names{i} = [field '_' names{i}];
end


                                   
% -------------------------------------------------------------------
% * for each slice...
% -------------------------------------------------------------------

for slicei = 1:dims(3)
    
    if vb > 1, t1 = clock;, fprintf(1,'\nSlice %3.0f \n------------>\n ',slicei),end
    
    % returns image names
    [x2,x2p,eff,z,p] = process_slice(slicei,Vall,task_indicator,field,names,vb,mask(:,:,slicei));

    if ~isempty(x2), x2f = x2;, end     % final versions
    if ~isempty(x2p), x2pf = x2p;, end 
    if ~isempty(eff), efff = eff;, end 
    if ~isempty(z), zf = z;, end
    if ~isempty(p), pf = p;, end 
    
    if vb > 1, fprintf(1,'\t%6.0f s for slice',etime(clock,t1)),end
end

OUT.x2 = x2f;
OUT.x2p = x2pf;
OUT.eff = efff;
OUT.z = zf;
OUT.p = pf;

str = ['save STATS_VARS_' field ' OUT'];
eval(str)

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
    
    
function [x2,x2p,eff,z,p] = process_slice(slicei,Vall,task_indicator,field,names,vb,varargin)

x2 = []; x2p = []; eff = []; z = []; p = [];
mask = []; if length(varargin) > 0, mask = varargin{1};, end
sl = []; if length(varargin) > 1, sl = varargin{1};, end

global dims
global cols

% -------------------------------------------------------------------
% * load the slice
% -------------------------------------------------------------------

if isempty(sl)
    
if vb > 0, fprintf(1,'\tLoading data >'), end 
et = clock;
if ~isempty(mask) & ~(sum(sum(sum(mask))) > 0)
    % skip it
    fbetas = NaN * zeros([dims(1:2) cols]);
    ntrimmed = NaN;
    if vb > 1, fprintf(1,'...Empty slice...'), end
    return
else
    sl = timeseries_extract_slice(Vall,slicei);
end

if vb > 1, fprintf(1,'loaded in %3.2f s.',etime(clock,et)),end

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

fprintf(1,'\n\tCalculating P(Task|Activity) for %3.0f voxels > ',length(i))


% -------------------------------------------------------------------
% * compute p(T|A) for each in-mask voxel
% -------------------------------------------------------------------
nsize = size(sl); 
nsize = [nsize(1:2) size(task_indicator,2)];    % vox x vox x task condition

% output images -- initialize
x2field = zeros(nsize(1:2));     % chi2 overall for test field
pfield  = zeros(nsize(1:2));     % p-value

levels = zeros(nsize);           % effect size for each level
levelz = zeros(nsize);           % z-score   
levelp = zeros(nsize);           % p-value

et = clock;
for k = 1:length(i)

    y = squeeze(sl(i(k),j(k),:)); 
    
    % could use optional arg here to use %, but now treats everything as
    % activated or not
    [pt_given_a,eff,z,p,stat] = p_task_given_activation(task_indicator,y);
    
    x2field(i(k),j(k)) = stat.x2;
    pfield(i(k),j(k)) = stat.x2p;
    
    levels(i(k),j(k),:) = eff';
    levelz(i(k),j(k),:) = z';
    levelp(i(k),j(k),:) = p';
   
    if k == 1000, fprintf(1,'%3.0f s per 1000 vox.',etime(clock,et)), end
end


clear sl



% -------------------------------------------------------------------
% * write output images
% ------------------------------------------------------------------- 

emptyimg = zeros(dims);     % in case we need to create a new volume, used later as well
    
warning off
if vb > 1, fprintf(1,'\n\tWriting f (filtered) plane > '), end

x2 = write_beta_slice(slicei,Vall(1),x2field,emptyimg,[field '_chi2_']);
x2p = write_beta_slice(slicei,Vall(1),pfield,emptyimg,[field '_chi2p_']);
eff = write_beta_slice(slicei,Vall(1),levels,emptyimg,[field '_p(t|a)_eff_'],names);
z = write_beta_slice(slicei,Vall(1),levelz,emptyimg,[field '_p(t|a)_z_'],names);
p = write_beta_slice(slicei,Vall(1),levelp,emptyimg,[field '_p(t|a)_p_'],names);

warning on

    
    
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
    V.descrip = ['Created by whole_brain_ptask_givena '];

    % create volume, if necessary
    if ~(exist(V.fname) == 2), spm_write_vol(V,emptyimg);,end
        
    spm_write_plane(V,betas(:,:,voli),slicei);
    
    if ~(exist('Pw')==1), Pw = which(V.fname);, 
    else Pw = str2mat(Pw,which(V.fname));
    end
    
end



warning on
return

