function cl = pmap_threshold(pimg,searchmask,pt,varargin)
% cl = pmap_threshold(pimg,mask image,p-threshold,[image with values to extract])
%
% Threshold a p-value image (pimg) and display suprathreshold blobs
% Only display significant values included in the mask image (img2)
% 
% p-threshold is the p-value threshold to use (e.g., .005)
% if no threshold is specified, the default is to use FDR correction
%
% If you enter a 4th [optional] argument, the function extracts data 
% from the image specified in this argument, e.g., 't-value_image.img'
% This is useful because if p-values can correspond to increases or
% decreases (i.e., the test is 2-tailed), then you can visualize positive
% effects in red and negative ones in blue.
%
% pmap_threshold('Activation_thresholded_pmap.img',DB.maskname);
% pmap_threshold('Activation_thresholded_pmap.img',DB.maskname,.005);
% cl = pmap_threshold('chi2p_Omnibus.img',DB.maskname,.05,'chi2_Omnibus.img');

if nargin == 0
    % interactive mode
    pimg = spm_get(1,'*img','Select p-value image.',pwd);
    searchmask = spm_get(Inf,'*img','Select mask image or DONE.',pwd);
    pt = spm_input(['Enter p-value threshold or 1 for FDR: ']);
    if pt == 1, pt = [];,end
    varargin{1} = spm_get(Inf,'*img','Select image with data values or DONE.',pwd);
    if isempty(varargin{1}), varargin = {};, end
end


cl = [];

if ~(exist('searchmask') == 1), searchmask = [];, end
if ~(exist('pt') == 1), pt = [];, end

V = spm_vol(pimg);
p = spm_read_vols(V);

if ~isempty(searchmask)
    m = spm_read_vols(spm_vol(searchmask));
    m(isnan(m)) = 0;
    p = p + double(~m);                     % make out of mask values too big
end
    
tmp = p(:); tmp(p > 1-10*eps | isnan(p))=[];   % eliminate out-of-analysis voxels

nvox = length(tmp);

if isempty(pt)
    % FDR threshold
    pt = FDR(tmp,.05);                  % threshold  

    if isempty(pt), pt = -Inf;, end
    fdrstr = 'FDR';
else
    fdrstr = 'UNC';
end
    
% thresholded p image
p(p > pt) = NaN;

minp = min(tmp(tmp>0));      % replace zeros with min p-value so it extracts OK (just technical for display)
p(p<eps) = minp;

p = -log10(p);               % convert for display

sigvox = sum(tmp<pt);

fprintf(1,'Voxels: %3.0f, Threshold: %s p < %3.4f, Sig. Voxels = %3.0f\n',nvox,fdrstr,pt,sigvox);

if sigvox == 0, return, end


[d,f,e] = fileparts(V.fname); V.fname = fullfile(d,[f '_' fdrstr e]); V.descrip = ['-log10(pvalue)'];
spm_write_vol(V,p);

if length(varargin) > 0
    cl = mask2clusters(V.fname,varargin{1});
else
    cl = mask2clusters(V.fname);
end

cluster_orthviews(cl);

% set up interactive bar plotting
%bar_interactive;

return

