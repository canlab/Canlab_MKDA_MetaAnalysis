function [cl,ptname,pname,pt] = dens_fdr_thresh(fname,histo,varargin)
% [cl,ptname,pname,pt] = dens_fdr_thresh(fname,histo,[doplot],[pthresh],[count data])
% 
% fname is name of image to threshold
% histo is vector of Ho (null hypothesis) hvalues of image
%   distribution of histo is used to convert image intensity values to
%   p-values.  Then fdr is applied, and a results image is written.
%
% needs at least 1000 samples in histo to work well.
%
% optional inputs:
% doplot:   make a plot
% pthresh:  FDR p-value threshold
% count data:  image of absolute study counts; threshold = at least 2
% contrasts

doplot = 0;
pthresh = .05;  % default
if length(varargin) > 0, doplot = varargin{1};,end
if length(varargin) > 1, pthresh = varargin{2};,end
if length(varargin) > 2, cdat = varargin{3};,end

% load image
V = spm_vol(fname); v = spm_read_vols(V);

% make histogram and smooth it some
[h,x]= hist(histo,100);
h = smooth_timeseries(h',10);

if doplot, tor_fig; plot(x,h./sum(h)); title('Null hypothesis density');drawnow,end

h2 = 1 - (cumsum(h) ./ sum(h)); % cum dist -- p-values

% look up p-values in h2
pmap = ones(size(v));

for i = 1:length(x)
    pmap(v >= x(i)) = h2(i);
end

pmap(v==0 | isnan(v)) = 1;

% write output image
[dd,ff,ee]=fileparts(fname);
pname = fullfile(dd,[ff '_pmap' ee]);
ptname = fullfile(dd,[ff '_pmap_fdr_thresh' ee]);
clname = fullfile(dd,[ff '_pmap_fdr_cl']);

V.fname = pname;
spm_write_vol(V,pmap);

% FDR threshold
tmp = pmap; tmp(v==0 | isnan(v))=[];   % eliminate out-of-analysis voxels
pt = FDR(tmp,pthresh);  clear tmp;

if isempty(pt), pt = 0;, end    % if no sig results

wh = logical(pmap <= pt);
mask = zeros(V.dim(1:3));
mask(wh) = 1;
    
% additional optional count masking.
if exist('cdat') == 1, mask(cdat < 2) = 0;, end

% write image
V.fname = ptname;
spm_write_vol(V,mask);

% get clusters
cl = mask2clusters(V.fname);
eval(['save ' clname ' cl']);

return
