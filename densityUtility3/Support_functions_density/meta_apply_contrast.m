function [icon,ctxtxi,betas,num_act,prop_act] = meta_apply_contrast(data,varargin)
%
% Fast form for use in Monte Carlo:
% [icon,ctxtxi] = meta_apply_contrast(data,ctxtxi)
%
% Full form used in setup: (meta_prob_activation)
% [icon,ctxtxi,betas,num_act,prop_act] = meta_apply_contrast(data,Xi,wts,contrasts)
% betas = OLS fits of indicators to data
%       weighted proportion of studies activating IF Xi are orthogonal.
% num_act = weighted number of studies activating 
% prop_act = weightd proportion of studies activating
%
% Used in Meta_activation_FWE and Meta_Task_Indicator
%
if isempty(varargin)
    error('Enter Xi, wts, contrasts OR ctxtxi matrix');
elseif length(varargin) == 1
    ctxtxi = varargin{1};
elseif length(varargin) == 3
    Xi = varargin{1};
    wts = varargin{2};
    contrasts = varargin{3};
else
    error('Enter Xi, wts, contrasts OR ctxtxi matrix');
end

if ~exist('ctxtxi','var')
    wts = wts ./ mean(wts);     % weights have mean 1
    wXi = diag(wts) * Xi;       % weighted indicator matrix
    
    % Note: the LM formulation gives the same result as scaling X to yield
    % proportion of activations by condition.  thus, betas are proportions by
    % condition IF Xi are orthogonal.
    
    xtxi = inv(wXi' * wXi) * wXi';  % for writing beta map
    if nargout > 2, betas = (xtxi * data')'; end
    
    if nargout > 3, num_act = (wXi' * data')'; end
    
    if nargout > 4, prop_act = num_act ./ repmat(sum(wXi),size(num_act,1),1); end
    
    ctxtxi = contrasts * xtxi;
end

icon = (ctxtxi * data')';      % icon: one col per contrast, rows are voxels

return
    
    % X is 1 -1 contrast matrix
    % Normalize so that contrast weights sum to 1
    % so that dat * X gives proportion of + vs. - studies at each voxel
    %X2 = X; X2(X>0)=X2(X>0)./sum(X(X>0));
    %X2(X<0)=X2(X<0)./-sum(X(X<0));
    