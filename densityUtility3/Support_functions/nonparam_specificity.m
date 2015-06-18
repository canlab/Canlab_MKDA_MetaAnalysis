function [eff,p,sig,pt_given_a,pt,success,yp] = nonparam_specificity(y,X,iter,varargin)
% [eff,p,sig,pt_given_a,pt,success,yp] = nonparam_specificity(y,X,iter,[w])
%
% weights should be mean = 1

if length(varargin) > 0, w=varargin{1};, else, w=ones(size(y));, end

[n,k] = size(X);    % number of obs. and tasks

% weights
y = y .* w;

[eff,pt,pt_given_a,success] = ptask(y,X,n);


% initialize output
nulleff = zeros(iter,k);
yp = zeros(size(y,1),iter);

for i = 1:iter
    yp(:,i) = y(randperm(n));
    
    nulleff(i,:) = ptask(yp(:,i),X,n);
end

p = 1 - (sum(repmat(eff,iter,1) > nulleff) ./ iter);
    
sig = p <= .05;


return



function [eff,pt,pt_given_a,success] = ptask(activ_indicator,task_indicator,n)

% activ_indicator can be weighted

task_count = sum(task_indicator > 0);           % number of contrasts for each task type
pt = task_count ./ n;                           % prob of task x overall -- a fixed quantity depending on studies included


a = sum(activ_indicator);                       % activation count
success = activ_indicator' * task_indicator;    % actual number of activations for each task, weighted or not
                                                %  joint p(task, activ) *
                                                %  activ

pt_given_a = success ./ a;                      % prob of task x given activity in region i

eff = pt_given_a - pt;                          % effect size for specificity info given activation

return
    