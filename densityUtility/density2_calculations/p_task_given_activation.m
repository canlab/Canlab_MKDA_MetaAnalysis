function [pt_given_a,eff,z,p,stat] = p_task_given_activation(task_indicator,activ_indicator,varargin)
% [pt_given_a,eff,z,p,STAT] =
% p_task_given_activation(task_indicator,activ_indicator,varargin)
%
% Binomial or one-proportion z-test for a higher proportion than expected
% given activation in a region
%
% inputs: task_indicator:   matrix of contrasts x tasks, with nonzero entries signifying task
%                           identity.
%         activ_indicator:  column vector of contrasts, with nonzero
%                           entries indicating activation in a region
%         weighting:        optional 1/0 flag indicating to weight by
%                           entry values in task_indicator
%                           these may, for example, contain sample sizes
%                           for each contrast
%         flag: suppress binarizing activ_indicator (good for whole brain
%         analysis)
%
%
% ouputs:
% 
% eff: pt_given_a - pt, the increase (or decrease) in likelihood of being
%   in task state t given observed activation in region r
%   positive values for eff are an increase in p(t|a) by x% over
%   the base-rate for the task (p(t)).  
%   negative values mean that the task is LESS likely than expected
%   from the base rate by x%
%
% Examples:
% [pt_given_a,eff,z,sd,p] = p_task_given_activation(newti,rmatx(:,i));
%
% Used in:
% dbcluster_contrast_table and whole_brain_ptask_givena


% ---------------------------------------------
% values we need to calculate preliminarily
% ---------------------------------------------

if length(varargin) > 1
    % do nothing
else
    activ_indicator = real(activ_indicator > 0);          % otherwise we weight by # points as well!
end

a = sum(activ_indicator);                       % activation count
activ = a ./ length(activ_indicator);           % proportion active contrasts overall; expected proportion under Ho
  
% set weighting and adjust task indicator to make sure it's 1's and 0's if
% no weighting is desired
if length(varargin) > 0, doweight = varargin{1};, else, doweight = 0;, task_indicator = real(task_indicator > 0);, end

task_count = sum(task_indicator > 0);           % number of contrasts for each task type

normf = sum(task_indicator) ./ task_count;      % normalize by average; so use RELATIVE sample sizes
success = activ_indicator' * task_indicator ./ normf;    % actual number of activations for each task, weighted or not




pt = task_count ./ size(task_indicator,1);     % prob of task x overall -- a fixed quantity depending on studies included

pt_given_a = success ./ a;                      % prob of task x given activity in region i

eff = pt_given_a - pt;                          % effect size for specificity info given activation
    

expected_s = task_count .* activ;                       % expected number of activations for each task
expected_f = task_count .* (1-activ);                   % expected non-activations

% ---------------------------------------------
% omnibus test for effects of indicators
% assumes contrasts can belong to only one task (independent tasks)
% ---------------------------------------------
x2 = sum( (success - expected_s) .^ 2 ./ expected_s );
x2df = length(expected_s) - 1;
x2p = 1 - chi2cdf(x2,x2df);

% ---------------------------------------------
% tests on individual task indicators
% ---------------------------------------------

% to use normal approx,
% # of expected successes and failures must each => 10.  
              
if any(expected_s) < 10 | any(expected_f) < 10

    % use binomial distribution, n is task_count, p is successes/task count
    meth = 'binomial';
    
    for i = 1:length(task_count)        % for each task type
        % compare successes (act | task) to total attempts with task
        % against Ho proportion of overall successes blind to task
        pp = binocdf(success(i),task_count(i),activ);
        z(i) = norminv(pp);                             % z-score equivalent for p-value 
        p(i) = min(1, 2 * min(pp,1-pp));
        
        sd(i) = (task_count(i) * success(i)./ task_count(i) * (1-success(i)./ task_count(i))) .^ .5;
        %[m,v] = binostat(task_count(i),success(i)./ task_count(i));  % to get standard error
        %sd(i) = sqrt(v);                        % standard deviation of sampling dist.
    end                  

else
    
    % use Normal approximation
    meth = 'Normal';
    % observed pt given act compared to expected, which is just the pt overall


    sd = (pt .* (1-pt) ./ a) .^ .5;             % st dev of sampling dist for normal approx to binomial
    z = eff ./ sd;                              % z-score for difference between obs and expected
    p = 2*(1 - normcdf(abs(z)));                % p-values, 2-tailed

end

% additional outputs
stat.success = success;
stat.expected_s = expected_s;
stat.expected_f = expected_f;
stat.task_count = task_count;
stat.pt = pt;
stat.pa = activ;

if ~exist('sd') == 1, return, end

stat.sd = sd;
stat.meth = meth;
if doweight, stat.weighted = 'Yes';, else, stat.weighted = 'No';,end
stat.x2 = x2;
stat.x2df = x2df;
stat.x2p = x2p;

return