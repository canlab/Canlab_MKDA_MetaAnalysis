% [chi2,df,p,sig,warn,freq_table,expected_table,isnonparametric] = chi2test(counts,datatype,[obs. weights],[nonpar flag])
%
% Weighted or unweighted Chi-square test 
% from frequency (contingency) table or rows of observations
% Optional nonparametric estimation for questionable results
%
%
% Takes either tabular (frequency) data or rows of observations
% (integer data, i.e., 1s and 0s or dummy codes for categories)
% Makes table:
% rows are No (top) then Yes (bottom)
% columns are variables in counts
%
% Datatype: is 'table' or 'obs' for tabular data or rows of observations
% if obs, then 1st column is response (1/0 for yes/no) and 2nd col. is
% condition assignment (coded by intergers)
%
% Table:
% rows are conditions
% columns are measures
% e.g., rows are students/nonstudents, columns are number reading and
% number watching TV
% OR
% rows are tasks, columns are "Yesses" and "Nos"
%
% Weights option: Weights should have a mean of one to retain
% interpretation of weighted counts
% Weights option only works for 'obs' datatype right now!
%
% % nonparametric option: if chi2 assumptions violated (any expected < 5)
% and p < .2, do nonparametric chi2
%
% WARNING:
% Weighted chi-square alpha level is not correct without nonparametric
% test; degrees of freedom adjustment due to weights is not made.
%
% E.g.
% Two columns in counts, each is either 0 or 1
% Runs chi-sq test for independence between v1 and v2
%
% Two columns in counts, one is response (0 or 1), one is dummy coded (1 2
% 3)
% Runs chi-sq test for homogeneity in the response across levels of the
% dummy-coded variable.
%
% tor wager
% tested simple examples against SPSS, Jan 23, 06
%
% [chi2,df,chi2p,sig,warn,tab,e] = chi2test([y condf],'obs'); 
% 
% Example of nonparametric chi2:
% [chi2,df,p,sig,warn,counts,e,isnonparametric] = chi2test([1 1 1 1 0 0 0 0; 1 2 1 1 2 2 2 2]', 'obs', [], 1);
%
% another example, with text print-out of results:
% [chi2,df,p,sig,warn,counts,e,isnonparametric] = chi2test([x diso], 'obs', w, 1);
% fprintf('Chi2 output: chi2: \t%3.2f\t df \t%3.2f\t p \t%3.6f\t \n', chi2, df, p);


% SIMULATION CODE
% -------------------------------------------------------------------------
% nonparametric option: if chi2 assumptions violated (any expected < 5) or
% p < .15, do nonparametric chi2

% Simulations: Test FPR
% %%
% p = .2; n = 50;
% 
% niter = 1000;
% pp = zeros(niter, 1);
% 
% sig = zeros(niter, 1);
% 
% 
% for i = 1:niter
%     data = binornd(1, p, n, 2);
%     [chi2,df,pp(i),sig(i)] = chi2test(data,'obs',[],0);
% 
% end
% 
% sum(sig) ./ niter
% 
% %%
% p = .2; n = 50;
% 
% niter = 1000;
% pp = zeros(niter, 1);
% ppols = zeros(niter, 1);
% 
% sig = zeros(niter, 1);
% sigols = zeros(niter, 1);
% 
% for i = 1:niter
%     data = binornd(1, p, n, 2);
%     w = rand(n, 1);
%     w = w ./ repmat(mean(w), n, 1);
%     
%     [chi2,df,ppols(i),sigols(i)] = chi2test(data,'obs',w,0);
%     
%     [chi2,df,pp(i),sig(i)] = chi2test(data,'obs',w,1);
% 
% end

function [chi2,df,p,sig,warn,counts,e,isnonparametric] = chi2test(counts,datatype,varargin)

w = [];
dononpar = 0;
if length(varargin) > 0
    w = varargin{1}; 
    
end      % weights for observations; optional

if length(varargin) > 1, dononpar = varargin{2}; end

% This warning cannot be done here because chi2test is called recursively
% during nonparam_chi2.m
%if ~dononpar, warning('Weighted chi-square alpha level is not correct without nonparametric test; degrees of freedom adjustment due to weights is not made.'); end
    
if isempty(w)
    % Only used for 'obs' below
    w = ones(size(counts(:,1)));
end
        
% cross-tabulate, if necessary
switch datatype
    case 'obs'

        u1 = unique(counts(:,1));
        u2 = unique(counts(:,2));
        for i = 1:length(u1)                   % rows are 0 then 1 on the first var, "Nos" then "Yesses"
            for j = 1:length(u2)               % for each column
                tab(i,j) = sum( (counts(:,1) == u1(i) & counts(:,2) == u2(j)) .* w );
            end
        end
        
        % save original data if nonparametric option is asked for
        if dononpar, y = counts(:,1); condf = counts(:,2); end
        
        counts = tab;

    case 'table'
        % do nothing,  counts already == crosstabs table
        if dononpar, error('Nonparametric option only works with datatype ''obs'''); end

    otherwise
        error('Datatype must be ''obs'' or ''table''');
end


s = sum(counts(:));

rowmarg = sum(counts,2) ./ s;   % marginal proportions for rows
colmarg = sum(counts,1) ./ s;   % marginal proportions for columns

e = rowmarg * colmarg .* s;     % expected counts

d = (counts - e).^2 ./ e;
chi2 = sum(d(:));

if nargout > 1

    warn = 0;   % warning for low expected counts (< 5)

    df = (size(counts,1)-1) .* (size(counts,2)-1);
    if any(e(:) < 5), warn = 1; end

    p = 1 - chi2cdf(chi2,df);
    sig = p < .05;

    isnonparametric = 0;
    if dononpar && p < .15    %% && warn && p < .15
        [chi2,p,sig] = nonparam_chi2(y,condf,1000,w);
        isnonparametric = 1;
    end
    
end

return