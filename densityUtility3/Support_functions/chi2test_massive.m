% [chi2,df,p,sig,warn,isnonparametric] = chi2test_massive(seedcounts, counts,[obs. weights],[nonpar flag])
%
% Weighted or unweighted Chi-square test
% from observations
% Optional nonparametric estimation for questionable results
% SEE CHI2TEST.M
% THE PURPOSE OF THIS FUNCTION IS TO IMPLEMENT A CHI2 TEST EFFICIENTLY WITH
% A MASSIVE NUMBER OF OUTCOME VARIABLES
%
% Inputs:
% seedcounts = a vector of 1's and 0's, N obs x 1, for yes and no counts from seed variable.
% counts = a vector of 1's and 0's for N obs x k variables (e.g., brain voxels)
%
% Takes rows of observations
% (integer data, i.e., 1s and 0s or dummy codes for categories)
% Makes table:
% rows are No (top) then Yes (bottom) in "seed counts" data vector
% columns are variables in counts
%
% Weights option: Weights should have a mean of one to retain
% interpretation of weighted counts
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
%
% Examples:
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

function [chi2,df,p,sig,warn,isnonparametric] = chi2test_massive(seedcounts, counts,varargin)

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

disp('Checking data')
n = size(seedcounts, 1);
if n < length(seedcounts), error('seed variable must be column vector'); end

if length(unique(seedcounts)) > 2, error('Seed counts must be 1 or 0 values'); end
if length(unique(counts(:))) > 2, error('Data matrix counts must be 1 or 0 values'); end

% data
% rows are 0 then 1 on the first var, "Nos" then "Yesses"

disp('Getting counts')

mydat = counts(seedcounts == 0, :);
nono = sum(mydat == 0);
noyes = sum(mydat == 1);

mydat = counts(seedcounts == 1, :);
yesno = sum(mydat == 0);
yesyes = sum(mydat == 1);

%totalcount = nono + noyes + yesno + yesyes; = n

% expectations
disp('Calculating expectations')

seedpropno = 1 - (sum(seedcounts) ./ n);
seedpropyes = sum(seedcounts) ./ n;

propno = sum(counts == 0) ./ n;
propyes = sum(counts == 1) ./ n;

expnono = seedpropno .* propno .* n; %expected value for all 4 cells of indepndence chi-sq
expnoyes = seedpropno .* propyes .* n;

expyesno = seedpropyes .* propno .* n;
expyesyes = seedpropyes .* propyes .* n;

% chi-square
disp('Calculating chi-square stats')

chinono = (nono - expnono) .^ 2 ./ expnono;
chinoyes = (noyes - expnoyes) .^ 2 ./ expnoyes;
chiyesno = (yesno - expyesno) .^ 2 ./ expyesno;
chiyesyes = (yesyes - expyesyes) .^ 2 ./ expyesyes;

chi2 = sum([chinono; chinoyes; chiyesno; chiyesyes]);
chi2 = full(chi2);

if nargout > 1
        
    df = ones(1, size(counts, 2));
    
    warn = full(any([expnono; expnoyes; expyesno; expyesyes] < 5));
    
    p = 1 - chi2cdf(chi2, df);
    sig = p < .05;
    
    isnonparametric = zeros(1, size(counts, 2));
    badchi2 = warn & p < .15;
    
    fprintf('Variables with close to sig chi-square that violate assumptions: %3.0f\n', sum(badchi2));
    
    if dononpar && any(badchi2)
        
        fprintf('Calculating nonparametric values for variables that violate assumptions: 00000\n');
        wh = find(badchi2);
        
        for i = 1:length(wh)
            
            v = badchi2(i);
            
            [chi2(v),p(v),sig(v)] = nonparam_chi2(y,condf,1000,w);
            isnonparametric(v) = 1;
            
            fprintf('\b\b\b\b\b%3.5f', i);
            
        end
        
        fprintf('\n')
    end
    
end % if nargout > 1

return