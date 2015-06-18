function [OUT] = bootstrap_mtx(data,e,varargin)
% [OUT] = bootstrap_mtx(data,e,[meth],[verbose],[niter])
%
% tor wager
%
% BOOTSTRAP TEST FOR SIG. DIFFERENCE IN COVARIANCE FROM e (Ho)
%
% This function takes as input two matrices,
% an actual matrix of data
% and an expected / null hypothesis
% covariance or correlation matrix.
%
% the algorithm provides significance levels for
% the cov / correlation of the columns of data,
% based on the expected correlation e
%
% it computes a test statistic, which is the
% average squared deviation from the expected
% values, element by element, over the matrix.
% This statistic provides an "omnibus" test
% of whether there are deviations from the 
% expected values in the matrix.
%
% for each element, a statistic is also computed -
% the squared deviation from the expected -
% which provides a test of significance for each
% element.  
%
% a permutation test is used to create an Ho 
% distribution.  Columns of the matrix are
% randomly permuted (independently), and the test
% statistics assessed over iterations.
%
% called in db_cluster_burt_table.m
%
% example:
% data = bmatx1(:,1:5); ss = data'*data; e = diag(diag(ss));
% OUT = permute_mtx(bmatx1(:,1:5),e,'ss',1,1000);

if length(varargin) > 0, meth = varargin{1};, else, meth = 'corr';, end
if length(varargin) > 1, verbose = varargin{2};, else, verbose = 1;, end
if length(varargin) > 2, niter = varargin{3};, else, niter = 5000;, end

switch meth
case 'cov', m = cov(data);
case 'ss', m = data'*data;
case 'corr', m = corrcoef(data);
otherwise, error('Unknown method!  OK methods are cov, ss, and corr')
end

if verbose,
    fprintf(1,'\nBootstrap_mtx.m\n--------------------------\n')
    fprintf(1,'Test that columns of data have no systematic relationship\n')
    fprintf(1,'\nActual %s matrix for data',meth),m,
    fprintf(1,'\nExpected Ho %s matrix for data',meth),e,
end

% get index of lower triangular matrix
% to compute stats only on these values.
mask = triu(Inf*eye(size(m,2))); maskind = find(mask==0);


[s1,s2] = getstats(m,e,maskind);

for i = 1:niter
    
    if mod(i,1000) == 1, fprintf(1,'.');, end
    [mi] = getcov(data,meth);
    [s1n(i),s2n(i,:)] = getstats(mi,e,maskind);
    meanmi(:,:,i) = mi;
    
end
fprintf(1,'\n')
meanmi = mean(meanmi,3);

% p is prob that test statistic lies below 0
p = 1 - sum(s1n < 0) ./ niter;
for i = 1:length(s2), p2(i) = 1 - sum(s2n(:,i) < 0)./ niter;, end
p2c = p2 .* length(p2);

if verbose,
    fprintf(1,'\nMean Ho %s matrix\n',meth),meanmi
    fprintf(1,'\nOmnibus test for differences from expected on %s\n',meth)
    fprintf(1,'D2 %3.3f, Ho mean D2 = %3.3f, ste = %3.3f, p = %3.4f\n',s1,mean(s1n),std(s1n),p)
    fprintf(1,'\nSignificant individual tests for differences from expected on %s\n',meth)
    sig = find(p2 <= .05);
    for i = 1:length(sig)
        [row,col] = ind2sub(size(m),maskind(sig(i)));
        fprintf(1,'[%3.0f,%3.0f], D2 = %3.3f, Ho D2 = %3.3f, ste = %3.3f, p = %3.4f, bonf_p = %3.4f\n',row,col,s2(sig(i)),nanmean(s2n(:,sig(i))),nanstd(s2n(:,sig(i))),p2(sig(i)),p2c(sig(i)))
        %figure; hist(s2n(:,sig(i)))    
    end
    if isempty(sig), disp('No significant results.'), end
end

OUT.m = m; OUT.e = e; OUT.s1 = s1; OUT.s2 = s2;
OUT.p = p; OUT.p2 = p2; OUT.p2c = p2c;OUT.meanmi = meanmi;
OUT.sig = sig; OUT.niter = niter; OUT.meth = meth; OUT.s1n = s1n; OUT.s2n = s2n;

return




function [m] = getcov(data,meth)

% this would be a bootstrap
% resample data - tested to evenly sample rows
len = size(data,1);
wh = ceil(rand(len,1) .* len);
d2 = data(wh,:);

% permute columns
%for i = 1:size(data,2)
%    d2(:,i) = getRandom(data(:,i));
%end
    
switch meth
case 'cov', m = cov(d2);
case 'ss', m = d2'*d2;
case 'corr', m = corrcoef(d2); 
otherwise, error('Unknown method!  OK methods are cov, ss, and corr')
end

return



function [s1,s2] = getstats(m,e,maskind)
% input is cov matrix of whatever form - corr, ss, cov

m = m(maskind); m = m(:);
e = e(maskind); e = e(:);

s2 = ((m - e) .^ 2)';
s1 = nanmean(s2);

return


