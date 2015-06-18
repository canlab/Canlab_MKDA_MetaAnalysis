function [p,pobs] = fishers_exact(tab)
%
% m x n generalization of Fisher's exact test
%
% See Agresti, 1992.  Exact inference for contingency tables.  Statistical
% Science
%
% Eric W. Weisstein. "Fisher's Exact Test." From MathWorld--A Wolfram Web
% Resource. http://mathworld.wolfram.com/FishersExactTest.html
%
% Input: contingency table
% Output: p-value for Fisher's exact test
%
% This was abandoned in development because listing possibilities for more than 2 x 2 tables
% with even modest sample sizes seems computationally infeasible.
%
% see nonparam_chi2.m
% Tor Wager, May 2006

p = [];

rowsum = sum(tab,2); colsum = sum(tab,1);
N = sum(colsum);

pobs = get_observed_p(tab,N,rowsum,colsum);


% now we need to find tables whose chi2 values are at least as high as the
% observed one, and sum pobs for these tables.

%tables = list_possible_tables(tab,rowsum,colsum);





function p = get_observed_p(tab,N,rowsum,colsum)
% hypergeometric distribution given a parameter value of 1 (independence)

numer = prod([factorial(rowsum') factorial(colsum)]);

denom = prod([factorial(N) factorial(tab(:)')]);
p = numer/denom;

return



function tables = list_possible_tables(tab,rowsum,colsum);

[r,c] = size(tab);

totaln = prod(rowsum);
for i = 1:r
   % for each row, get possible values
   % num values = rowsum
   c = [0:rowsum(i); rowsum(i):-1:0]';  % combinations of this row, assuming 2 columns only!
   
   % fill in 2nd row, holding constant
   
end



return
