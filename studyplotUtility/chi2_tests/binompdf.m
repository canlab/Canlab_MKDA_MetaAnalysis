function y = binompdf(x,n,p)
% BINOPDF Binomial probability density function.
%   Y = BINOPDF(X,N,P) returns the binomial probability density 
%   function with parameters N and P at the values in X.
%   Note that the density function is zero unless X is an integer.
%
%   The size of Y is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.    

%   Reference:
%      [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964, 26.1.20.

%   Copyright (c) 1993-98 by The MathWorks, Inc.
%   $Revision: 2.7 $  $Date: 1997/11/29 01:44:56 $


if nargin < 3, 
    error('Requires three input arguments');
end

[errorcode x n p] = distchk(3,x,n,p);

if errorcode > 0
    error('Requires non-scalar arguments to match in size.');
end

% Initialize Y to zero.
y = zeros(size(x));
 
% Binomial distribution is defined on positive integers less than N.
k = find(x >= 0  &  x == round(x)  &  x <= n);
if any(k),
   nk = gammaln(n(k) + 1) - gammaln(x(k) + 1) - gammaln(n(k) - x(k) + 1);
   lny(k) = nk + x(k).*log( p(k)) + (n(k) - x(k)).*log(1 - p(k));
   y(k) = exp(lny(k));
end
