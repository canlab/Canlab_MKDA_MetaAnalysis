function p = chsqcudf(x,v)
%CHI2CDF Chi-square cumulative distribution function.
%   P = CHI2CDF(X,V) returns the chi-square cumulative distribution
%   function with V degrees of freedom at the values in X.
%   The chi-square density function with V degrees of freedom,
%   is the same as a gamma density function with parameters V/2 and 2.
%
%   The size of P is the common size of X and V. A scalar input   
%   functions as a constant matrix of the same size as the other input.    

%   References:
%      [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964, 26.4.

%   Copyright (c) 1993-98 by The MathWorks, Inc.
%   $Revision: 2.6 $  $Date: 1997/11/29 01:45:04 $

if   nargin < 2, 
    error('Requires two input arguments.');
end

[errorcode x v] = distchk(2,x,v);

if errorcode > 0
    error('Requires non-scalar arguments to match in size.');
end
    
% Call the gamma distribution function. 
p = gammcdf(x,v/2,2);

% Return NaN if the degrees of freedom is not a positive integer.
k = find(v < 0  |  round(v) ~= v);
if any(k)
   tmp = NaN;
    p(k) = tmp(ones(size(k)));
end
