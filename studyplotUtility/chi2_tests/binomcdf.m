function y = binomcdf(x,n,p)
%BINOCDF Binomial cumulative distribution function.
%   Y=BINOMCDF(X,N,P) returns the binomial cumulative distribution
%   function with parameters N and P at the values in X.
%
%   The size of Y is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.    
%
%   The algorithm uses the cumulative sums of the binomial masses.

%   Reference:
%      [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964, 26.1.20.

%   B.A. Jones 1-12-93
%   Revised by ZP You 12-5-97
%   Copyright (c) 1993-98 by The MathWorks, Inc.
%   $Revision: 2.9 $  $Date: 1998/09/11 19:28:59 $

if nargin < 3, 
    error('Requires three input arguments.'); 
end 

scalarnp = (prod(size(n)) == 1 & prod(size(p)) == 1);

[errorcode x n p] = distchk(3,x,n,p);

if errorcode > 0
    error('Requires non-scalar arguments to match in size.');
end

% Initialize Y to 0.
y=zeros(size(x));

% Y = 1 if X >= N
k = find(x >= n);
y(k) = 1;


% assign 0 to p==0 indices.
k2 = (p == 0);
y(k2) = 0;

% Return NaN if any arguments are outside of their respective limits.
% this may overwrite k2 indices.
k1 = (n < 0 | p < 0 | p > 1 | round(n) ~= n | x < 0); 
y(k1) = NaN;

% Compute Y when 0 < X < N.
xx = floor(x);
k = find(xx >= 0 & xx < n & ~k1 & ~k2);
k = k(:);

% Accumulate the binomial masses up to the maximum value in X.
if any(k)
   val = max(xx(:));
   i = (0:val)';
   if scalarnp
      tmp = cumsum(binompdf(i,n(1),p(1)));
      y(k) = tmp(xx(k) + 1);
   else
      compare = i(:,ones(size(k)));
      index = xx(k);
      index = index(:);
      index = index(:,ones(size(i)))';
      nbig = n(k);
      nbig = nbig(:);
      nbig = nbig(:,ones(size(i)))';
      pbig = p(k);
      pbig = pbig(:);
      pbig = pbig(:,ones(size(i)))';
      y0 = binompdf(compare,nbig,pbig);
      indicator = find(compare > index);
      y0(indicator) = 0;
      y(k) = sum(y0,1);
   end
end

% Make sure that round-off errors never make P greater than 1.
k = find(y > 1);
y(k) = 1;
