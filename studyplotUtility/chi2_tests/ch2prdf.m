function y = ch2prdf(x,v)
if nargin < 2, 
    error('Requires two input arguments.'); 
end

[errorcode x v] = distchk(2,x,v);

if errorcode > 0
    error('Requires non-scalar arguments to match in size.');
end

% Initialize Y to zero.
y=zeros(size(x));

y = gammpdf(x,v/2,2);

