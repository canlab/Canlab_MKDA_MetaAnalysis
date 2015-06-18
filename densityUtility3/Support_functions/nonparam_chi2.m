function [realchi2,chi2p,sig,df,tab,yp] = nonparam_chi2(y,condf,nperms,varargin)
% [realchi2,chi2p,sig,df,tab,yp] = nonparam_chi2(y,condf,nperms,[w])
%
% weights should be mean = 1

if length(varargin) > 0, w=varargin{1}; else w=ones(size(y)); end

[realchi2,df,tmp,tmp,warn,tab,e] = chi2test([y condf],'obs',w); 

n = size(y,1);

% weights
%y = y .* w;

% initialize output
chi2 = zeros(nperms,1);
yp = zeros(n,nperms);

% setup up valid permutations (no repeats), perms x observations
permindx = permute_setupperms(n,nperms); % approximate test with nperms obs.

for i = 1:nperms
    yp(:,i) = y(permindx(i,:));
    
    [chi2(i)] = chi2test([yp(:,i) condf],'obs',w);
end

chi2p = sum(chi2>=realchi2) ./ nperms;
    
chi2p = max(chi2p, 1000 * eps);

sig = chi2p <= .05;


return

