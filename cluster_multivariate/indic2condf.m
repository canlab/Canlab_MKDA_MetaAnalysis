function condf = indic2condf(indic)
% condf = indic2condf(indic)
%
% indic is an m x n matrix valued [1, 0]
% Entries are designed to index which of n categories each element of m is
% a member of, i.e., n columns are n sets.
% 
% indic2condf converts these into an integer-valued condition function 
% of m x 1 elements, with integers indicating which set each element belongs to 
% i.e., which of the n columns.
% 
% If indic has non-zero elements for multiple columns, this function
% returns a warning and values of zero.
%
% Written by tor wager, ages ago...

[m,n] = size(indic);
condf = indic .* repmat(1:n, m, 1);

wh = sum(condf,2) - max(condf,[],2);

if any(wh), warning('Task categories compared are not mutually exclusive.')
    condf(wh,:) = 0;
end

condf = sum(condf,2);

end 

