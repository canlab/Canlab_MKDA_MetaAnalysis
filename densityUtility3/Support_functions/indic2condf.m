function condf = indic2condf(indic)
[m,n] = size(indic);
 condf = indic .* repmat(1:n,m,1);
 wh = sum(condf,2) - max(condf,[],2);
    if any(wh), warning('Task categories compared are not mutually exclusive.')
condf(wh,:) = 0; 
end
    
condf = sum(condf,2);
return

