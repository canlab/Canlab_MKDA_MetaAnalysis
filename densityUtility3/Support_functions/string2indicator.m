function [indic, nms, condf] = string2indicator(str,varargin)
%[indic, nms, condf] = string2indicator(str,varargin)
% 
% Takes a cell vector of string labels and returns an indicator matrix
% Optional argument is a cell vector of strings for which values to match.
%
% Examples:
%[indic,nms] = string2indicator(CL{1}(1).valence);
%[indic,nms] = string2indicator(CL{1}(1).valence,{'neg' 'pos'});
%
% Note: This function changed 11/1/2023 to return names in stable order
% (the order which they appear in the original str input).
% Also, previously, the function would accept a non-cell string matrix but
% return potentially incorrect output; now enforces cell array


if ~iscell(str) 
    str = cellstr(str);
end

if ~isempty(varargin)
    nms = varargin{1};
else
    nms = unique(str, 'stable');  % tor edited 11/1/2023, 'rows' would fix edge case with non-cell string matrix, but now enforce cell anyway 
end

if ~iscell(nms) %, for i = 1:length(nms), nms2{i} = nms(i); end, nms = nms2; end  % this was returning only the first character of a string matrix
    nms = cellstr(nms);
end

nms(strcmp(nms, 'NaN')) = [];

indic = zeros(length(str), 1);

for k = 1:length(nms)
    indic(:,k) = strcmp(str, nms{k});
end

if size(nms,1) > size(nms,2)
    nms = nms';
end

if nargout > 2
% make condition function
	condf = indic2condf(indic);
end

end
