function [titl,tm,volinfo] = titlemarker(s)
%
% given cells with citation info (from citation_mca.m), gives title string
% and title marker index
%
% title must begin with (year)

volinfo = 'None.';

for i = 1:length(s)
    tmp = findstr(s{i},'('); if length(tmp) > 1, tmp = tmp(1); end
    tmp2 = findstr(s{i},')');if length(tmp2) > 1, tmp2 = tmp2(1); end
    if tmp == 1
    elseif s{i}(tmp - 1) == ' ' & ~isempty(tmp2 - tmp)
        % start of title
        titl = s{i};
        tm = i;
    end
end
       
if ~exist('tm','var')
    % we couldn't find a year in ( )
    disp('Could not find year in ( )')
    tm = 1;
    titl = s{tm};
end

for t = tm+1:length(s)
    tmp = findstr(s{t},'('); if length(tmp) > 1, tmp = tmp(1);, end
    tmp2 = findstr(s{t},')');if length(tmp2) > 1, tmp2 = tmp2(1);, end
    if tmp == 1
    elseif s{t}(tmp) - 1 ~= ' ' & ~isempty(tmp - tmp2)
                % we have found volume info; stop
                volinfo = cat(2,s{t:end});
    else
                titl = [titl ', ' s{t}];
    end
end

if isempty(titl), titl = 'CANNOT FIND.'; tm = length(s);, end

return
