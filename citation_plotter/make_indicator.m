function indic = make_indicator(aus,aulistr)
% indic = make_indicator(aus,aulistr)
%
% aus is list of cell arrays with words (rows of indic are cells)
% aulistr is string matrix of words for columns of indic

indic = zeros(length(aus),size(aulistr,1));
% Make study table
for i = 1:length(aus)
    for j = 1:length(aus{i})  
        % find index of this author
        tmp = strmatch(aus{i}{j},aulistr,'exact');
        if ~isempty(tmp)
            indic(i,tmp) = 1;
        end
    end
end


return
