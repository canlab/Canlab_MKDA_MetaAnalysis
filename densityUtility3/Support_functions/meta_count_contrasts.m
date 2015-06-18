% [num, wh, weighted_num, conindx] = meta_count_contrasts(DB, testfield, fieldvalue)
%
% e.g., testfield = 'Stimuli'
% fieldvalue = 'faces';
% [num, wh] = meta_count_contrasts(DB, testfield, fieldvalue)

function [num, wh, weighted_num, conindx] = meta_count_contrasts(DB, testfield, fieldvalue)

    wh = strcmp(DB.(testfield), fieldvalue);
    [u, i] = unique(DB.Contrast(wh));
    num = length(u);

    conindx = zeros(length(DB.connumbers), 1);
    %conindx(i) = 1;

    if nargout > 2
        % weighted total
        for i = 1:length(u)
            tmp = find(DB.connumbers == u(i));   % the index in contrast list for this con number
            if isempty(tmp), disp('Warning! No contrast-list entry for existing contrast number.  Was DB modified?');, end
            whcons(i) = tmp;

            % save list of indices
            conindx(tmp) = 1; %DB.connumbers(tmp)) = 1;
        end

        if isempty(u), weighted_num = 0;, return, end

        w = DB.studyweight;
        w = w ./ mean(w);
        w = w(whcons);
        weighted_num = sum(w);
    end
end