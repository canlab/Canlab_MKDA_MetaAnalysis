function [contrasts, first_peak_in_each_con] = get_contrast_indicator_improved(DB, myfields)
    % [contrasts, first_peak_in_each_con] = get_contrast_indicator_improved(DB, myfields)
    %
    % [contrasts, first_peak_in_each_con] = get_contrast_indicator_improved(LTMDB, {'Memtype' 'Control'});
    %
    % Simple function for counting contrasts, which are assumed to be unique
    % values across a set of variables
    % In meta-analysis, coordinate entries are peak coordinates
    % If peaks come from the same study and they have the same values across
    % all levels of all factors specified in myfields (a cell array of field
    % names), then they are counted as unique contrasts.

    [indic,nms,condf] = string2indicator(DB.Study);

    for i = 1:length(myfields)

        [indic2,nms2,condf2] = string2indicator(DB.(myfields{i}));

        condf = [condf condf2];

    end

    % unique integers for each combination of levels in each study
    [uniquecons, first_peak_in_each_con] = unique(condf, 'rows', 'stable');

    n = size(condf,1); % number of rows

    contrasts = zeros(n, 1);

    for i = 1:size(uniquecons,1)
        wh = ~any(condf - repmat(uniquecons(i,:), n, 1), 2);

        contrasts(wh) = i;
    end


end


