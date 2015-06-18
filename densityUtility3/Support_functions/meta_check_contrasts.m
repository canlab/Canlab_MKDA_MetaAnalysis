function meta_check_contrasts(DB,testfield)
% meta_check_contrasts(DB,testfield)
% 
% Checks to see if each unique Contrast has one and only one level of the
% selected variable.

for i = 1:length(DB.Contrast)
    tmp = DB.(testfield)(DB.Contrast == i);

    if length(unique(tmp)) > 1, 
        fprintf(1,'Warning! Contrast %3.0f has more than one value for %s.\n',i,testfield);
    end
end




return