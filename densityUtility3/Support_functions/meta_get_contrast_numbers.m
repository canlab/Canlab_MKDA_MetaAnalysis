function DB = meta_get_contrast_numbers(DB,varargin)
% function DB = meta_get_contrast_numbers(DB,testfields)
%
% e.g., for Pain meta-analysis
% DB = meta_get_contrast_numbers(DB,'Study','Subjects','Kind_of_Pain','Limb_v_Head_or_Chest','Right_vs_Left');

% concatenate fields of interest
% ------------------------------------
for i=1:length(varargin)
    if i == 1, tmp = DB.(varargin{i});,
    elseif length(DB.(varargin{i})) == length(tmp)
        
        myfield = DB.(varargin{i});
        if iscell(myfield)
            tmp = [tmp myfield];
        else
            % convert to cell if it's not
            for i = 1:length(myfield)
                myfield2{i,1} = num2str(myfield(i));
            end
            myfield = myfield2;
            tmp = [tmp myfield];
        end
    else
        error('testfields are not all same length.');
    end
end

clear tmp2
for i=1:size(tmp,1),
    tmp2{i,1} = cat(2,tmp{i,:});
end

% concatenate fields of interest
% ------------------------------------

%str = ['[levels,i3,j] = unique(DB.' testfield ');']; eval(str)

[levels,i3,j] = unique(tmp2);

% preserve order
[i2,wh]=sort(i3); levels=levels(wh);

contrast_numbers = NaN .* zeros(size(tmp2,1),1);

for i = 1:length(levels)
    
    wh = find(strcmp(tmp2,levels{i}));
    
    contrast_numbers(wh) = i;
    
end

DB.Contrast = contrast_numbers;

return
