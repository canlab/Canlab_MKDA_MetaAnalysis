function task_indicator_table(connames,taskstr,taskindic,tasknames)

% header 1
fprintf(1,'Task types in analysis\n\nContrast Name\tTask type\tTask indicators with sample sizes\n')

% header 2

fprintf(1,'\t\t')
    for i = 1:length(tasknames)
        fprintf(1,['%s\t'],tasknames{i})
    end
fprintf(1,'\n')
    
% print for for each contrast
for j = 1:size(connames,1)
    
    fprintf(1,'%s\t%s\t',connames{j},taskstr{j})
    
    str = [repmat('%3.0f\t',1,size(taskindic,2))];
    fprintf(1,str,taskindic(j,:))
    fprintf(1,'\n')
    
end

return