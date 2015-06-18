textlabel = {'valence'};
plots = {'surface';'medial'};

filenames = {'emotrev4_posneg_t.txt'};
plotdims= [1 2];

conditions = { ...
    % plot 1
	'strcmp(valence,''asdf'') | strcmp(valence,''fdas'')' ...
    'strcmp(colors,''yv'')' ...
    'strcmp(colors,''wo'')' ...
	'strcmp(cogdemand,''y'') | strcmp(cogdemand,''n'')' ...
    ; ...
    % plot 3  
	'(strcmp(valence,''asdf''))' ...
    '~strcmp(colors,''yv'') & ~(strcmp(region,''amy'') | strcmp(region,''hpc'') | strcmp(region,''phpc'') | strcmp(region,''bg'') | strcmp(region,''tc'') | strcmp(region,''sl''))' ...
    '~strcmp(colors,''wo'') & ~(strcmp(region,''amy'') | strcmp(region,''hpc'') | strcmp(region,''phpc'') | strcmp(region,''bg'') | strcmp(region,''tc'') | strcmp(region,''sl''))' ...
 	'(strcmp(cogdemand,''y'') | strcmp(cogdemand,''n'')) & ~(strcmp(region,''amy'') | strcmp(region,''hpc'') | strcmp(region,''phpc'') | strcmp(region,''bg'') | strcmp(region,''tc'') | strcmp(region,''sl''))' ...   
}

switch sum(plotdims)
case 2
	figpos = [123   575   934   374];
case 3
	figpos = [270    43   568   900];
case 4
	figpos = [140 53 1020 879];
otherwise
	figpos = [360   514   560   420];
end

for i = 1:size(textlabel,2)
	disp(['***************************************************************'])
	disp(['Running studyplot for ' textlabel{i}])
	disp(['***************************************************************'])
    fh = figure; set(gcf,'Position',figpos)
    plotindex = 1;
	clear h
    % ---- define figure handles ----
    for j = 1:plotdims(1) * plotdims(2)
        figure(fh); h(j) = subplot(plotdims(1),plotdims(2),j);
    end
    drawnow

    % ---- make individual plots ----
    for j = 1:size(plots,1)
        
        studyplot(plots{j},filenames{i},conditions{j,i},'handles',h(plotindex:end));
        
        switch plots{j}
        case 'surface'
			if j == 1, view(105,20),elseif j == 2,view(255,20),end
			camzoom(1.5)
            plotindex = plotindex + 1;
		case 'surface2'
			if j == 1, view(105,20),elseif j == 2,view(255,20),end
			ch = fill3([0 0 0 0],[-125 90 90 -125],[-70 -70 80 80],[.95 .84 .84]);
			camzoom(1.5)
            plotindex = plotindex + 1;
		case 'surface3'
            plotindex = plotindex + 2;
        case 'medial'   
            plotindex = plotindex + 2;
        case 'glass'   
            plotindex = plotindex + 3;
        end
     
    end
    
    % ---- save plot figure ----
    eval(['saveas(gcf,''multi_' textlabel{i} ''',''fig'');']);
    eval(['saveas(gcf,''multi_' textlabel{i} ''',''jpg'');']); close 
     
end



    %switch plots{j}
    %case 'surface'
    %    eval(['saveas(gcf,''surface_' textlabel{i} ''',''fig'');']); close 
    %case 'medial'
    %    eval(['saveas(gcf,''Rmedial_' textlabel{i} ''',''fig'');']); close 
    %    eval(['saveas(gcf,''Lmedial_' textlabel{i} ''',''fig'');']); close
    %case 'glass'
    %    eval(['saveas(gcf,''glass1_' textlabel{i} ''',''fig'');']); close 
    %    eval(['saveas(gcf,''glass2_' textlabel{i} ''',''fig'');']); close
    %    eval(['saveas(gcf,''glass3_' textlabel{i} ''',''fig'');']); close