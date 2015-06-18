function plot_all_surface(testregion,coordfile,testcondition,region,namecode,donumbers,varargin)
% plot_all_surface(testregion,coordfile,cdef,testcondition,region,namecode,donumbers,[opt] cdef)
%
% Here's how this works:
% You need to make a text file that has fields for x, y, and z elements, and a column to define 
% colors (cdef).  check read_database.m for text file setup.
% try 'all' for testregion.  See also studyplot.m, called by this function.
% studyplot.m can be used to generate individual plots.
%
% testregion    a cell array, whose n entries create n multiplots for separate regions
%               if 'all', this function creates 1 multiplot (4 panels) with all regions.
%
% coordfile     the name of the text file to load
% cdef          OPTIONAL: text string, name of column in database to define colors by
%               if colors is not a column in the database.
%
% namecode      the name of the output code to save files under
% donumbers     plot study numbers with studyplot (1) or just points (0)
%               right now, numbers only works for glass, so this won't do anything
% region        the cell array of region values, or empty to use those in the file
%
% examples of how some variables are organized:
%
% plot_all_surface({'all'},{'emotrev7_emotion2_t.txt'},'colors',{'~strcmp(colors,''wo'') & ~(strcmp(region,''amy'') | strcmp(region,''hpc'') | %strcmp(region,''phpc'') | strcmp(region,''bg'') | strcmp(region,''tc'') | strcmp(region,''sl''))'},[],'newemomulti',0)
%
% textlabel = {'valence2','method2','emotion2','cogdemand2'};
% plots = {'surface';'surface2'};
% filenames = {'emotrev4_posneg_t.txt' 'emotrev5_stimtype_t.txt' 'emotrev6_emotion_t.txt' 'emotrev9_cog_t.txt'};
%conditions = { ...
% plot 1
%	'(strcmp(valence,''pos'') | strcmp(valence,''neg''))' ...
%    '~strcmp(colors,''yv'')' ...
%    '~strcmp(colors,''wo'')' ...
%	'(strcmp(cogdemand,''y'') | strcmp(cogdemand,''n''))' ...
%    ; ...
%    % plot 2  
%	'(strcmp(valence,''pos'') | strcmp(valence,''neg'')) & ~(strcmp(region,''amy'') | strcmp(region,''hpc'') | strcmp(region,''phpc'') | strcmp(region,''bg'') | strcmp(region,''tc'') | strcmp(region,''sl''))' ...
%    '~strcmp(colors,''yv'') & ~(strcmp(region,''amy'') | strcmp(region,''hpc'') | strcmp(region,''phpc'') | strcmp(region,''bg'') | strcmp(region,''tc'') | strcmp(region,''sl''))' ...
%    '~strcmp(colors,''wo'') & ~(strcmp(region,''amy'') | strcmp(region,''hpc'') | strcmp(region,''phpc'') | strcmp(region,''bg'') | strcmp(region,''tc'') | strcmp(region,''sl''))' ...
% 	'(strcmp(cogdemand,''y'') | strcmp(cogdemand,''n'')) & ~(strcmp(region,''amy'') | strcmp(region,''hpc'') | strcmp(region,''phpc'') | strcmp(region,''bg'') | strcmp(region,''tc'') | strcmp(region,''sl''))' ...   
%
%   ; ...
%   plot 3    
%   '~strcmp(colors,''yv'') & ~(strcmp(region,''amy'') | strcmp(region,''hpc'') | strcmp(region,''phpc'') | strcmp(region,''bg'') | strcmp(region,''tc'') | strcmp(region,''sl''))' ...
%   '~strcmp(colors,''wo'') & ~(strcmp(region,''amy'') | strcmp(region,''hpc'') | strcmp(region,''phpc'') | strcmp(region,''bg'') | strcmp(region,''tc'') | strcmp(region,''sl''))' ...
%}

if length(varargin) > 0, cdef = varargin{1};, end

plotdims= [2 2];
plots = {'surface';'surface';'medial'};

if strcmp(testregion,'all')
    % set up to create one multiplot 
    % -------------------------------------------------
    textlabel = []; textlabel{1} = [namecode '_all'];
    for j = 1:length(plots)
        conditions{j,1} = testcondition;
    end
    filenames{1} = coordfile;
else    
    % set up to create many multiplots for each region
    % -------------------------------------------------
    for i = 1:length(testregion)
        for j = 1:length(plots)
            textlabel{j,i} = [namecode '_' testregion{i}];
            conditions{j,i} = [testcondition ' & strcmp(region,''' testregion{i} '''))']; 
        end
        filenames{i} = coordfile;
    end
end

switch plotdims(1)*plotdims(2)
case 2
	figpos = [123   575   934   374];
case 3
	figpos = [140    53   400   879];  %[270    43   568   900];
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
        
        if isempty(region)
            if donumbers
                if exist('cdef') == 1
                    studyplot(plots{j},filenames{i},'select',conditions{j,i},'handles',h(plotindex:end),'numbers','cdef',cdef);
                else
                    studyplot(plots{j},filenames{i},'select',conditions{j,i},'handles',h(plotindex:end),'numbers');
                end
               
            else
                if exist('cdef') == 1
                    studyplot(plots{j},filenames{i},'select',conditions{j,i},'handles',h(plotindex:end),'cdef',cdef);
                else
                    studyplot(plots{j},filenames{i},'select',conditions{j,i},'handles',h(plotindex:end));
                end
            end
        else
            if donumbers
                if exist('cdef') == 1
                    studyplot(plots{j},filenames{i},'select',conditions{j,i},'handles',h(plotindex:end),'newregion',region,'numbers','cdef',cdef);
                else
                    studyplot(plots{j},filenames{i},'select',conditions{j,i},'handles',h(plotindex:end),'newregion',region,'numbers');
                end
                
            else
                if exist('cdef') == 1
                    studyplot(plots{j},filenames{i},'select',conditions{j,i},'handles',h(plotindex:end),'newregion',region,'cdef',cdef);
                else
                    studyplot(plots{j},filenames{i},'select',conditions{j,i},'handles',h(plotindex:end),'newregion',region);
                end
            end
        end
        
        switch plots{j}
        case 'surface'
            plotindex = plotindex + 1;
			view(105,20)
			camzoom(1.5)
            material dull
            
            if j == 2,  % turn brain to the other side
                view(255,20)
                lh = lightangle(300,20); set(lh,'Color',[.5 .5 .5])
            end
            
            
		case 'surface2'
			view(105,20)
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
    eval(['saveas(gcf,''surf_' textlabel{i} ''',''fig'');']);
    eval(['saveas(gcf,''surf_' textlabel{i} ''',''jpg'');']); close 
     
end


return

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