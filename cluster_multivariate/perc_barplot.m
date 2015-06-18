function out = perc_barplot(studycount,totalstudycount,code,varnamecode,regionnames)
% out = perc_barplot(studycount,totalstudycount,code,varnamecode,regionnames)
%
% given counts, does barplot and CHI2 analysis
% easy to use with output from dbcluster_contrast_table.m
%
% perc :    counts or percent studies to plot; row is region, column is
%           condition
% code :    cell string of legend for conditions

colors = {'bo' 'gs' 'r^' 'yv' 'cs' 'mo'};

format compact

% get things in the right transposition

studycount = studycount';
if length(regionnames) > size(regionnames,1), regionnames = regionnames';,end
if length(totalstudycount) > size(totalstudycount,2), totalstudycount = totalstudycount';,end

perc = 100 * studycount ./ repmat(totalstudycount,size(studycount,1),1);


eval(['diary ' varnamecode '_chisq.out']); diary off      %   varnamecode,testregion, diary off
mytitle = [];
myylabel = ['Percentage of studies'];
mylegend = code;         %{'Spatial','Object','Verbal'};

for i = 1:length(mylegend), mylegend{i}(mylegend{i} == '_') = ' ';,end

disp('--- --- --- --- --- --- --- --- ---')

% make color map for bar graphs
% ======================================================
figure;
for i = 1:length(colors)
    h = plot(0,0,colors{i});
    cmap(i,:) =  get(h,'Color');
end
cmap = cmap(1:length(code),:);
close


fig2 = figure('Position',[122.0000  610.0000  921.0000  369.0000]);
set(gca,'FontSize',18)
myplot2 = subplot(1,1,1);
plotbarcounts(perc,1,regionnames,'title',mytitle,'ylabel',myylabel,'legend',mylegend,'plothandle',myplot2,'overallcounts',totalstudycount,'arealabels',0);
xlabel('')
%ylabel('% of total studies','FontSize',16)
colormap(cmap)
drawnow


drawnow


% --------------------------------------------------------------------------------
% chi2 computation - effects of CONDITION
% --------------------------------------------------------------------------------

[chi2,actual,expected] = computeynchi2(studycount,totalstudycount);
%[chi2pts,actualpts,expectedpts] = computeynchi2(count,totalcount);

diary on

disp(['Analysis for: ' varnamecode])
disp('--------------------------------------------------------')
disp('CHI2 RESULTS: STUDY COUNT FOR effects of condition')
disp('--------------------------------------------------------')
fprintf(1,'region\tchi2\tdf\tp\tcramV\tsig\twarn\tpyes\n')
disp('--------------------------------------------------------')
for trow = 1:size(chi2,1)
    if chi2(trow,3) <= .05, mystr='*';, elseif chi2(trow,3) <= .10, mystr='+';, else mystr = ' ';, end
    fprintf(1,'%s\t%3.2f\t%3.0f\t%3.3f\t%3.3f\t%3.0f\t%3.0f\t%3.2f\t%s\n',regionnames{trow},chi2(trow,1),chi2(trow,2),chi2(trow,3),chi2(trow,4),chi2(trow,5),chi2(trow,6),chi2(trow,7),mystr)
end

%disp('--------------------------------------------------------')
%disp('CHI2 RESULTS: POINT COUNT FOR effects of condition')
%disp('--------------------------------------------------------')
%fprintf(1,'region\tchi2\tdf\tp\tcramV\tsig\twarn\tpyes\n')
%disp('--------------------------------------------------------')
%for trow = 1:size(chi2pts,1)
%    if chi2pts(trow,3) <= .05, mystr='*';, elseif chi2pts(trow,3) <= .10, mystr='+';, else mystr = ' ';, end
%    fprintf(1,'%s\t%3.2f\t%3.0f\t%3.3f\t%3.3f\t%3.0f\t%3.0f\t%3.2f\t%s\n',testregion{trow},chi2pts(trow,1),chi2pts(trow,2),chi2pts(trow,3),chi2pts(trow,4),chi2pts(trow,5),chi2pts(trow,6),chi2pts(trow,7),mystr)
%end
%diary off


% ------------------------------------------------
% * this is for making stars on the graph
% ------------------------------------------------
figure(fig2), subplot(1,1,1);
myssig{1} = find(chi2(:,5)); mysplus{1} = find(chi2(:,3) <= .1 & ~chi2(:,5));

%mysig{1} = find(chi2pts(:,5)); myplus = find(chi2pts(:,3) <= .01);
myperc = max(perc')';

shift = 0; % [-.35,-.15,0,.2];
ind =1;
if ~isempty(myssig)
for analy = 1:1
    for i = 1:length(myssig{analy})
        h(ind) = text(myssig{analy}(i) + shift(analy),myperc(myssig{analy}(i))+sign(myperc(myssig{analy}(i)))*.07*max(max(myperc)),'*','FontSize',24,'FontWeight','b');
        ind = ind+1;    
    end
    
    for i = 1:length(mysplus{analy})
            hh(ind) = text(mysplus{analy}(i) + shift(analy),myperc(mysplus{analy}(i))+sign(myperc(mysplus{analy}(i)))*.07*max(max(myperc)),'+','FontSize',18,'FontWeight','b');
            ind = ind+1; 
    end        
end
drawnow
end

%subplot(2,1,2);
%myperc = max(percentcount')';
%shift = 0; % [-.35,-.15,0,.2];
%ind =1;
%if ~isempty(mysig)
%for analy = 1:1
%    for i = 1:length(mysig{analy})
%        h(ind) = text(mysig{analy}(i) + shift(analy),myperc(mysig{analy}(i))+sign(myperc(mysig{analy}(i)))*.10*max(max(myperc)),'*','FontSize',14,'FontWeight','b');
%        ind = ind+1;    
%    end
%end
%drawnow
%end

% keyboard

saveas(gcf,[varnamecode '_bar'],'fig')
saveas(gcf,[varnamecode '_bar'],'tif')

out.perc = perc;
out.chi2 = chi2;
out.myssig = myssig;
tmp = diff(out.perc'); tmp = tmp(:,out.myssig{1});
if size(tmp,1) > 1, 
    % we have more than 2 conditions
    % pick the one w/the greatest response
    tmp = perc' == repmat(max(perc'),size(perc,2),1);
    [tmp,dummy] = find(tmp);
    out.Greatest = tmp;
    out.GreatNames = 'vector of wh conditions show largest %; indices of regions are in out.myssig{1}';
    out.AvsB = [];
    out.BvsA = [];
else
    out.AvsB = out.myssig{1}(tmp<0,1);  % index of regions with sig > A than B
    out.BvsA = out.myssig{1}(tmp>0,1);
    out.AB = 'indices of sig. clusters, valid for 2-condition comparisons only';
    out.Greatest = [];
end

return




% EXTRA STUFF

% --- plot ---
	
% no longer plotting laterality of studycount
%figure(fig2)
%myylabel = ['R <- L-R studies -> L'];
%mytitle = ['Lateralization of regional comparisons'];
%myplot2 = subplot(2,1,2);
%plotbarcounts(diffstudies,1,'title',mytitle,'ylabel',myylabel,'legend',mylegend,'colors',colors,'plothandle',myplot2,'overallcounts',totalstudycount,'arealabels',0);
%legend off
%colormap(cmap)
%saveas(gcf,[varnamecode '_bar3'],'fig')
%saveas(gcf,[varnamecode '_bar3'],'jpg')



%figure
%myylabel = ['L <-       lateralization of peaks       -> R'];
%mytitle = ['Lateralization of regional comparisons'];
%plotbarcounts(diff,1,'title',mytitle,'ylabel',myylabel,'legend',mylegend,'colors',colors,'overallcounts',totalcount,'arealabels',0);
%legend off
%colormap(cmap)
%view(90,90)
%set(gcf,'Position',[495   197   654   776])
%saveas(gcf,[varnamecode '_bar4'],'fig')
%saveas(gcf,[varnamecode '_bar4'],'jpg')
warning on

dostuff = 0;
if dostuff

% for plotting glass brains with each area - pos and neg points
% taken from evolvedrunstudyplot3
% ==========================================================================================
plots = {'glass'};
plotdims= [3 1];
textlabel{1} = [];
%conditions{1} = [];

for i = 1:length(testregion)
    textlabel{i} = [oldvarnamecode '_' testregion{i}];
    % put this at the top!
    %conditions{i} = ['((strcmp(valence,''pos'') | strcmp(valence,''neg'')) & strcmp(region,''' testregion{i} '''))']; 
    filenames{i} = coordfile;
end

switch plotdims(1)*plotdims(2)
case 2
	figpos = [123   575   934   374];
case 3
	figpos = [140    53   400   879];  %[270    43   568   900];
case 4
	figpos = [140 53 1020 879];
otherwise
	figpos = [123   575   934   374];
end

% cd ..   % get up to main TICS directory where database is stored

if doglass
    diary on
    
for i = 1:length(textlabel)
	disp(['***************************************************************'])
	disp(['Running studyplot for ' textlabel{i} ' / conditions are: ' conditions{i}])
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
    for j = 1:length(plots)
        
        if donumbers
            studyplot(plots{j},filenames{i},'cdef',mycondition,'select',conditions{j,i},'handles',h(plotindex:end),'newregion',region,'numbers');
        else
            studyplot(plots{j},filenames{i},'cdef',mycondition,'select',conditions{j,i},'handles',h(plotindex:end),'newregion',region);
        end
    
        switch plots{j}
        case 'surface'
            plotindex = plotindex + 1;
		case 'surface2'
			if j == 1, view(105,20),elseif j == 2,view(255,20),end
            plotindex = plotindex + 1;
		case 'surface3'
            plotindex = plotindex + 2;
        case 'medial'   
            plotindex = plotindex + 2;
        case 'glass'   
            plotindex = plotindex + 3;
            subplot 311
            material dull, camlight left
            subplot 312
            material dull
            subplot 313
            material dull
        end
     
    end

    %saveas(gcf,[oldvarnamecode '_' testregion{i}],'jpg')
    %title([oldvarnamecode '_' testregion{i}],'FontSize',14,'Color','r')
    %saveas(gcf,[oldvarnamecode '_' testregion{i}],'fig')
    saveas(gcf,[oldvarnamecode '_' conditions{i}],'jpg')
    title([oldvarnamecode '_' conditions{i}],'FontSize',14,'Color','r')
    saveas(gcf,[oldvarnamecode '_' conditions{i}],'fig')
    %pause(5)
    %close
end   

end % end if doglass

% dosurface
% multiplot of surface
% ==========================================================================================
if dosurface
    eval(['cdef = ' mycondition ';'])
    % make sure mycondition is global!
    plot_all_surface('all',coordfile,conditions{1},region,oldvarnamecode,donumbers,mycondition)
    %plot_all_surface(testregion,coordfile,testcondition,region,namecode,donumbers,varargin)
end

end
