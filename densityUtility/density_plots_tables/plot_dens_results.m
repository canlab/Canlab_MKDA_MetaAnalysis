function [f1,f2,OUT] = plot_dens_results(OUT)
% function [f1,f2,OUT] = plot_dens_results(OUT)
%
% ----------------------------------------------------------------
% * read input fields into variables
% ----------------------------------------------------------------
N = fieldnames(OUT);
for i = 1:length(N)
    eval([N{i} ' = OUT.' N{i} ';'])
end


% ----------------------------------------------------------------
% * read max density results if they don't exist
% ----------------------------------------------------------------
if ~isfield(OUT,'xpdf')
    [u,xpdf,ypdf,xcdf,ycdf,h] = d_pdf(maxd,crit_p);
end

% ----------------------------------------------------------------
% * plot max density results
% ----------------------------------------------------------------
f1 = figure;f2 = figure;

figure(f1)
bar(xpdf,ypdf)
mylegend{1} = 'max density';
title('Probability density function of maximum point density under the null hypothesis')

figure(f2); 
bar(xcdf,ycdf)
mylegend{1} = 'max density';
title('Cumulative density function of maximum point density under the null hypothesis')

% ----------------------------------------------------------------
% * plot cluster size
% ----------------------------------------------------------------

f1 = figure; f2 = figure;
index = 0;
for j = 1:size(numc,2)
    for i = 1:size(numc,1)
        mycolor = rand(1,3);
        
        index = index + 1;
        vec = str2num(num2str(numc(i,j,:)));
        [cl_u(index),d,d,d,d,h] = d_pdf(vec,crit_p,2,f1,f2);
        
        mylegend{index} = ['u = ' num2str(th(i)) ' k = ' num2str(cls(j))];
        set(h,'Color',mycolor)
        hold on;
        plot(cl_u(index),.95,'d','Color',mycolor)

    end
end
cl_u = reshape(cl_u,size(numc,1),size(numc,2));

figure(f1); legend(mylegend);
title('Probability density function of expected number of clusters under the null hypothesis')
figure(f2); legend(mylegend)
title('Cumulative density function of expected number of clusters under the null hypothesis')

OUT.legend = mylegend;

return
