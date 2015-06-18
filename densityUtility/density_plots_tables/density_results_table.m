function density_results_table(OUT)
% density_results_table(OUT)
%     
% prints table of results for density_pdf or density_diff_pdf
% Tor Wager


fprintf(1,'Density results for %s\n',OUT.fname)
fprintf(1,'Brain mask file: %s with %3.2f x %3.2f x %3.2f mm voxels',OUT.mask_file,OUT.voxsize(1),OUT.voxsize(2),OUT.voxsize(3))
fprintf(1,'\nAlpha = %3.3f\n',OUT.crit_p)
fprintf(1,'\n--------------------------------------------------------------------------------\n')

fprintf(1,'Main effect results for mask 1, %3.0f points\n',OUT.n)
sub_table(OUT.u, OUT.th, OUT.cls, OUT.sz_u,OUT.cl_u)

if isfield(OUT,'cl_u2')     % check for a field made only by density_diff_pdf
    
    fprintf(1,'Main effect results for mask 2, %3.0f points\n',OUT.n2)
    sub_table(OUT.u2, OUT.th, OUT.cls, OUT.sz_u2,OUT.cl_u2)
    
    fprintf(1,'Difference effect for mask 1-2\n')
    sub_table(OUT.du, OUT.th, OUT.cls, OUT.sz_du,OUT.cl_du)
    
    fprintf(1,'Difference effect for mask 2-1\n')
    sub_table(OUT.du2, OUT.th, OUT.cls, OUT.sz_du2,OUT.cl_du2)
    
end

fprintf(1,'\nSize threshold: %3.0f%% of the time, under Ho, no cluster will exceed this number of voxels',100*(1-OUT.crit_p))
fprintf(1,'\nNum  threshold: %3.0f%% of the time, under Ho, there will be fewer than this number of \n\t\t\t\t clusters of this size',100*(1-OUT.crit_p))



function sub_table(u,th,cls,sz,numc)
% numc  : matrix of number of clusters
%         rows index density thresholds
%         columns index cluster size thresholds
%         with (th x cls) elements
%
% sz and numc are threshold values - not the original OUT.numc
% so they're the 95% random pdf value of numc and cl_size

fprintf(1,'\tCritical height threshold for any density\t%3.4f\n',u)
fprintf(1,'\tAt prespecified threshold')
for i = 1:length(th)
    fprintf(1,'\t%3.4f',th(i))
end
fprintf(1,'\n\t\t                         ---------------------------------\n')
fprintf(1,'\tSize threshold             ')

for i = 1:length(th)
    fprintf(1,'\t%3.0f\t',sz(i))
end

for j = 1:length(cls)
fprintf(1,'\n\tNum. of sz > %3.0f      \t',cls(j))
  
    for i = 1:length(th)
        fprintf(1,'\t%3.0f\t',numc(i,j))
    end
    
end

fprintf(1,'\n--------------------------------------------------------------------------------\n')

return
