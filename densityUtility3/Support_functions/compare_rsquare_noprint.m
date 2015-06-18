function [F,p,df1,df2] = compare_rsquare_noprint(r2full,r2red,dfful,dfred)

df1 = dfred - dfful;
df2 = dfful;

num = (r2full - r2red) ./ df1;   % percentage increase per df lost

if num < 0, error('Full model error df should be smaller than reduced!');,end
    
denom = (1 - r2full) ./ df2;               % unexplained variance per df in full model

F = num ./ denom;


p = 1 - fcdf(F,df1,df2);


%fprintf(1,'R^2 change = %3.2f, F(%3.0f,%3.0f) = %3.2f, p = %3.4f\n',r2full-r2red,df1,df2,F,p);

    
return


