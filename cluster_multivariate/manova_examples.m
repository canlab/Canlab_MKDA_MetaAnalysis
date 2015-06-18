% ------------------------------------------------
% Crossing data in multivariate space
% ------------------------------------------------

% make 2 vars with positive relationship
C = [1 .7; .7 1];
x = randn(100,2);

x2 = x * C;

% add another group with negative relationship
C2 = [1 -.7; -.7 1];
x2 = [x2; x * C2];

% make grouping variable
group = ones(200,1);  group(101:end) = 2;

% plot the data
figure;plot(x2(group==1,1),x2(group==1,2),'bo')
hold on; plot(x2(group==2,1),x2(group==2,2),'rs')

[d,p,stats] = manova1(x2,group)
p = 1 - chi2cdf(stats.chisq,stats.chisqdf)

% from C:\tor_scripts\Matlab_toolboxes\RegrToolbox
[F,error,R2] = mlr([group==1 ones(size(group))],x2); F, R2

% ------------------------------------------------
% Classic ANCOVA situation
% ------------------------------------------------

% make 2 vars with positive relationship
C = [1 .7; .7 1];
x = randn(100,2);
x2 = x * C;

% add another group with shift
x = randn(100,2); 
x = x * C; x(:,2) = x(:,2) + .5;
x2 = [x2; x];

% make grouping variable
group = ones(200,1);  group(101:end) = 2;

% plot the data
figure;plot(x2(group==1,1),x2(group==1,2),'bo')
hold on; plot(x2(group==2,1),x2(group==2,2),'rs')

[d,p,stats] = manova1(x2,group)
p = 1 - chi2cdf(stats.chisq,stats.chisqdf)

