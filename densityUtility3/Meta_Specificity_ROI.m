function Meta_Specificity_ROI(XYZlist, varargin)
    % Meta_Specificity_ROI(XYZlist, [overlay])
    %XYZlist = [36 -18 3] %[39 -15 21] % [34 -14 19]; % dpins %[60 -31 25] SII/op1
    %
    overlay =which('spm2_single_subj_T1_scalped.img');

    for i = 1:length(varargin)
        if ischar(varargin{i})
            switch varargin{i}

                % functional commands
                case 'overlay', overlay = varargin{i+1};
               
                otherwise, warning(['Unknown input string option:' varargin{i}]);
            end
        end
    end
    
    
    cluster_orthviews();
spm_check_registration(overlay);

spm_orthviews('reposition', XYZlist)
sphere_roi_tool_2008([], 10, 'useexisting');
%%

basedir = pwd;
dirname = fullfile(basedir, sprintf('Meta_Specificity_ROI_output_%3.0f_%3.0f_%3.0f', XYZlist(1, :)));
mkdir(dirname); cd(dirname)

filename = fullfile(dirname, sprintf('Meta_Specificity_ROI_output_%3.0f_%3.0f_%3.0f.txt', XYZlist(1, :)));

diary(filename)
OUT = meta_analysis_table(XYZlist, 10); 

%% Figure

OUT.dbname, 
whdb = [2 7 3 4 8 10], 
names = OUT.dbname(whdb)

names = {'Emotion' 'Pain' 'Inh' 'Mem' 'Switch' 'WM'};    
p = OUT.percentage_by_database(whdb);
create_figure('bar'); h = bar(p);
set(gca, 'Xtick', 1:length(whdb), 'XTickLabel', names, 'FontSize', 24)
set(h, 'Facecolor', [.5 .5 .5])

ylabel('% of studies activating')
p = p ./ 100;
sebars = 100 * (p .* (1 - p) ./ OUT.nstudies(whdb)) .^ .5;

p = p*100;
p(p == 0) = .5;
tor_bar_steplot(p, sebars, {'k'});

for i = 1:length(whdb)
    str{i} = sprintf('%1.0f/%3.0f\nstudies', OUT.nyes(whdb(i)), OUT.nstudies(whdb(i)));
    text(i-.5, p(i)+ 5, str{i}, 'FontSize', 24);
end

%% Bayesian stuff

p = OUT.percentage_by_database ./ 100; 
n = OUT.nstudies;
p = p(whdb); 
n = n(whdb);

pa = sum(p .* n) ./ sum(n); % prior on activation
ps = n ./ sum(n)            % prior on study type
L = p;                      % likelihood of act given study type
                            % Bayes' rule
pp = L .* ps ./ pa;         % Posterior prob, p(study type | act), p(ps | act)

% assume flat priors on study type
n = ones(1, length(n));
pa = sum(p .* n) ./ sum(n); % prior on activation
ps = n ./ sum(n);            % prior on study type

pp = L .* ps ./ pa;         % Posterior prob, p(study type | act), p(ps | act)

% sensitivity, prob act given task category
sens = p(2);  

% specificity, 1 - prob of act given not this task category
px = p; nx = n;
px(2) = []; nx(2) = [];
spec = 1 - (sum(px .* nx) ./ sum(nx)); 
fprintf('Positive predictive value: (p task given act) for PAIN\n%3.2f\n', pp(2))

fprintf('Positive predictive value: (p task given act)\n')
disp(pp);

%% CHi2 test

counts = [OUT.nstudies(whdb) - OUT.nyes(whdb); OUT.nyes(whdb)]
[chi2,df,p,sig,warn,freq_table,expected_table,isnonparametric] = chi2test(counts,'table')

diary off

cd(basedir)
%%
end
