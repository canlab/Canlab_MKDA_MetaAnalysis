function PRED = Meta_SVM_from_mkda(DB, fieldname, names_cell, varargin)
% Take output from Meta_Setup (DB) and
% runs a one-vs-one nonlinear SVM  classifier, on labels of your choice.
%
% [nbc, data] = meta_SVM_from_mkda(MC_Setup, DB, fieldname, names_cell)
%
% Note: You need the spider package on your matlab path.
%
% Inputs:
% DB: Meta-analysis (MKDA) database structure 
%     - required fields: DB.x, y, z, Contrast, Study, named field to classify (cell array of strings)
%     - see Meta_Setup and read_database.m
% fieldname : Field in database to use for classification
%
% Optional inputs:
% -----------------------------------------------------------------
% Enter variable name in quotes followed by new value for any of the options below:
% nfolds = 10;              % 'nfolds' followed by number of folds
% rbf_sigma = 15;           % 'rbf_sigma' followed by st dev in mm
% balanced_ridge = 8;       % arbitrary; 0 is no penalty for unbalanced classes
% holdout_meth = 'kfold';   % or 'study' for leave one study out
% test_method = 'distance'; % 'vote' or 'distance'; method for multi-class classifications
%                           % vote picks lowest class # in case of ties so
%                           is less preferred, but very similar to
%                           'distance' in practice.
%                           % 'distance' sums the dist from hyperplane
%                           across all pairs to end up with a dist measure
%                           for each class, and the max is selected.
%
% Examples:
% -----------------------------------------------------------------
% [nbc, data] = Meta_SVM_from_mkda(DB, MC_Setup, 'Valence', {'positive' 'negative'})
%
% fieldname = 'Stimuli';
% [nms, contrastcounts] = meta_explore_field(DB, fieldname);
% labels = nms(contrastcounts > 20)';
% [nbc, data] = Meta_SVM_from_mkda(DB, MC_Setup, fieldname, labels);
%
% fieldname = 'Emotion';
% names_cell = {'sad' 'happy' 'anger' 'fear' 'disgust'};
% PRED = Meta_SVM_from_mkda(DB, fieldname, names_cell, 'balanced_ridge', balridgeval, 'holdout_meth', 'study');
%
%PRED = Meta_SVM_from_mkda(DB, fieldname, names_cell, 'balanced_ridge', 8, 'holdout_meth', 'study');
%PRED = Meta_SVM_from_mkda(DB, fieldname, names_cell, 'balanced_ridge', balridgeval, 'nfolds', 3);
%
% Tor Wager, Feb 2013

% Parameters -- hard-coded
% ---------------------------------------------------------------------

nfolds = 10;
rbf_sigma = 15; % in mm
balanced_ridge = 8; % arbitrary; 0 is no penalty for unbalanced classes
holdout_meth = 'kfold'; % or 'study' for leave one study out
test_method = 'distance';  % 'vote' or 'distance'; method for multi-class classifications

% Optional inputs
% ---------------------------------------------------------------------
for opt_args = {'balanced_ridge' 'nfolds' 'rbf_sigma' 'holdout_meth' 'test_method'}
    
    wh = strcmp(varargin, opt_args{1});
    if any(wh)
        wh = find(wh); wh = wh(1);
        eval([varargin{wh} ' = varargin{wh+1};']);
    end
    
end

% set up empty svm object for training
% need spider package!
% svmobj = svm({kernel('rbf', rbf_sigma), 'C=1', 'optimizer="andre"'});
% balanced ridge=1 only shifts things a little bit...
svmobj = svm({kernel('rbf', rbf_sigma), 'C=1', 'optimizer="andre"', ['balanced_ridge=' sprintf('%d', balanced_ridge)]});

% GET RELEVANT DATA
% Get list of coords matching each category of input
% ---------------------------------------------------------------------

cfun = @(name) strcmp(DB.(fieldname), name);
wh = cellfun(cfun, names_cell, 'UniformOutput', 0);

% wh : which coordinates to include
wh = cat(2, wh{:});
wh = any(wh, 2);

outcome = DB.(fieldname);
outcome = outcome(wh);
xyz = [DB.x DB.y DB.z];
xyz = xyz(wh, :);

contrastnum = DB.Contrast;
contrastnum = contrastnum(wh);

PRED.PARAMS = struct('nfolds', nfolds, 'rbf_sigma', rbf_sigma, 'svmobj', svmobj);

PRED.DATA = struct('wh_coords', wh, 'outcome', {outcome}, 'xyz', xyz, 'contrastnum', contrastnum);
PRED.DATA.outcome = PRED.DATA.outcome{1};

nclasses = length(names_cell);

% ---------------------------------------------------------------------
% true outcomes, by contrast
% ---------------------------------------------------------------------

u = unique(contrastnum);

for i = 1:length(u)
    truevals = outcome(contrastnum == u(i));
    true_val(i, 1) = truevals(1);
    
    % check for errors in data entry
    if any(~strcmp(truevals, truevals{1}))
        warning(['CONTRAST ' num2str(i) ' VALUES DO NOT ALL HAVE SAME OUTCOME VALUE.']);
    end
    
end

% Point counts
cn = PRED.DATA.contrastnum;
for i = 1:length(u)
    npoints(i, 1) = sum(cn == u(i));
end

PRED.CONTRASTDATA = struct('u', u, 'u_descrip', 'unique contrast numbers');
PRED.CONTRASTDATA.outcome = true_val;
PRED.CONTRASTDATA.point_counts = npoints;

% ---------------------------------------------------------------------
% holdout sets:
% ---------------------------------------------------------------------
switch holdout_meth
    case 'kfold'
        % stratified partition, leave out whole contrasts
        [trIdx, teIdx] = define_holdout_sets(true_val, nfolds, outcome, contrastnum);
        
    case 'study'
        [indic, nms, condf] =  string2indicator(DB.Study(wh));
        nfolds = size(indic, 2);
        
        fprintf('Leave one Study out cross-val: %d folds\n', nfolds);
        for i = 1:nfolds
            trIdx{i} = ~indic(:, i);
            teIdx{i} = logical(indic(:, i));
        end
        
    otherwise
        error('invalid entry for holdout set.');
end

% ---------------------------------------------------------------------
% ONE VS. ONE CLASSIFICATION
% ---------------------------------------------------------------------

% will classify each test COORDINATE
% then we aggregate the coordinates using a voting (or distance) method,
% by taking the sum of distance from the hyperplane across coodinates
% within each contrastnum, to obtain a classification and confidence for each
% contrastnum.

% initialize vars
% ---------------------------------------------------------------------

[cv_class_est, cv_class_weights, cv_wh_teIdx] = deal(cell(nclasses));  % for each pair of classes

for emo1 = 1:nclasses-1
    
    for emo2 = emo1+1:nclasses
        
        [cv_class_est{emo1, emo2}, cv_class_weights{emo1, emo2}] = deal(NaN .* zeros(size(outcome)));
        
        
    end
end

% ---------------------------------------------------------------------
% ---------------------------------------------------------------------
% TRAIN/TEST: Loop through folds and pairs
% ---------------------------------------------------------------------
% ---------------------------------------------------------------------

for fold = 1:nfolds
    
    fprintf('Fold %d\n', fold)
    
    % X and Y: all observations for this fold
    
    [~, nms, Y] = string2indicator(outcome(trIdx{fold}),names_cell);
    X = xyz(trIdx{fold}, :);
    
    [~, nms, Yt] = string2indicator(outcome(teIdx{fold}),names_cell);
    Xt = xyz(teIdx{fold}, :);
    
    % train each pair, one vs. one
    % Select data for pair
    % x and y: observations relevant for a specific pair
    
    % Loop through pairs and train/test
    % ------------------------------------------------------------------
    for emo1 = 1:nclasses-1
        
        for emo2 = emo1+1:nclasses
            
            fprintf('%s vs. %s\t', names_cell{emo1}, names_cell{emo2});
            
            % Get pair data
            % ----------------------------------------------------------
            y = zeros(size(Y));
            y(Y == emo1) = 1;
            y(Y == emo2) = -1;
            
            if all(y == 0)
                % no train points in this comparison
                fprintf('No training points.\n')
                continue
            end
            
            % wh_omit: points not part of this specific comparison emo1, emo2
            wh_omit = y == 0;
            y(wh_omit) = [];
            
            x = X;
            x(wh_omit, :) = [];
            
            dat = data(x, y); %<- for input into spider package
            
            
            y = zeros(size(Yt));
            y(Yt == emo1) = 1;
            y(Yt == emo2) = -1;
            
            % Now: Test ALL points in test set for every classifier.
            % We do not know for the test set what the true class is...
            % so we cannot compare only on the decisions involving the true
            % class.
            %                         if all(y == 0)
            %                             no test points in this comparison
            %                             fprintf('No test points.\n')
            %                             continue
            %                         end
            %
            %                          wh_omit = y == 0;
            %                          y(wh_omit) = [];
            %
            x = Xt;
            %                          x(wh_omit, :) = [];
            
            dattest = data(x, y); %<- for input into spider package
            
            wh_teIdx = find(teIdx{fold});
            %                          wh_teIdx(wh_omit) = [];
            
            
            % Train
            % ----------------------------------------------------------
            
            [restrain, svmtrained] = train(svmobj, dat);
            
            % Test and get loss <- res.X is predicted class
            % ----------------------------------------------------------
            
            res = test(svmtrained,dattest);
            
            % Save results in array
            % ----------------------------------------------------------
            %L = loss(res2, 'class_loss')
            cv_class_est{emo1, emo2}(wh_teIdx) = res.X;
            
            
            % NOTE:  try to get weights for class probability estimates and better
            % voting of contrastnum emotion type based on points. RBF is a really good idea here.
            % but i think we can do it using use_signed_output = 0
            
            svmtrained.use_signed_output = 0;
            res = test(svmtrained,dattest);
            
            cv_class_weights{emo1, emo2}(wh_teIdx) = res.X; % <- useful for multi-way classification
            
            %cv_wh_teIdx{emo1, emo2} = wh_teIdx;
            
            % cv_wh_teIdx
            % Save which estimates are directly relevant for classification
            % i.e., class1-2 or class2-3 for class 2. class 4-5 is not
            % relevant. Using this to estimate accuracy supposes that you
            % know something about the test class, i.e., that it's either 2
            % or something else, so is circular in straight-up multiclass
            % classification.
            %             wh_omit = y == 0;
            %             wh_teIdx(wh_omit) = [];
            %             cv_wh_teIdx{emo1, emo2} = wh_teIdx;
            
        end
    end
    
end % folds

PRED.PAIRS = [];
PRED.PAIRS.cv_class_est = cv_class_est;
PRED.PAIRS.cv_class_weights = cv_class_weights;
PRED.PAIRS.cv_wh_teIdx = cv_wh_teIdx;

% ---------------------------------------------------------------------
%  Get accuracy and cv weights (distances) by contrast (in contrastnum)
% ---------------------------------------------------------------------

% First, get list of cross-val classifier weights by contrast
% then, figure out whether sign of weights for each contrast matches true
% category
% 'weights' here is really distance from hyperplane

disp(' ');
[cvweights, cvN, cvstd] = deal(cell(nclasses));

for emo1 = 1:nclasses-1
    
    for emo2 = emo1+1:nclasses
        
        % initialize
        cvweights{emo1, emo2} = NaN .* zeros(size(u)); % <- cvweights is per contrast
        
        testw = cv_class_weights{emo1, emo2}; % <- testw is per peak coord
        
        % to select only relevant comparisons:
        % this would replicate earlier biased multiclass method for testing
        % purposes
        % comment out after testing. BUT THIS IS NOT WORKING RIGHT...DEBUG
        % IF TESTING IS NEEDED
        %         testw2 = NaN .* zeros(size(cv_class_weights{emo1, emo2}));
        %         testw2(cv_wh_teIdx{emo1, emo2}) = testw(cv_wh_teIdx{emo1, emo2});
        %         testw = testw2;
        
        wh = ~isnan(testw);
        
        numcons = length(unique(contrastnum(wh)));
        fprintf('%s vs. %s\t:\t%3.0f coords, %3.0f contrasts\n', names_cell{emo1}, names_cell{emo2}, sum(wh), numcons);
        
        for c = 1:length(u)
            
            this_con = contrastnum == u(c);
            
            this_con = this_con & wh;  % this contrastnum AND this pair
            
            if ~any(this_con)
                % empty
                [cvw, cvs] = deal(NaN);
                cvnum = 0;
            else
                % average weights for this contrastnum across peaks
                con_weights = testw(this_con);
                
                cvw = sum(con_weights);
                cvnum = length(con_weights);
                cvs = std(con_weights);
            end
            
            cvweights{emo1, emo2}(c, 1) = cvw;      % average weight (distance from hyperplane)
            cvN{emo1, emo2}(c, 1) = cvnum;          % number of peaks
            cvstd{emo1, emo2}(c, 1) = cvs;          % weight sd
            
        end % contrastnum
    end
    
end

PRED.CONTRASTPAIRS = [];
PRED.CONTRASTPAIRS.cvweights = cvweights;
PRED.CONTRASTPAIRS.cvN = cvN;
PRED.CONTRASTPAIRS.cvstd = cvstd;

% ---------------------------------------------------------------------
% Pairwise accuracy overall
% ---------------------------------------------------------------------

[cvtrue, cvcorrect] = deal(cell(nclasses));
cvaccuracy = zeros(nclasses);

for emo1 = 1:nclasses-1
    
    for emo2 = emo1+1:nclasses
        
        [cvtrue{emo1, emo2}, cvcorrect{emo1, emo2}] = deal(NaN .* zeros(size(true_val)));
        
        cvtrue{emo1, emo2}(strcmp(true_val, names_cell{emo1})) = 1;
        cvtrue{emo1, emo2}(strcmp(true_val, names_cell{emo2})) = -1;
        
        cvcorrect{emo1, emo2} = double(cvtrue{emo1, emo2} == sign(cvweights{emo1, emo2}));
        %cvcorrect{emo1, emo2}(isnan(cvtrue{emo1, emo2})) = NaN;
        
        % cvaccuracy is pairwise acc. 
        % cvtrue has NaNs where truth is neither of two classes being
        % compared.
        cvaccuracy(emo1, emo2) = sum(cvcorrect{emo1, emo2}) ./ sum(~isnan(cvtrue{emo1, emo2}));
        
    end
    
end

PRED.CONTRASTPAIRS.cvtrue = cvtrue;
PRED.CONTRASTPAIRS.cvcorrect = cvcorrect;
PRED.CONTRASTPAIRS.cvaccuracy = cvaccuracy;
PRED.CONTRASTPAIRS.cvaccuracy_descrip = 'Pairwise CV accuracy';

% ---------------------------------------------------------------------
% ---------------------------------------------------------------------
% get class estimates
% PRED.class_est
% ---------------------------------------------------------------------
% ---------------------------------------------------------------------

% 'vote': vote method to get predicted class in k-way classification
% for confusion matrix
% Max wins voting method: j friedman. Another approach to polychotomous  classification. Tech report, Stanford U, 1996
% One vs one: knerr, dreyfus, 1990 "Single-layer learning revisited" Neurocomputing
% ---------------------------------------------------------------------

% ---------------------------------------------------------------------
% 'distance': Sum distance from hyperplane method to get predicted class
% in k-way classification for confusion matrix
% Introduced for one-against-all in Vapnik 1998
% Vapnik, V. (1998). Statistical Learning Theory. Wiley-Interscience, NY.
%
% For one vs. one, could use "tournament method" or DAGSVM, but not
% necessarily better
% Kijsirikul, B. & Ussivakul, N. (2002) Multiclass support vector machines using adaptive directed acyclic graph. Proceedings of International Joint Conference on Neural Networks (IJCNN 2002), 980-985.
%
% ---------------------------------------------------------------------


% 2 vote methods: mode of pairwise, or count votes for vs. against
% should produce identical results (tor tested 3/2013, yes)
[votesbyclass, class_dist] = deal(zeros(length(u), nclasses));

indx = 1;

for emo1 = 1:nclasses-1
    
    for emo2 = emo1+1:nclasses
        
        w = cvweights{emo1, emo2};  % dist from hyplane really
        
        yesno = sign(w);
        
        % emo1 (row) is "on" class, emo2 is "off" class
        votesbyclass(:, emo1) = votesbyclass(:, emo1) + yesno;
        votesbyclass(:, emo2) = votesbyclass(:, emo2) - yesno;
        
        % convert to class number
        yesno(yesno == 1) = emo1;
        yesno(yesno == -1) = emo2;
        
        class_est(:, indx) = yesno;  % contrast x class, vote method
        
        % emo1 (row) is "on" class, emo2 is "off" class
        class_dist(:, emo1) = class_dist(:, emo1) + w;
        class_dist(:, emo2) = class_dist(:, emo2) - w;
        
        indx = indx + 1;

    end
    
end

% Final class estimates and confusion matrix
% ---------------------------------------------------------------------
if size(class_est, 2) > 1
    class_est1 = mode(class_est')';
else
    class_est1 = class_est;
end

[~, class_est2] = max(votesbyclass, [], 2); % <- should be identical to est1 above
[~, class_est3] = max(class_dist, [], 2);   % <- distance method

switch test_method
    case 'vote'
        PRED.class_est = class_est1;
    case 'distance'
        PRED.class_est = class_est3;
    otherwise error('unknown test method.')
end

PRED.class_est_pairwise = class_est;
PRED.class_est_pairwise_descrip = 'Estimated class, contrasts x pairs';

% ---------------------------------------------------------------------
% Print output
% ---------------------------------------------------------------------

[PRED.true_indic, PRED.true_names, PRED.true_class] = string2indicator(true_val, names_cell);

fprintf('\nPairwise SVM accuracy:\n');
print_matrix(cvaccuracy, PRED.true_names, PRED.true_names, '%3.2f');
pairacc = cvaccuracy(:);
meanpairacc = mean(pairacc(pairacc ~= 0));

fprintf('Mean pairwise SVM accuracy: %3.2f\n', meanpairacc);

fprintf('\nMulticlass SVM accuracy:\n');

[m,dprime,corr,far,misclass, mprop] = confusion_matrix(PRED.true_class, PRED.class_est);

PRED.confusion_mat = m;
PRED.confusion_mat_prop = mprop;
PRED.acc = corr;

disp('Accuracy by class')
print_matrix(PRED.acc, PRED.true_names, PRED.true_names, '%3.2f');

disp('Confusion matrix - counts')
print_matrix(PRED.confusion_mat, PRED.true_names, PRED.true_names, '%3.0f');

fprintf('\n')
disp('Confusion matrix - proportions')
print_matrix(PRED.confusion_mat_prop, PRED.true_names, PRED.true_names, '%3.2f');

fprintf('Mean balanced accuracy: %3.2f\n', mean(PRED.acc))

% n points x accuracy

% [h, hc] = hist(npoints, 7);
% n = length(h);
% hc = [-Inf hc]; % should adjust to get proper bin edges
% for i = 1:n
% wh = npoints >= hc(i) & npoints < hc(i+1);
% acc(i, 1) = 1 - mean(misclass(wh));
% end


%% Point maps:
% We need to generate point-by-point estimates for every coordinate in the
% brain, and its contribution to the decision.
%
% this could yield two maps:  A MLC (most likely class; analogue of maximum
% a posteriori probability)
% ...and a weight map for each point.


end % function



function [trIdx, teIdx] = define_holdout_sets(true_val, nfolds, outcome, contrastnum)

% holdout sets
% ---------------------------------------------------------------------
% stratified partition
% cvpartition preserves the proportions in each class in train/test
% but does not force the number of training obs to be equal.
%
% Need to stratify by CONTRAST, so that an entire contrastnum at a time is
% left out.
% ---------------------------------------------------------------------

u = unique(contrastnum);

[trIdx, teIdx] = deal(cell(1, nfolds));

cvpart = cvpartition(true_val,'k',nfolds); % <- do this first in contrast space

% display
%print_matrix([], [{'Fold'} names_cell])

for i = 1:cvpart.NumTestSets
    
    trIdx{i} = cvpart.training(i);
    teIdx{i} = cvpart.test(i);
    
    % display
    %                 [indic,nms] = string2indicator(outcome(trIdx{i}),names_cell);
    %                 print_matrix(sum(indic), [], {['Train Fold ' num2str(i)]});
    %
    %                 [indic,nms] = string2indicator(outcome(teIdx{i}),names_cell);
    %                 print_matrix(sum(indic), [], {['Test Fold ' num2str(i)]});
    
    
end

% Now re-assign trIdx and teIdx in coordinate space, leaving out whole
% contrasts
for i = 1:cvpart.NumTestSets
    
    [coord_trIdx{i}, coord_teIdx{i}] = deal(false(size(outcome)));
    
    whtrain = find(trIdx{i});
    whtest = find(teIdx{i});
    
    for c = 1:length(whtrain)
        
        wh = contrastnum == u(whtrain(c));
        
        coord_trIdx{i}(wh) = true;
        
    end
    
    for c = 1:length(whtest)
        wh = contrastnum == u(whtest(c));
        
        coord_teIdx{i}(wh) = true;
        
    end
    
    
end

trIdx = coord_trIdx;
teIdx = coord_teIdx;

fprintf('Cross-validated prediction with SVM, %3.0f folds\n', nfolds)

end


