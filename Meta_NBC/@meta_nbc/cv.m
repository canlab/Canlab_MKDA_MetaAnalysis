function obj = cv(obj, meta_data_obj)
% Cross-validated error for meta_nbc object
% on dataset stored in meta_data_obj
%
% obj = cv(obj, meta_data_obj)
%
% ----------------------------------------------------------------------
% Cross-validation on meta dataset
% This can be used to optimize parameters for testing on other holdout maps
% We'd like to define this so that we can plug in different algorithms
% as prediction functions.  Therefore, this should simply run predfun for
% each fold. Different predfuns can be specified for different algorithms.
% ----------------------------------------------------------------------
% 
% Here we will optimize the cross-validated error, so the 'test' class
    % is the holdout set.  We can use this to optimize parameters.
    % Later, we'll make predictions for the real test set, which is new
    % images. See the test method for that.
% 
% obj.class_pred is the predicted class. This is returned cross-validated,
% and it is used in the report method.  See help meta_dataset for examples.
%
    
% k-fold validation or predict new maps. for simplicity, just treat all
% the new maps as one fold.

t1 = tic;
fprintf('Cross-validating with %3.0f folds\n', obj.params.k);

% Matlab's cross-val object: Stratified samples (when done this way)
% see below for info on selecting training/test obs.
cvobj = cvpartition(meta_data_obj.classes,'Kfold', obj.params.k);

n_con = size(meta_data_obj.dat, 2);
n_class = max(size(meta_data_obj.terms));

probs = zeros(n_con, n_class);
class_pred = zeros(n_con, 1);

[p_errobj.b_err t_err] = deal(zeros(obj.params.k, 1));

for j=1:obj.params.k
    
    fprintf('Training fold %d...', j);
    
    % Training and test indices
    obj.trainIdx = cvobj.training(j); % Used in train.m method
    testIdx = cvobj.test(j);     % Used to select data here

    % train classifier
    fold_obj = train(obj, meta_data_obj); % uses trainIdx
    
    test_data = meta_data_obj.dat(:, testIdx);
    
    fold_obj = test(fold_obj, test_data);  
    
    fold_obj = get_error(fold_obj, meta_data_obj.classes(testIdx));
    
    p_err(j,1) = fold_obj.p_err;
    t_err(j,1) = fold_obj.t_err;
    b_err(j,1) = fold_obj.b_err;
    
    class_pred(testIdx) = fold_obj.class_pred;
        
    fprintf(' %3.0f sec\n', toc(t1)); t1 = tic;
    
end  % fold loop

% Average across folds
b_err = cvobj.TestSize * b_err ./ sum(cvobj.TestSize);
t_err = cvobj.TestSize * t_err ./ sum(cvobj.TestSize);
p_err = cvobj.TestSize * p_err ./ sum(cvobj.TestSize);

obj.p_err = p_err;
obj.t_err = t_err;
obj.b_err = b_err;

obj.class_pred = class_pred;
obj.test_class = meta_data_obj.classes;

fprintf('Updated p_err, t_err, b_err, class_pred with cross-validated error.\n');

end
