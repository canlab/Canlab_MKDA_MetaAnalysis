function obj = get_error(obj, test_class)
% obj = get_error(obj, test_class)
%
% This function returns three kinds of error rate.
% b_err is the balanced misclassification rate.
% t_err is the total (overall) misclassification rate.
% p_err is the class probability err
%       
% Misclassification rates
% -------------------------
% we might like the balanced misclassification rate or the total
% misclassification rate.  
% balanced error rates are preferred if we want to get each class right and
% there are unequal numbers of instances in each class.
%
% class probability err
%
% Misclassification rates are a coarse basis for selecting models 
% because they are not smooth.  What we'd rather do is pick the model with
% the highest posterior probability assigned to the correct class.
% class probability err is defined as the -sum(log(p(class|act)) across
% voxels, i.e., -posterior_probs for the correct class. (note: nbc_predfun
% now scales this by dividing by voxels.)

% truth x estimate confusion matrix
%cm = confusionmat(test_class, class_pred) % matlab's version

if isempty(obj.class_pred)
    error('You must test the nbc object before running: nbc_obj = test(nbc_obj, test_data)');
end

if length(obj.class_pred) ~= length(test_class)
    error('obj.class_pred and test_class are different lengths.');
end
[cm, dprime, hr, far, misclass]  = confusion_matrix(test_class, obj.class_pred);  

% hr = diag(cm) ./ sum(cm, 2)

% balanced misclass rate
b_err = mean(1 - hr);

% total misclass rate
t_err = sum(misclass) ./ length(misclass);

% Class probability err
for i = 1:length(test_class)
p_err(i) = obj.posterior_probs(i, test_class(i)); % maximizing this is good...
end
p_err = -mean(p_err); % minimizing this is good... -1 is perfect.

obj.test_class = test_class;
obj.p_err = p_err;
obj.b_err = b_err;
obj.t_err = t_err;

fprintf('Updated test_class, p_err, t_err, and b_err in meta_nbc object.\n');

end % function
