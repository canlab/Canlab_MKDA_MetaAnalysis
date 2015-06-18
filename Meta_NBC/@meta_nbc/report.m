function obj = report(obj)
% Print accuracy report for a meta_nbc object
% obj = report(obj)
%
% You should have already trained and tested the object using train and
% test with a meta_dataset data object and possibly test data.
% 
% The report method uses the last error values stored from either the test or
% cv methods.

if isempty(obj.b_err)
    error('You must test the object and attach err values (run obj = test(obj)) before reporting.');
end

[cm, aprime, hr, far]  = confusion_matrix(obj.test_class, obj.class_pred);
n_class = size(cm, 1);

fprintf('Post prob of correct class: %3.3f%%\n', 100*(-obj.p_err));

fprintf('Overall accuracy: %3.3f%%\n', 100*(1-obj.t_err));
fprintf('Balanced accuracy: %3.3f%%\n', 100*(1-obj.b_err));
fprintf('A-prime by class: %s\n', num2str(aprime));

fprintf('\nConfusion matrix:\nCounts:\n'); 

print_matrix(cm, obj.terms, obj.terms);%fprintf('\n');

fprintf('Proportions:\n');
print_matrix(cm ./ repmat(sum(cm, 2), 1, n_class), obj.terms, obj.terms);%fprintf('\n');

% significance (binomial on each class)
for i = 1:n_class
    iscorrect = obj.test_class == i & obj.class_pred == i;
    iscorrect = iscorrect(obj.test_class == i); % conditional on class
    
    RES{i} = binotest(iscorrect, 1 ./ n_class);
end

fprintf('\nSensitivity, specificity, binomial test on accuracy by class:\n')
for i = 1:n_class
    
    if length(obj.terms) < i || isempty(obj.terms{i})
        obj.terms{i} = num2str(i);
    end
    
    fprintf('Class %3.0f, %s:\t%3.3f\t%3.3f\tAcc = %3.3f +- %3.3f, P = %3.5f\n', i, obj.terms{i}, hr(i), 1-far(i), RES{i}.prop, RES{i}.SE, RES{i}.p_val);
end

end

