function [m,aprime,corr,far,misclass, mprop, statstable] = confusion_matrix(true,cla)
% [m,aprime,corr,far,missclass, mprop, statstable] = confusion_matrix(true,cla)
%
% calculates a confusion matrix (m, mprop) and a-prime scores (AUC
% effect size)
%
% from true class assignments (true) and
% estimated class assignments (cla)
%
% a-prime is Z(hr) - Z(far) for classification, like d-prime
% - corr/far of 0.0/1.0 are adjusted to 0.005/.995 to avoid unrealistic
% aprime values. Max a-prime is 5.15.
%
% corr is sensitivity, hit rate for each class, same as accuracy
% far is false alarm rate for each class. specificity = 1-far
%
% also returns vector of which observations were misclassified
%
% m, mprop:
% rows are actual (truth) values, columns are predicted (classified)
% values, in ascending order of numeric codes
% Accuracy is diag(mprop), same as sensitivity
%
% m = observation counts
% mprop = proportions of observations.
%
% corr: correct classification rate for true values
% True class = 1, 2, etc.
%
% far: false alarm rate for true values
% % of time a non-class i was classified as class i
%
% Other quantities in statstable output:
% ppv: Positive preditive value
% cases: Num cases in each condition, sum(m, 2)
% freq: Frequency of each condition, sum(m, 2) ./ sum(m(:))
%
% print output: e.g.,
% print_matrix(m, {'Pred. Innocent' 'Pred. Guilty'}, {'Actual Innocent' 'Actual Guilty'});
%
% Output table:
% mplus = [m corr'];
% total_class_rate = 1 - (sum(missclass) ./ length(missclass));
% mplus = [mplus; [far total_class_rate]];
% print_matrix(mplus, {'Pred. Innocent' 'Pred. Guilty' 'Correct Class Rate'}, {'Actual Innocent' 'Actual Guilty' 'False Alarm Rate'});

% Programmers' notes:
%
% tor wager, 7/08/04
% eliminate elements with true class of 0: mod 10/2/06
%
% tor : 4/2018.  corr/far of 0.0/1.0 are adjusted to 0.005/.995 to avoid unrealistic
% aprime values.

wh = (true ~= 0);
n = length(true);
misclass = false(n,1);

true = true(wh);
cla = cla(wh);

u = unique(true);

for i = 1:length(u)
    
    n(i) = sum(true == u(i));
    
    corr(i) = sum(true == u(i) & cla == u(i)) ./ sum(true == u(i));
    far(i) = sum(true ~= u(i) & cla == u(i)) ./ sum(true ~= u(i));
    
    % rows are actual (truth) values, columns are predicted (classified)
    % values
    for j = 1:length(u)
        m(i,j) = sum(true == u(i) & cla == u(j));
        
        mprop(i,j) = m(i,j) ./ n(i);
    end
    
end

% ppv: Positive preditive value
ppv = corr ./ (corr + far);

    

% a-prime: avoid perfect scores
corr_for_aprime = corr;
corr_for_aprime(corr_for_aprime == 0) = .005;
corr_for_aprime(corr_for_aprime == 1) = .995;

far_for_aprime = far;
far_for_aprime(far_for_aprime == 0) = .005;
far_for_aprime(far_for_aprime == 1) = .995;

aprime = norminv(corr_for_aprime) - norminv(far_for_aprime);

% avoid Inf values due to perfect solutions
% aprime(abs(aprime) > 10) = sign(aprime(abs(aprime) > 10)) * 4;

misclass(wh) = (true - cla) ~= 0;

statstable = table(sum(m, 2), sum(m, 2) ./ sum(m(:)), corr', (1-far)', ppv', aprime', 'VariableNames', {'Freq' 'Perc_cases' 'Sens' 'Spec' 'PPV' 'Aprime'});

end