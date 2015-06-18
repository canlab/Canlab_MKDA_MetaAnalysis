function [m,dprime,corr,far,misclass, mprop] = confusion_matrix(true,cla)
    % [m,dprime,corr,far,missclass, mprop] = confusion_matrix(true,cla)
    %
    % calculates a confusion matrix (m, mprop) and d-prime scores
    % from true class assignments (true) and
    % estimated class assignments (cla)
    %
    % d-prime is actually a-prime, the Z(hr) - Z(far) for classification
    % corr is hit rate for each class
    % far is false alarm rate for each class
    %
    % also returns vector of which observations were misclassified
    %
    % tor wager, 7/08/04
    % eliminate elements with true class of 0: mod 10/2/06
    %
    % m, mprop:
    % rows are actual (truth) values, columns are predicted (classified)
    % values, in ascending order of numeric codes
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
    % print output: e.g.,
    % print_matrix(m, {'Pred. Innocent' 'Pred. Guilty'}, {'Actual Innocent' 'Actual Guilty'});
%
% Output table:
% mplus = [m corr']; 
% total_class_rate = 1 - (sum(missclass) ./ length(missclass));
% mplus = [mplus; [far total_class_rate]];
% print_matrix(mplus, {'Pred. Innocent' 'Pred. Guilty' 'Correct Class Rate'}, {'Actual Innocent' 'Actual Guilty' 'False Alarm Rate'});

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

    dprime = norminv(corr) - norminv(far);

    % avoid Inf values due to perfect solutions
    dprime(abs(dprime) > 10) = sign(dprime(abs(dprime) > 10)) * 4;

    misclass(wh) = (true - cla) ~= 0;

    return