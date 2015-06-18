function [chi2,actual,expected,colnames] = computeynchi2(studycount,totalstudycount)
% function [chi2,actual,expected,colnames] = computechi2(studycount,totalstudycount)
%
% 07.22.01 by Tor Wager
%
% computes a yes/no chi2, where #nos is total - #yesses
%
% each row of studycount is the number of 'hits' for conditions 1, 2, 3, etc.
% so: [hc1 hc2 ch3, etc.]
%
% each row is another condition - condition 2
% a different region of the brain, or different comparison.
%
% this returns the chi2 for the effect of the first condition variable
% after adjusting for the overall proportion of hits across condition 2
% so it's a "hits to misses effect" at each level of condition 2
% 
%
% studycount is a 2D matrix: areas are rows, conditions are columns
% totalstudycount is a row vector: total counts for each row (i.e., in each condition)
%
% output: chi2:
% 1 chi2
% 2 df
% 3 p
% 4 cramer's V
% 5 sig flag
% 6 low expected value flag -- warning
% 7 prob of a yes overall
%
% examples
% do two tests, each across two conditions
% computeynchi2([12 20; 12 10],[20 20])
% do two tests, each across three conditions
% computeynchi2([12 20 25; 12 10 25],[20 20 40])

% --- loop through areas ----
for i = 1:size(studycount,1)
	
% overall prob of 'yes' across conditions
% this will work only for columns that cannot overlap - i.e., works for vis/aud/recall, but not for l/r/c
% because some studies may find both of those, and sum(studycount(i,:) would be > sum(totalstudycount)

	pyes(i,1) = sum(studycount(i,:)) ./ sum(totalstudycount);   % overall p of hit across all conditions
    % this is a normalization factor
	
	% actual and expected: 
	% cols are 'y' 'n', rows are the comparison groups (vis, aud, etc).
	% one cell for each region
	actual{i}(:,1) = studycount(i,:)'; % yesses
	actual{i}(:,2) = totalstudycount' - studycount(i,:)'; % nos
	expected{i}(:,1) = pyes(i,1) .* totalstudycount'; % p(yes) * row total
	expected{i}(:,2) = (1-pyes(i,1)) .* totalstudycount'; % p(no) * row total
	
	colnames = {'chi2' 'df' 'p' 'cramers-v' 'sig' 'warning' 'pyes'};
	
	%chi2 variable: chi2 value, df, p-value, effect size, significant, warning, pyes
	chi2(i,1) = sum(sum( ((actual{i} - expected{i}) .^ 2) ./ expected{i} ));
	chi2(i,2) = (length(totalstudycount) - 1) .* 1;		% #rows - 1 * #cols - 1
	chi2(i,3) = 1 - chsqcudf(chi2(i,1),chi2(i,2));		% p value
	
	if chi2(i,3) < .05, chi2(i,5) = 1; else  chi2(i,5) = 0;, end	% 1 if significant at a = .05
	if any(any(expected{i} <= 5)), chi2(i,6) = 1; else  chi2(i,6) = 0;, end	% 1 if warning for expected vals <= 5
	
	% effect size - cramer's V
	chi2(i,4) = sqrt( chi2(i,1) ./ (sum(totalstudycount) * min(length(totalstudycount),1)) );
	
end	% loop through areas

chi2(:,7) = pyes;


return