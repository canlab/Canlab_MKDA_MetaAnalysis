function [chi2,actual,expected,colnames] = simplechi2(studycount,totalstudycount)
% function [chi2,colnames] = computechi2(studycount,totalstudycount)
%
% 07.22.01 by Tor Wager
%
% computes an a vs b chi2, where #nos is total - #yesses
% can generalize to other things, L/R for example, but must give L or R and total
%
% studycount is a column vector of counts for first condition
% totalstudycount is a col. vector of total (a + b) counts across the two conditions

warning off

	if size(totalstudycount,1) > 1
		overalltotalstudycount = totalstudycount;
		if ~(size(totalstudycount,1) == size(studycount,1) )
			error('Error in computechi2: totalstudycount rows must equal 1 or number of regions (rows in studycount)')
		end
	end
	
% --- loop through areas ----
for i = 1:size(studycount,1)
	
	% if different totals for each region, then get the totals for that region.
	if exist('overalltotalstudycount') == 1
		totalstudycount = overalltotalstudycount(i,:);
	end
		
% overall prob of 'yes' across conditions
% this will work only for columns that cannot overlap - i.e., works for vis/aud/recall, but not for l/r/c
% because some studies may find both of those, and sum(studycount(i,:) would be > sum(totalstudycount)

	pyes(i,1) = sum(studycount(i,:)) ./ sum(totalstudycount);
	
	% actual and expected: 
	% cols are 'y' 'n', rows are the comparison groups (vis, aud, etc).
	% one cell for each region
	actual{i}(:,1) = studycount(i,:)'; % yesses
	actual{i}(:,2) = totalstudycount' - studycount(i,:)'; % nos
	expected{i}(:,1) = totalstudycount' ./ 2;
	expected{i}(:,2) = totalstudycount' ./ 2;
	
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

warning on

return