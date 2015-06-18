function [OUT] = indic2classify2(indic,n,nms)
%
% indic is matrix of indicators, studies/contrasts are rows, 1/0 codes yes/no for
% task and regional activation.  integers in indic above 1 are treated as
% 1.
%
% first n indicator vars are task codes, rest are regional activation
% nms is cell array of strings, first tasks, then regions
%
% You can use output of db_cluster_burt_table:
% indic2classify1(OUT.bmatx,OUT.n,OUT.nms)
%
% Outputs:
% perc: percent studies activating regions, rows are regions and combos of
% regions, columns are tasks;
% assumes you have one row for each study/contrast!!  if you remove rows
% corresponding to non-activating studies, the output here will not be the
% % of studies.
%
% prob: probability of being task x (row) given region (column) or combo of
% regions (column)
%
% Examples:
% indic2classify2(clnew.bmatx,clnew.n,nms);

% make sure all indic entries are 1 / 0
indic = double(indic > 0);

% remove rows (studies) with no task class
indic(find(all(indic(:,1:n)==0,2)),:) = [];

nreg = size(indic,2) - n;
reg = indic(:,n+1:end);
regnms = nms(n+1:end);

singlereg = reg;

% separate condition and data, make condition function (class vector)
% make sure that fuzzy categories are eliminated properly - remove them or
% something.
for i = 1:n
    indic(:,i) = indic(:,i) * i;
end
tmp = indic(:,1:n)';
om = find(sum(tmp) ~= max(tmp));    % omit rows with more than one class
clas = sum(tmp)';
clas(om) = NaN;                     % may cause knn to crash if NaNs?

indic = indic(:,n+1:end);

% eliminate contrasts with no activation in any region?
osize = size(indic,1);
whsav = find(any(indic ~= 0,2));
wh = find(all(indic == 0,2));
indic(wh,:) = []; clas(wh) = [];

%indic(indic == 0) = -1;
%indic = scale(scale(indic,1)',1)';  % double-center
%indic = scale(indic',1)';   % center rows only

OUT = tor_knn(indic,clas,3);

% save values in full length vector of all contrasts
classest = zeros(osize,1); misclass = classest;
classest(whsav) = OUT.class;
misclass(whsav) = OUT.misclass;
OUT.misclass = misclass;
OUT.class = classest;

try
figure;hold on
mean(indic(clas==1,:)); plot(ans,'bo-','LineWidth',2)
mean(indic(clas==2,:)); plot(ans,'go-','LineWidth',2)
mean(indic(clas==3,:)); plot(ans,'ro-','LineWidth',2)
mean(indic(clas==4,:)); plot(ans,'yo-','LineWidth',2)
catch
end

return

return



