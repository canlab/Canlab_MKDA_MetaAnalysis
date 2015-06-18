function [prob,perc] = indic2classify1(indic,n,nms)
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
% indic2classify1(clnew.bmatx,clnew.n,nms);

% make sure all indic entries are 1 / 0
indic = double(indic > 0);

% remove rows (studies) with no task class
indic(find(all(indic(:,1:n)==0,2)),:) = [];

nreg = size(indic,2) - n;
reg = indic(:,n+1:end);
regnms = nms(n+1:end);

singlereg = reg;


% code for combinations of regions (pairs of 2)
sz = size(reg,2);
for i = 1:sz - 1
    for j = (i+1):sz
        
        clear tmp
        tmp = reg(:,i) & reg(:,j);
        %tt(i,j) = 1; imagesc(tt);drawnow
        if any(tmp)
            regnms{end+1} = [regnms{i} '+' regnms{j}];
            reg(:,end+1) = tmp;
        end
    end
end

indic = [indic(:,1:n) reg];
nms = [nms(:,1:n) regnms];

% perc: percent studies activating regions, rows are regions and combos of
% regions, columns are tasks
%perc = burt_profile([indic' * indic],nms,n);
burt2 = [indic' * indic];
perc = diag(1./diag(burt2)) * burt2; % this normalizes cov by row var %corrcoef(imatx);
perc = perc(1:n,n+1:end)';

scount = burt2(1:n,n+1:end);   % count studies, so that we can ensure that at least n studies activated, 
% to constrain predictions to combinations with reasonable counts

% prob: probability of being task x (row) given region (column) or combo of
% regions (column)
tmp = perc' ./ repmat(sum(perc'),size(perc,2),1);  
prob = tmp;

% plot this
%figure('Color','w'); set(gca,'FontSize',18)
wh = 1 : size(singlereg,2);
v = prob(:,wh)';
%plot(v,'s-','LineWidth',2); 
%legend(nms{1:n})
%set(gca,'XTick',1:length(wh))
%set(gca,'XTickLabel',nms(wh+n),'XLim',[.5 length(wh)+.5])


stackedbar2(v,[],nms(wh+n))
ylabel('Brain region'),title('Predictions of task given brain activity')
xlabel('Probability') 
set(gcf,'Position',[1363         621        1053         490])
saveas(gcf,'meta_meta_predict_single','fig')
saveas(gcf,'meta_meta_predict_single','tif')


% eliminate low-count combinations and save highly predictive ones

for i = 1:size(prob,2)
    tmp = find(prob(:,i) == max(prob(:,i)));
    if isempty(tmp), tmp(1) = 0;,maxpcount(i) = 0;
    else
 
        maxp(i) = tmp(1);
        maxpcount(i) = scount(maxp(i),i);
    end
end

% select combos
wh = size(singlereg,2)+1 : size(prob,2);
v = prob(:,wh)';
combnms = nms(wh+n);
combmp = maxpcount(wh);

omit = find(combmp < 3);
wh(omit) = [];
v(omit,:) = [];
combnms(omit) = [];
combmp(omit) = [];

stackedbar2(v,[],combnms)
ylabel('Brain region'),title('Predictions of task given brain activity')
xlabel('Probability') 
set(gcf,'Position',[1363         621        1053         490])
%figure('Color','w'); set(gca,'FontSize',18)

%plot(v,'s-','LineWidth',2); 
%legend(nms{1:n})
%set(gca,'XTick',1:length(wh))
%set(gca,'XTickLabel',combnms,'XLim',[.5 length(wh)+.5])
%ylabel('Probability') 
%xlabel('Probability'),title('Predictions of task given brain activity')
saveas(gcf,'meta_meta_predict_combos','fig')
saveas(gcf,'meta_meta_predict_combos','tif')



disp('Prob: columns show each region, rows are prob of being task x')
tmp;

% select for display and print table 
%
thr = .1;
tmp = sort(tmp .* perc');
wh = find(tmp(end,:) - tmp(end-1,:) > thr);
tmp = prob(:,wh) .* perc(wh,:)';
if isempty(tmp), tmp = zeros(n,1);,end
tnms = regnms(wh);
maxtmp = max(tmp);

for i = 1:n 
    fprintf(1,'Regions in which activation predicts task: %s\n\t', nms{i})
    fprintf(1,'\tPercent studies\t\tProb task | region\n\t\t')
    % header
    for k = 1:n
        fprintf(1,'%s\t', nms{k}) 
    end
            
    for k = 1:n          
        fprintf(1,'%s\t', nms{k})
    end
            
    fprintf(1,'\n')

    if ~any(tmp(i,:) == maxtmp), fprintf(1,'\t\tNo significant regions.\n');,
    else
        for j = 1:size(tmp,2),  % for each significant element
            if tmp(i,j) == max(tmp(:,j))
                % print row for this one
                fprintf(1,'\t%s\t', tnms{j})
                for k = 1:n
                   fprintf(1,'%3.2f\t', perc(wh(j),k))
                end
                for k = 1:n
                   fprintf(1,'%3.2f\t', prob(k,wh(j)))
                end
                fprintf(1,'\n')
            end
        end
    end
end


% TREE

if ~all(sum(indic(:,1:n),2) == 1)
    disp('Not all task indicators are unique, so categories are fuzzy.  Skipping tree.')
else
    
    for i = 1:size(indic,1),names{i} = nms{find(indic(i,1:n))};,end

    t = treefit(reg(:,1:nreg), names);
    treedisp(t,'names',regnms(1:nreg));

    [c,s,n,best] = treetest(t,'cross',reg(:,1:nreg),names);
    tmin = treeprune(t,'level',best);
    treedisp(tmin,'names',regnms(1:nreg));

end



return



