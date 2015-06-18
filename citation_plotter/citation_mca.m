function OUT = citation_mca(file)
% OUT = citation_mca(file)
% Citations are output from Endnote in APA 5th format
% Save as MS-DOS text
% example: 'citation_list1.txt'
%
% Multiple people with the same name, e.g., Zhang, are not counted
% correctly!

% read text file
s = textread(file,'%s','delimiter',',\n');

% find & and use as marker of title
% remove & from last author name

clear refmarker titlemarker, alltitle = [];

% refmarker marks start of each references, indexes cell BEFORE reference uses blank cells (carr. return)
% titlemarker 

for i = 1:length(s), 
    if isempty(s{i}), 
        refmarker(i) = 1;, 
        %titlemarker(i-3) = 1;, 
    elseif s{i}(1) == '&', 
        s{i} = deblank(s{i}(2:end));
        
    else
        s{i} = deblank(s{i});
    end
end

refmarker = find(refmarker);
%refmarker = [1 refmarker];  % was [0 refmarker]


% eliminate extra line breaks from refmarker list
tmp = find(diff(refmarker) == 1);
refmarker(tmp+1) = [];


% Authors and title of each study

aulist = [];
% get alphabetical author list
for i = 1:length(refmarker)-1 
    
    % get title and titlemarker
    [titles{i},tm] = titlemarker(s(refmarker(i)+1:refmarker(i+1)));
    
    aus{i} = s(refmarker(i)+1:refmarker(i)+tm-1);  aus{i} = aus{i}(1:2:end);  % every other, avoid initials
    
    % double-check to remove existing initial spaces
    for j = 1:length(aus{i}), 
        aus{i}{j} = deblank(aus{i}{j});,
        if aus{i}{j}(1) == ' ', aus{i}{j} = aus{i}{j}(2:end);,end   
    end
    
    aulist = str2mat(aulist,aus{i}{:});
end
aulist2 = (unique(aulist,'rows'));
aulist2 = aulist2(2:end,:);
titles2 = str2mat(titles{:});


% make Burt table
burt = zeros(length(aulist2));

fprintf(1,'Making BURT...\n')
% First do diagonals
for i = 1:length(aus)
    for j = 1:length(aus{i})     
        % find index of this author
        tmp = strmatch(aus{i}{j},aulist2,'exact');
        if length(tmp) > 1, keyboard, end
        burt(tmp,tmp) = burt(tmp,tmp) + 1;
        
    end    
end

for i = 1:length(aus), ausize(i) = length(aus{i});,end
fprintf(1,'%3.0f studies, average # authors = %3.2f, avg. papers per author = %3.2f\n',length(refmarker),mean(ausize),mean(sum(burt)))
tmp = find(sum(burt) == max(sum(burt)));
for i = 1:length(tmp),
    fprintf(1,'Most prolific author: %s with %3.0f papers.\n',aulist2(tmp(i),:),max(sum(burt)))
end

few = find(sum(burt) < 2);
burtr = burt; burtr(:,few) = []; burtr(few,:) = [];
aulistr = aulist2; aulistr(few,:) = [];

%%% OLD CODE - NOT NEEDED %%%%%
%%% OLD CODE - NOT NEEDED %%%%%
%%% OLD CODE - NOT NEEDED %%%%%

go = 0;
if go
    
fprintf(1,'Making BURT reduced (with off-diagonals)...\n')
% Then off-diagonals
for i = 1:length(aus)
    for j = 1:length(aus{i})  
        % find index of this author
        tmp = strmatch(aus{i}{j},aulistr,'exact');
        
        for k = j+1:length(aus{i})
            tmp2 = strmatch(aus{i}{k},aulistr,'exact');
            if ~isempty(tmp) & ~isempty(tmp2) & tmp ~= tmp2
                burtr(tmp,tmp2) = burtr(tmp,tmp2) + 1;
                burtr(tmp2,tmp) = burtr(tmp2,tmp) + 1;
            end
        end
    end    
end

end  % go

%%% END - OLD CODE - NOT NEEDED %%%%%


indic = zeros(length(aus),size(aulistr,1));
% Make study table
for i = 1:length(aus)
    for j = 1:length(aus{i})  
        % find index of this author
        tmp = strmatch(aus{i}{j},aulistr,'exact');
        if ~isempty(tmp)
            indic(i,tmp) = 1;
        end
    end
end

[indic2,scount,tnames,words] = title_keywords(titles2);
isize = size(indic,2);
indic = [indic indic2];

% double-center
burtr = indic' * indic;
burtr2 = scale(burtr')'; burtr2 = scale(burtr2);


for i = 1:size(aulistr,1), nms{i} = deblank(aulistr(i,:));,end
for i = 1:size(tnames,1), nms{end+1} = deblank(tnames(i,:));,end
%tor_mca(indic,'rowcol',isize,nms)


% run MCA


[pc,score,latent] = princomp(burtr2);

figure('Color','w'), set(gca,'FontSize',18),hold on
plot(score(:,1),score(:,2),'ws','MarkerFaceColor','w','LineWidth',0.5)
v = std(score(:,1:2));
for i = 1:size(score,1)
    if i <= isize, fcolor = 'k';, else, fcolor = 'r';,end
    text(score(i,1),score(i,2)+(rand(1)-.5)*v(2)/5,nms{i},'FontSize',round(14+burtr(i,i)./2),'Color',fcolor)
end



% output save
OUT.aulist2 = aulist2; OUT.titles2 = titles2;
%OUT.burt = burt;
OUT.aus = aus;
OUT.burtr = burtr2;
OUT.aulistr = aulistr;
OUT.pc = pc;
OUT.score = score;
OUT.indic = indic;

% sorted list
fprintf(1,'Top 20 most frequent\n')
tmp = aulistr; tmp2 = diag(burtr(1:isize,1:isize));
[dummy,i] = sort(1./tmp2);
tmp2 = tmp2(i);
tmp = tmp(i,:);
OUT.ausort = tmp;

for i = 1:min(20,length(tmp2))
    fprintf(1,'%s\t:\t%3.0f papers\n',tmp(i,:),tmp2(i));
end

return

