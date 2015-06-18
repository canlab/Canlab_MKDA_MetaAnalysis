function [indic,scount,su,words] = title_keywords(titles2)
% [indic,scount,su,words] = title_keywords(titles2)
%
% given a string matrix, stores individual words in each row, and finds the most frequent words.
% Then stores a matrix of which rows have which of the most frequent words.

t = titles2;

% get all words in titles

for i = 1:size(t,1)
    ti = titles2(i,:);
    
    try
        w1 = find(ti == ')'); w1 = w1(1); w2 = find(ti == '.'); w2 = w2(w2 > w1); w3 = w2(end-1)-1; w2 = w2(1)+1;
        % w2 is first period following ), w3 is 2nd to last period

        ti2 = ti(w2:w3);
        words{i} = getwords(ti2);
    catch
        words{i} = {'Cannot_parse_title.'};
        disp(['title keywords.m, problem at row ' num2str(i) ': probably not a title with two periods, one after year one at end.'])
    end
end

for i = 1:length(words), if isempty(words{i}), wh(i) = 1;, else, wh(i) = 0;, end, end
wh = find(wh);
words1 = words;
words(wh) = [];

% concatenate
allw = str2mat(words{1});
for i = 2:length(words)
    allw = str2mat(allw,str2mat(words{i}));
end

[u,a,b] = unique(allw,'rows');  % elements of b are counts for each word
for i = 1:length(u),cnt(i) = sum(b==i);,end % frequency of each word
[scount,wh] = sort(1./cnt);
su = u(wh,:);                   % sorted words in terms of frequency
scount = cnt(wh);          

% check accuracy of above code
%c = 0;for i=1:length(allw),if strcmp(deblank(allw(i,:)),'of'),c = c+1;,end,end


% save top n words - you pick
n = 25;

su = su(1:n,:); scount = scount(1:n);
for i = 1:n
    fprintf(1,'%3.0f   %3.0f   %s\n',i,scount(i),su(i,:))
end

wh = input('Enter vector of numbers of words to save: ');
su = su(wh,:); scount = scount(wh);

% indicator matrix
indic = make_indicator(words1,su);
        

return


function w = getwords(ti)

w4 = find(ti == ' ');
for j = 1:length(w4)-1, 
    w{j} = lower(deblank(ti(w4(j)+1:w4(j+1))));, 
    
    % exclude uninteresting characters
    w{j}([findstr(w{j},'(') findstr(w{j},')') findstr(w{j},'.') findstr(w{j},',') findstr(w{j},'''') findstr(w{j},':')]) = [];
    
    % exclude uninteresting words
    if strcmp(w{j},'a') |  strcmp(w{j},'to') | strcmp(w{j},'and') | strcmp(w{j},'then') | strcmp(w{j},'for') | ...
            strcmp(w{j},'by') | strcmp(w{j},'then') |  strcmp(w{j},'in') |  strcmp(w{j},'the') |  strcmp(w{j},'of') | ...
             strcmp(w{j},'activation') | strcmp(w{j},'responses') | strcmp(w{j},'central') | ...
             strcmp(w{j},'with') | strcmp(w{j},'an') | strcmp(w{j},'study') | strcmp(w{j},'during') | strcmp(w{j},'on')               

        w{j} = [];
    end
        
end
w = w';
for i = 1:length(w), if isempty(w{i}), wh(i) = 1;, else, wh(i) = 0;, end, end
wh = find(wh);
w(wh) = [];

return
