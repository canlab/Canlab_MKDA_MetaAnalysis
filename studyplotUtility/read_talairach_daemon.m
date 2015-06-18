function [BA,gray,prange,num,tx,ty,tz,subc] = read_talairach_daemon(infile,varargin)
% function [BA,gray,prange,num,tx,ty,tz,subc] = read_talairach_daemon(infile,[opt] print text table)
%
% reads output from Talairach Daemon Database (v1.0)
% and outputs BA (Brodmann area), range from nearest label, and gray matter status
% by Tor Wager

[num,tx,ty,tz] = textread(infile,'%d%d%d%d%*[^\n]','delimiter',',');
res = textread(infile,'%s','delimiter','\n');

for i = 1:size(res,1)
    
    % get range of this coordinate from the nearest Brodmann gray matter label
    
    r = findstr('Range=',res{i});
    if ~isempty(r)
        prange(i,1) = str2num(res{i}(r+6:end));
    else
        prange(i,1) = NaN;
    end
    
    % get gray matter status (within search range, depending on options specified in TDD
    
    r = findstr('Gray',res{i});
    if ~isempty(r)
        gray(i,1) = 1;
    else
        gray(i,1) = 0;
    end
    
    % get Brodmann area, if any
    r = findstr('Brodmann area',res{i});
    if ~isempty(r)
        BA(i,1) = str2num(res{i}(r+13:r+15));
    else
        BA(i,1) = NaN;
    end
    
    % get subcortical label
    r = findstr('Amygdala',res{i});
    if isempty(r), r = findstr('Caudate',res{i});,end
    if isempty(r), r = findstr('Putamen',res{i});,end
    if isempty(r), r = findstr('Globus',res{i});,end
    if isempty(r), r = findstr('Substantia Nigra',res{i});,end
    if isempty(r), r = findstr('Midbrain',res{i});,end
    if isempty(r), r = findstr('Thalamus',res{i});,end
    if isempty(r), r = findstr('Hippocampus',res{i});,end
    if isempty(r), r = findstr('Cerebellum',res{i});,end
    if ~isempty(r)
        subc(i,1:3) = (res{i}(r(1):r(1)+2));
    else
        subc(i,1:3) = 'NaN';
    end
end



% print table, if specified
if nargin > 1
    
    fprintf(1,'Results for %s\n', infile)
    fprintf(1,'Index\tx\ty\tz\tGray\tBA\tRange\n')
    
    for i = 1:size(res,1)
        if gray(i) == 1, mystr = 'Gray matter';, else, mystr = 'White/ventricle?';,end
        fprintf(1,'%3.0f\t%3.0f\t%3.0f\t%3.0f\t%s\t%3.0f\t%3.0f\n', ...
            num(i),tx(i),ty(i),tz(i),mystr,BA(i),prange(i));
    end
    
end