function [xyzall,codesall,orderall] = average_nearby_xyz(XYZ,thresh,codes)
% [xyzall,codesall,orderall] = average_nearby_xyz(XYZ,thresh,codes)
%
% for plotting peak coordinates on the brain, for example, eliminating
% redundant nearby coordinates from the same study
%
% Fixed bug for SINGLE coord with unique code, Jan 2011: Tor Wager

% if codes are not a cell array of strings, make so...
if ~iscell(codes) 
    codes = mat2cell(codes, ones(size(codes, 1), 1), 1);
end

if ~isempty(codes{1}) && ~ischar(codes{1})
    for i = 1:length(codes)
        codes{i} = num2str(codes{i});
    end
end

u = unique(codes);


codesall = [];
xyzall = [];
orderall = [];

for i = 1:length(u)
    
    studywh = find(strcmp(codes,u{i}));
    wh = XYZ(studywh,:);
    
    dis = squareform(pdist(wh));

    if size(wh, 1) == 1  % fix for one unique point...
        dis = 0;
    end
    
    dis = dis < thresh;
    
    %howmany = sum(dis,2);
    %[dummy,order] = sort(1./howmany);
    %howmany = howmany(order,:);
    %dis = dis(order,:);
    
    [u2,order2] = unique(dis,'rows');
    
    xyzout = [];
    for j = 1:size(u2,1)
        
        tmpxyz = wh(find(u2(j,:) == 1),:);
        xyzout(j,:) = mean(tmpxyz,1);
        
    end
    
    xyzall = [xyzall; xyzout];
    codesall = [codesall; repmat(u(i),size(xyzout,1),1)];
    orderall = [orderall; studywh(order2)];
end



return