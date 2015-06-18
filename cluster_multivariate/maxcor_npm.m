function [mc,stats] = maxcor_npm(X,perms,varargin);
% [mc,stats] = maxcor_npm(X,perms,[c],[MCD robust outlier removal])
%
% X     data matrix, columns are variables, rows observations
% c     contrast matrix for rotation of X prior to correlation
%       contrasts should be rows



% ----------------------------------------------------------
% outliers
% ----------------------------------------------------------

if length(varargin) > 1
        [res]=fastmcd_noplot(X);
        
        % remove n most extreme outliers and recompute correlation
        wh = res.flag==0; nout = sum(res.flag==0);
        
        %res = fastmcd(X);
        %stats.rdthresh = input('Enter robust distance threshold for outliers:');
        %wh = res.robdist > stats.rdthresh;
        %nout = sum(wh);
        
        X(wh,:) = [];
        stats.nout = nout;
        fprintf(1,'Removing %3.0f outliers\n',stats.nout)
end

% ----------------------------------------------------------
% max cor
% ----------------------------------------------------------
if length(varargin) > 0, c = varargin{1};,else, c=eye(size(X,2));,end
mc = maxcor(X,c);
stats.n = size(X,1);
stats.cor = corrcoef(X*c');

for i = 1:perms,   % permute columns and test H0: no relation
    
    xtst=X;
    for j = 1:size(X,2)
        xtst(:,j) = getRandom(X(:,j));
    end
    [cctst(i,:),ccall(i,:)] = maxcor(xtst,c);

end

stats.mc = mc;

if perms
    stats.cctst = cctst;
    stats.perms = perms;
    stats.thresh = prctile(cctst(:,1),95);
    stats.thresh(2) = prctile(cctst(:,2),5);
    stats.bias = squareform(mean(ccall));
    stats.allupper = prctile(ccall,95);
    stats.alllower = prctile(ccall,5);
    
    
    % table and sig
    l = squareform(stats.alllower);
    u = squareform(stats.allupper);
    stats.sig = ~(stats.cor==1) & (stats.cor > u | stats.cor < l);
    stats.adjcor = stats.cor - stats.bias;
    
    str = sprintf(['Adjusted correlations\n']);
    for i = 1:size(stats.cor,1),
        for j=1:size(stats.cor,2)    %i,
            if stats.sig(i,j),t='*';,
            else,t='';,
            end,
            str=[str sprintf('%3.2f%s\t',stats.adjcor(i,j),t)];,
        end,
        str=[str sprintf('\n')];,
    end

    disp(str)           
end

return






function [mc,c] = maxcor(X,c)

X = X * c';
c = tril(corrcoef(X)); c(c==0 | c==1) = [];
mc = [max(c) min(c)];

return
