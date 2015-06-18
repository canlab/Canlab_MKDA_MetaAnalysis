function MCA = tor_mca(imatx,smeth,n,nms,decompmethod,niter)
    % MCA = tor_mca(imatx,smeth,n,nms,decompmethod,niter)
    %
    % multiple correspondence analysis
    % enter an indicator matrix and a scaling method
    %
    % *** NOW ONLY WORKS COMPLETELY FOR NMDS / CORR ***
    % Choices:
    % 1) Whether to center rows, columns, or both (default is both)
    % 2) Prin. comp., Metric or nonmetric multidimensional scaling (def. is nonmetric)
    %
    % imatx: indicator matrix of 1's and 0's
    % smeth: 'col' 'row' 'rowcol' 'corr' 'correspondence' 'phi' 'tau'
    %        phi is recommended for 1/0 data, tau good too for
    %        non-continuous; phi and tau appear identical for 1/0 data
    % n:     number of variables in first set
    %        now hard-coded for two sets of variables in diff colors
    % nms:   cell array of names of varables
    % decompmethod: decomposition method: pca or nmds (default) for
    %        nonmetric MDS
    % niter: iterations for clustering permutation test 
    %
    % I think the standard scaling is column; see Michalidis and Statsoft text
    % Not sure whether to analyze Burt table or correlation matrix;
    % sources seem to suggest Burt directly.
    %
    % MCA.v = Normalized values for cross correlations between
    % Set1 and Set2, rows are regions, cols are task conds
    % normed so that sums of squares for each diagonal is 1.
% 
% to make nmdsfig connectivity figure after running this:
% f1 =
% nmdsfig(MCA.coords,'classes',MCA.ClusterSolution.classes,'names',names,'sig',MCA.p < .05,'legend',{'Pos' 'Neg'},'sizescale',[4 16]);

    
    %imatx2 = imatx;
    %burt2 = imatx'*imatx;   % for last plot

    %imatx = imatx - mean(imatx(:));

    if nargin < 5, decompmethod = 'nmds'; end

    if nargin < 6, niter = 2000; end
    
    % -----------------------------------------------------
    % Input matrix scaling
    % -----------------------------------------------------

    imatx = full(double(imatx));

    [nobs,nvars] = size(imatx);
    if isempty(n), n = nvars; end
    
    if isempty(smeth), smeth = 'col'; end
    %if findstr(decompmethod,'mds'), smeth = 'corr'; end

    if findstr(smeth,'row')
        imatx = imatx - repmat(mean(imatx,2),1,nvars);  % rows
        disp('Centering rows.');
    end

    if findstr(smeth,'col')
        imatx = imatx - repmat(mean(imatx),nobs,1);    % cols
        disp('Centering columns.');
        
    elseif findstr(smeth,'none')
    elseif findstr(smeth,'corr')
        burt = corrcoef(imatx);
        D = (1 - burt) ./ 2;
     
    elseif findstr(smeth,'tau')   
        
        [corr,t,p] = correlation('tau',imatx);  %[corr,t,p] = c
        D = .5 .* (1 - corr);
        
    elseif findstr(smeth,'phi')   
        
        [corr,t,p] = correlation('phi',imatx);  %[corr,t,p] = c
        D = .5 .* (1 - corr);
        
    elseif findstr(smeth,'correspondence')
        
        % counts number both 1 + both 0 divided by n
        g = ones(nobs,nvars);
        sim = ((X'*X) + (g-X)' * (g-X)) ./ nobs;
        D = 1 - sim; % dissim, from 1 to zero

    else
        error('Unknown method in tor_mca.m')
    end

    % check for already scaled and adjust smeth
    %tol = 10E-6;
    %if all(mean(imatx) < tol), smeth = 'col'; end
    %if all(mean(imatx') < tol), smeth = ['row' smeth]; end
    %if ~all(mean(imatx) < tol) && ~all(mean(imatx') < tol) && ~strcmp(smeth,'corr'), smeth = 'none'; end

    if findstr(smeth,'corr')
        
    else
        burt = imatx'*imatx;
        
    end

    % this would be a different kind of scaling: normalized co-activations
    % not the same as correlations
    %rmat = repmat(sum(imatx),n,1);
    %cmat = repmat(sum(imatx)',1,n);
    %tmp = burt ./ (  rmat .* cmat ).^.5;

    disp(['Scaling method is ' num2str(smeth)])
    disp(['Decomposition method is ' num2str(decompmethod)])

    MCA.smeth = smeth;
    MCA.smeth_descrip = 'Scaling method for input';
    MCA.decompmethod = decompmethod;

    
    if exist('t','var')
        MCA.t = t;
    end

    if exist('p','var')
        MCA.p = p;
    end

    if exist('corr','var')
        MCA.corr = corr;
    end
    
    switch decompmethod
        case 'pca'

            [pc,score,latent] = princomp(burt);

            % Eigenvalue plot
            figure('Color','w'), set(gca,'FontSize',18),bar(latent),
            xlabel('Components'),ylabel('Variance')

            MCA.pc = pc;
            MCA.coords = score;
        case 'nmds'
            
            if ~exist('D','var'), error('Must have distance matrix D; current smethods for this are phi, corr and correspondence.'); end
            
            if all(all(D >= 0 & abs(D - D') <= 10*eps*max(max(D))))
                % positive, symmetric
            else
                error('Distance matrix is not positive and symmetric.');
            end

            if all(diag(D) == 0)
                % self-distances are zero
            else
                error('Distance matrix diagonals are not zero.');
            end

            %D = simtodis(burt);

            % dims to test: min of nvars ./ 2 or n. obs / 10, but not more than 20
            ntest = min([round(nvars./2) round(nobs./10) 20]);
            % test at least 5 if we have them
            ntest = max(ntest,5);
            ntest = min(ntest,round(nvars./2));
            
            disp(['Testing stress over ' num2str(ntest) ' dimensions.']);

            [Y,obs,imp,stress] = shepardplot(D,ntest);

            MCA.coords = Y;
            MCA.obs = obs;
            MCA.imp = imp;
            MCA.stress = stress;

        otherwise error('Invalid decomposition method.')
    end

    % Component plot
    % figure('Color','w'), set(gca,'FontSize',18)
    % nmdsfig(score(:,1:2),ones(size(pc,1),1),nms);

    % Profile plot

    if isempty(nms), clear nms; for i=1:length(burt),nms{i}=['V' num2str(i)]; end, end

    % only make the plot if we have 2 sets of variables to associate
    if n ~= nvars && ~isempty(n) && n ~= 1
        v = burt_profile(burt,nms,n);  % task x region interaction if double-centered
        if strcmp(smeth,'rowcol')
            title('Task specificity across regions'),ylabel('Marginal study counts')
        elseif strcmp(smeth,'row')
            title('Row-centered profile across regions'),ylabel('Marginal study counts')
        elseif strcmp(smeth,'col')
            title('Column-centered profile across regions'),ylabel('Marginal study counts')
        elseif strcmp(smeth,'none')
            title('Activations by task and brain region'),ylabel('% Studies')
        elseif strcmp(smeth,'corr')
            title('Correlations between tasks and brain regions'),ylabel('Correlation')
        end

        MCA.v = v;
        MCA.v_descrip = 'Normalized values for cross correlations between Set1 and Set2, rows are regions, cols are task conds';
    end
    
    % -----------------------------------------------------------------------
    % testcluster:
    %
    % Permute objects in space, get null hypothesis cluster quality, and test
    % observed cluster solution against this.
    % -----------------------------------------------------------------------
    nclust = min(30,round(nvars./3));  %ntest;
    nclust = max(nclust,3);
    nclust = min(nvars,nclust);
    
    MCA.ndims = size(MCA.coords,2);

    % pval is p-values for each number of clusters tested
    % classes = class ("network") assignments for best clustering solution
    [MCA.ClusterSolution.pvals,MCA.ClusterSolution.classes, MCA.ClusterSolution.classnames ...
        X,MCA.ClusterSolution.include,names]= ...
        testclustnew(MCA.coords,2:nclust,MCA.ndims,niter,nms,'keep','average');

    if ~any(MCA.ClusterSolution.pvals < .05)
        % nothing significant
        MCA.ClusterSolution.classes = ones(size(MCA.ClusterSolution.classes));
    end

    % now get new colors based on superclusters from testclust

    MCA.colors = nmdsfig_tools('class_colors',MCA.ClusterSolution.classes);


    MCA.names = nms;




    % Chi-square: Actual vs. expected
    % for the table as a whole, and for each row
    % do for each sub-matrix (correspondence analysis on each indicator)
    % and for the cross-products (indicator 1 vs. 2)
    %keyboard
    %chi2 = get_chi2(burt)

    %tmp = diag(burt2); e = tmp(1:n) * tmp(n+1:end)'
    %e = e ./ sum(e(:))

    % is there some combination of weights on columns (Regions) that
    % allows us to discriminate among switch types?

    % instead, we could permute the columns independently and get
    % univariate significance levels for each pair of columns
    % but this may be the same as standard correlation?
    % ...and also multivariate sig. for chi2 value or something.

    return

