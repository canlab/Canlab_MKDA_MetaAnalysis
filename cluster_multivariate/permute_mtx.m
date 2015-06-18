function [OUT] = permute_mtx(data,varargin)
    % [OUT] = permute_mtx(data,[meth],[verbose],[niter],[separator])
    %
    % tor wager
    %
    % PERMUTATION TEST FOR STOCHASTIC ASSOCIATION BETWEEN COLUMNS
    % OF DATA (Ho)
    %
    % This function takes as input two matrices,
    % an actual matrix of data
    % and an expected / null hypothesis
    % covariance or correlation matrix.
    %
    % the algorithm provides significance levels for
    % the cov / correlation of the columns of data,
    % based on the expected correlation e
    %
    % it computes a test statistic, which is the
    % average squared deviation from the expected
    % values, element by element, over the matrix.
    % This statistic provides an "omnibus" test
    % of whether there are deviations from the
    % expected values in the matrix.
    %
    % for each element, a statistic is also computed -
    % the squared deviation from the expected -
    % which provides a test of significance for each
    % element.
    %
    % a permutation test is used to create an Ho
    % distribution.  Columns of the matrix are
    % randomly permuted (independently), and the test
    % statistics assessed over iterations.
    %
    % R2 is average squared off-diagonal value in
    % lower triangular matrix
    %
    % called in db_cluster_burt_table.m
    %
    % example:
    % data = bmatx1(:,1:5); ss = data'*data; e = diag(diag(ss));
    % OUT = permute_mtx(bmatx1(:,1:5),e,'ss',1,1000);
    %
    % data = bmatx1(:,length(fnames)+1:end); ss = data'*data; e = diag(diag(ss));
    % OUT = permute_mtx(data,'ss',1,1000);
    %
    % If we scale rows of data so that they sum to 0, we can remove
    % the effects of some rows having higher values overall than others
    % Since columns are permuted, there's no problem with those in the stats,
    % but we may want to scale to make visualization more interpretable (less
    % misleading)
    %
    % for Multiple Correspondence Analysis (MCA)
    % enter a Separator value as the last argument
    % - this is an integer that tells it how to partition the matrix
    % - indicates number of columns in primary partition
    % - works only for two variables, with 2 sets of indicator columns
    % - tests omnibus, 1st part, 2nd part, and covariance between 1 and 2
    %
    % Example:
    % bmatx1 is a dataset with 2 sets of indicator variables
    %
    % data = bmatx1; ss = data'*data; e = diag(diag(ss));
    % OUT = permute_mtx(data,'ss',1,1000,length(fnames));

    if length(varargin) > 0, meth = varargin{1};, else, meth = 'corr';, end
    if length(varargin) > 1, verbose = varargin{2};, else, verbose = 1;, end
    if length(varargin) > 2, niter = varargin{3};, else, niter = 5000;, end

    switch meth
        case 'cov', m = cov(data);
        case 'ss', m = data'*data;
        case 'corr', m = corrcoef(data);
        otherwise, error('Unknown method!  OK methods are cov, ss, and corr')
    end

    if length(varargin) > 3
% -----------------------------------------------------------------------        
        % Data is partitioned into two sets; run on each set
% -----------------------------------------------------------------------
        sep = varargin{4};

        disp(' ');disp('RESULTS FOR FIRST PARTITION')
        OUT.part1 = do_permute(m(1:sep,1:sep),data(:,1:sep),meth,verbose,niter);
        disp(' ');disp('RESULTS FOR SECOND PARTITION')
        OUT.part2 = do_permute(m(sep+1:end,sep+1:end),data(:,sep+1:end),meth,verbose,niter);
        disp(' ');disp('RESULTS FOR CROSS-CORRESPONDENCE')
        OUT.part12 = do_permute(m,data,meth,verbose,niter,sep);

        OUT.sep = sep;
        OUT.sep_descrip = 'Last column in first data partition.';
        OUT.p = combine_parts(OUT,'p2');
        OUT.p_descrip = 'p-vals for correspondence between each pair in superindicator';

        OUT.sign_matrix = sign(combine_parts(OUT,'actual_vs_expected'));
        OUT.sig05 = OUT.p < .05 .* OUT.sign_matrix;
        [OUT.fdr_pthr,OUT.sigfdr] = fdr_correct_pvals(OUT.p,OUT.sign_matrix);

    else
% -----------------------------------------------------------------------        
        % No data partition
% -----------------------------------------------------------------------
        % for regular MCA, do it with the whole matrix
        OUT = do_permute(m,data,meth,verbose,niter);

    end


    return





function OUT = do_permute(m,data,meth,verbose,niter,varargin)
    % if varargin is entered, this function takes the cross-product part of MCA
    % where the first set of vars is defined by the integer n in varargin{1}
    %
    % permutes only the SECOND set of columns, after sep, to preserve the
    % row structure of the first set, assumed to represent tasks.
    % the second set, after sep, is assumed to represent brain regions

    if verbose,
        fprintf(1,'\nPermute_mtx.m\n--------------------------\n')
        fprintf(1,'Test that columns of data have no systematic relationship\n')
        fprintf(1,'\nActual %s matrix for data',meth),m,
    end


    % mask: which elements of cov matrix to use
    %
    % for MCA: take cov ones
    if length(varargin) > 0
        % -----------------------------------------------------------------------        
        % Data is partitioned into two sets; Entering sep means we're
        % asking for the correlations/covs between sets only
        % -----------------------------------------------------------------------
        sep = varargin{1};      % separator
        mm = ones(size(m));     % full mtx size
        rows = 1:sep;
        nrows = sep;
        cols = sep+1:size(m,1);
        ncols = length(cols);
        mm(rows,cols) = 0; maskind = find(mm==0);
        disp('Mask matrix: Elements with zeros are used in computation of stats'),mm
    else
        % -----------------------------------------------------------------------        
        % Setup for both separate partitions A and B, and single matrix
        % -----------------------------------------------------------------------
        %regular correspondence
        % get index of lower triangular matrix
        % to compute stats only on these values.
        rows = 1:size(m,1);
        cols = rows;
        disp('Using lower triangular matrix for stats computation')
        mask = triu(Inf*eye(size(m,2))); maskind = find(mask==0);
    end

    % -----------------------------------------------------------------------        
        % s1 and s2 are statistics for the correct permutation
    % -----------------------------------------------------------------------
    [s1,s2] = getstats(m,maskind);

    % -----------------------------------------------------------------------        
        % Iterate and get null hypothesis stats
    % -----------------------------------------------------------------------
    for i = 1:niter

        if mod(i,1000) == 1, fprintf(1,'.');, end

        if length(varargin) > 0
            [mi] = getcov(data,meth,sep);               % get Ho realization permuting columns AFTER sep index
        else
            [mi] = getcov(data,meth);                     % get Ho realization
        end
        [s1n(i),s2n(i,:)] = getstats(mi,maskind);     % Ho R-square values
        meanmi(:,:,i) = mi;

    end
    fprintf(1,'\n')
    meanmi = mean(meanmi,3);

    % -----------------------------------------------------------------------        
        % Get p-values
    % -----------------------------------------------------------------------
    p = 2 * min(sum(s1 >= s1n) ./ niter,sum(s1 <= s1n) ./ niter);
    for i = 1:length(s2), p2(i) = sum(s2(i) < s2n(:,i))./ niter;, end
    p2c = p2 .* length(p2);

    % -----------------------------------------------------------------------        
        % Print output
    % -----------------------------------------------------------------------
    if verbose
        fprintf(1,'\nExpected (Mean Ho) %s matrix\n',meth),meanmi
        fprintf(1,'\nOmnibus test for differences from expected on %s\n',meth)
        fprintf(1,'Obs. R2 %3.3f, Expected (Ho) R2 = %3.3f, p = %3.4f\n',s1,mean(s1n),p)
        %fprintf(1,'\nSignificant individual tests for differences from expected on %s\n',meth)
        sig = find(p2 <= .05);
        %for i = 1:length(sig)
        %[row,col] = ind2sub(size(m),maskind(sig(i)));
        %fprintf(1,'[%3.0f,%3.0f], R2 = %3.3f, Ho R2 = %3.3f, p = %3.4f, bonf_p = %3.4f\n',row,col,s2(sig(i)),nanmean(s2n(:,sig(i))),p2(sig(i)),p2c(sig(i)))
        %figure; hist(s2n(:,sig(i)))
        %end
        if isempty(sig), disp('No significant results.'), end
    end

    % -----------------------------------------------------------------------        
        % Save output vars in structure
    % -----------------------------------------------------------------------
    OUT.actual_vs_expected = m(rows,cols) - meanmi(rows,cols);
    if strcmp('meth','ss')
        signvals = sign(OUT.actual_vs_expected);
    else
        signvals = sign(m(rows,cols));
    end

    OUT.m = m; OUT.s1 = s1; OUT.s2 = s2;
    OUT.p = p; OUT.p2 = p2; OUT.p2c = p2c;OUT.meanmi = meanmi;
    OUT.sig = sig; OUT.niter = niter; OUT.meth = meth; OUT.s1n = s1n; OUT.s2n = s2n;

    % add fnames, and fnames for 2nd set!

    % -----------------------------------------------------------------------        
        % Get indicator matrices for corrected/uncorrected significance
    % -----------------------------------------------------------------------
    % matrices of significant results
    if exist('sep','var')

        OUT.sig05 = reshape(OUT.p2 <= .05,nrows,ncols) .* signvals;

        OUT.bonf_pthr = .05 ./ length(OUT.p2);
        OUT.sigbonf = reshape(OUT.p2 <= OUT.bonf_pthr,nrows,ncols) .* signvals;

        [OUT.fdr_pthr,OUT.sigfdr] = fdr_correct_pvals(OUT.p2,signvals(:)');
        OUT.sigfdr = reshape(OUT.sigfdr,nrows,ncols) .* signvals;

    else
        OUT.sig05 = (squareform(OUT.p2) <= .05);
        OUT.sig05 = double(OUT.sig05 - eye(size(OUT.sig05)) .* signvals);

        OUT.sigbonf = (squareform(OUT.p2) <= (.05 ./ length(OUT.p2)));
        OUT.sigbonf = double(OUT.sigbonf - eye(size(OUT.sigbonf)) .* signvals);
    end

    return





function [m] = getcov(data,meth,varargin)

    % this would be a bootstrap
    % resample data - tested to evenly sample rows
    %len = size(data,1);
    %wh = ceil(rand(len,1) .* len);
    %d2 = data(wh,:);

    startat = 1;
    if length(varargin) > 0
        startat = varargin{1}+1;
        d2(:,1:varargin{1}) = data(:,1:varargin{1});
    end

    % permute columns
    for i = startat:size(data,2)
        d2(:,i) = getRandom(data(:,i));
    end

    switch meth
        case 'cov', m = cov(d2);
        case 'ss', m = d2'*d2;
        case 'corr', m = corrcoef(d2);
        otherwise, error('Unknown method!  OK methods are cov, ss, and corr')
    end

    return



function [s1,s2] = getstats(m,maskind)
    % input is cov matrix of whatever form - corr, ss, cov

    % this is for MCA, with separator integer indicating last col. in 1st set
    %if length(varargin) > 0,
    %    sep = varargin{1};,
    %    m2 = m(1:sep,sep+1:end);
    %    s2 = (m2 .^ 2)';
    %    s1 = nanmean(ssep2);
    %else

    m = m(maskind); m = m(:);

    s2 = (m .^ 2)';
    s1 = nanmean(s2);
    %end

    return


function [s1,s2] = getvec(m,maskind)
    % input is cov matrix of whatever form - corr, ss, cov

    m = m(maskind); m = m(:);

    s2 = m';
    s1 = nanmean(s2);

    return



function [pthr,sig] = fdr_correct_pvals(p,r)

    % if p is a square correlation matrix, vectorize
    if diff(size(p)) == 0
        psq = p; psq(find(eye(size(p,1)))) = 0;
        psq = squareform(psq);
    end
    pthr = FDR(p,.05);
    if isempty(pthr), pthr = 0; end

    sig = sign(r) .* (p < pthr);

    return



function outmtx = combine_parts(OUT,myfield)
    nrows = length(1:OUT.sep);
    ncols = size(OUT.part2.m,2);

    a = OUT.part1.(myfield);
    if any(size(a) - nrows)
        a = squareform(a);
    end

    b = OUT.part2.(myfield);
    if any(size(b) - ncols)
        b = squareform(b);
    end

    c = OUT.part12.(myfield);
    if any(size(c) - [nrows ncols])
        c = reshape(c,nrows,ncols);
    end
    outmtx = [a c; c' b];

    return
