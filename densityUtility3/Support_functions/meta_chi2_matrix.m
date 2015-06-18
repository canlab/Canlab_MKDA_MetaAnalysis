function out = meta_chi2_matrix(dat,w)

w = w ./ mean(w);
myalpha = .05;

 [N,npairs] = size(dat);


    [rows,cols,ncorr] = corrcoef_indices(npairs);
    
    str = sprintf('Computing differences among correlations %04d',0); fprintf(1,str);

    chi2 = zeros(ncorr,1);
    chi2p = zeros(ncorr,1);
    diffr = zeros(ncorr,1);
    asym = zeros(ncorr,1);
    %sig = zeros(ncorr,2);

    for cc = 1:ncorr

        fprintf(1,'\b\b\b\b%04d',cc);

        d = [dat(:,rows(cc)) dat(:,cols(cc))];
        [chi2(cc),df,chi2p(cc),sig,warn,tab,e] = chi2test(d,'obs',w,1);

        % weighted # activations for each column
        rsum = sum(tab);
        csum = sum(tab,2);
        yesses = [csum(2) rsum(2)];

        % signed measure of association; diff from expected
        diffr(cc) = (tab(2,2) - e(2,2)) ./ (min(yesses) - e(2,2));  
        
        % asymmetry: 
        a = [tab(2,1) tab(1,2)] ./ yesses;     % proportion unique (nonshared) in each col
        asym(cc) = diff(a) ./ sum(a);
    end

    erase_string(str);


    chi2 = reconstruct(chi2,npairs,ncorr,rows,cols);
    p = reconstruct(chi2p,npairs,ncorr,rows,cols);

    diffr = reconstruct(diffr,npairs,ncorr,rows,cols);
    asym = reconstruct(asym,npairs,ncorr,rows,cols);
    
    sig = (p <= myalpha - eye(size(p))) .* sign(diffr);

    % FDR corrected
    pthr = FDR(chi2p,.05);
    if isempty(pthr), pthr = 0; end

    sigfdr = (p <= pthr) .* sign(diffr);

    out = struct('alpha',myalpha,'diffr',diffr,'asym',asym,'chi2',chi2, ...
        'p',p,'sig',sig,'pthr',pthr,'sigfdr',sigfdr);
    
    
    
    
    
function [rows,cols,ncorr] = corrcoef_indices(npairs)
    % upper triangle only
    tmp = triu(ones(npairs));
    tmp = tmp - eye(npairs);
    [rows,cols] = find(tmp);
    ncorr = length(rows);
    return
    
    
    
function valmat = reconstruct(vals,npairs,ncorr,rows,cols)

    valmat = zeros(npairs);
    for i = 1:ncorr
        valmat(rows(i),cols(i)) = vals(i);
    end
    valmat = valmat + valmat';

    return



function erase_string(str1)
    fprintf(1,repmat('\b',1,length(str1))); % erase string
    return
