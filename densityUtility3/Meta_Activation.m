function DB = Meta_Activation(DB)
% DB = Meta_Activation(DB)
% 
% NEEDS:
%DB.PP       % list of contrast image names
%DB.XYZmm    % mm coords for each study
%DB.radius   % in voxels

spm_defaults; defaults.analyze.flip = 0;

% -----------------------------------------------------
% across all density images, compute prob(activation), p(a)
% -----------------------------------------------------
fprintf(1,'Assuming images are not yet weighted by studyweights.\n');
fprintf(1,'Creating image of activation probability.\n'); tic

% Old code: When weighting was implemented in Meta_Setup.
% Now weighting is implmented in this step, with weighted_sum_image.m
%warning off, v = spm_read_vols(spm_vol(DB.PP)); warning on
%tor_spm_mean_ui(DB.PP,'Activation.img');
%V = spm_vol('Activation.img'); v = spm_read_vols(V);
%v = v .* length(DB.studyweight);   % go from average to weighted sum
%spm_write_vol(V,v); pa = v;

weighted_sum_image(DB.PP,'Activation.img',DB.studyweight);
V = spm_vol('Activation.img'); pa = spm_read_vols(V);


toc,fprintf(1,'\n');

%pa = sum(v,4);

%V = DB.maskV; V.fname = 'Activation.img'; spm_write_vol(V,pa);


% -----------------------------------------------------
% get expected null-hypothesis probability
% uncorrected version
% -----------------------------------------------------
fprintf(1,'Test of overall activation.\n'); tic
[p0,pa_h0] = p0_monte_carlo(DB);
toc,fprintf(1,'\n');
DB.pa_h0 = pa_h0;

% -----------------------------------------------------
% load count image for absolute number of studies threshold
% -----------------------------------------------------
Vc = spm_vol('Activation_counts.img'); cdat = spm_read_vols(Vc);

% -----------------------------------------------------
% Threshold (uncorrected) and print results
% -----------------------------------------------------

pa(pa <= p0) = 0;       % threshold based on Monte Carlo p-value
pa(cdat < 2) = 0;       % threshold additionally based on counts (must have 2 contrasts)

V.fname = 'Activation_thresholded.img'; spm_write_vol(V,pa);

% -----------------------------------------------------
% get p-map and FDR corrected maps, display results
% -----------------------------------------------------

numvox = pa>0; numvox = sum(numvox(:));
fprintf(1,'Uncorrected p < .05 and at least 2 contrasts: %3.0f voxels.\n',numvox);
if numvox < 40000
    % Do FDR correction here
    [cl,ptname,pname,pt] = dens_fdr_thresh(V.fname,pa_h0,1,p0,cdat);
    cluster_orthviews(cl,{[1 0 0]});
    fprintf(1,'FDR corrected: threshold = p < %3.4f, %3.0f voxels.\n',pt,sum(cat(1,cl.numVox)));
    
else
    fprintf(1,'Too many voxels to display meaningfully.\n');
end


docorr = 0;
if docorr
    % -----------------------------------------------------
    % get expected null-hypothesis probability
    % Corrected for whole-brain search
    % -----------------------------------------------------

    fprintf(1,'Brain-wise FWER correction.\n'); tic
    [p0fwe,pa_h0fwe] = p0_monte_carlo_FWE(DB,100);
    DB.pa_h0fwe = pa_h0fwe;
    toc,fprintf(1,'\n');

    % -----------------------------------------------------
    % threshold activation map
    % -----------------------------------------------------
    
    pa2 = pa;
    pa2(pa2 <= p0fwe) = 0;
    V.fname = 'Activation_thresholded_FWE.img'; spm_write_vol(V,pa2);
    
    % -----------------------------------------------------
    % display results
    % -----------------------------------------------------
    
    numvox = pa2>0; numvox = sum(numvox(:));
    fprintf(1,'FWE corrected: threshold = p < %3.4f, %3.0f voxels.\n',pt,sum(cat(1,cl.numVox)));
    
    clfwe = mask2clusters(V.fname);
    
    if ~isempty(clfwe)
        save Activation_FWE_clusters clfwe
        cluster_orthviews(clfwe,{[1 0 0]});
        try, 
            DB2 = DB; clear DB; global DB, DB = DB2;, clear DB2
            Meta_interactive_table;
            set(gcf,'WindowButtonUpFcn','Meta_interactive_table;')
       catch, disp('Cannot display tables.'), 
       end
       
    else
        fprintf(1,'No results to display.\n');
    end

end


return



% -----------------------------------------------------
% get variance and p-values
% -----------------------------------------------------
%varp0 = sum(DB.studyweight.^2) * p0 * (1-p0)

%Z = pa ./ sqrt(varp0);
%p = 2*(1-normcdf(Z));

% problem: assumptions are violated.  successes/failures.  

