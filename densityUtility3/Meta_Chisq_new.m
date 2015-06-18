% Multi-mode function for performing voxel-wise chi-squared analysis on
% meta-analysis data (peak activations)
%
% R = Meta_Chisq_new('compute', MC_Setup, ['mask', maskimg])
% cl = Meta_Chisq_new('write', R);
%
% Note: Meta_Chisq works with list of image names for all study maps
% Meta_Chisq_new works with output of 'setup' function in
% Meta_Activation_FWE, and does not require writing study contrast maps as
% image files (all info is contained in MC_Setup struct)
%
% new: October 2008: 
% if there are two conditions in the chi-square analysis, this function can
% write maps and return clusters for "A greater than B" effects--areas in 
% which the proportion of activations for the first condition entered is greater
% than for the second one--and vice versa. 
% Get clusters for directional effects by typing:
% [cl, clavsb, clbvsa] = Meta_Chisq_new('write', R);
% clavsb and clbvsa are the "A>B" and "B>A" regions, respectively
% 
% To get plots of activation proportions in each region, use the following
% code (This is done automatically if you have an MC_Info.mat file in the chi2 directory!):
% load MC_Info
% [studybyroi,studybyset] = Meta_cluster_tools('getdata', clavsb, MC_Setup.unweighted_study_data, MC_Setup.volInfo);
% [prop_by_condition,se,num_by_condition,n] = Meta_cluster_tools('count_by_condition',studybyroi, R.Xi, R.w, 1, {'Diff'}, MC_Setup.Xinms, {[0 0 1] [1 1 0]});
%
% [studybyroi,studybyset] = Meta_cluster_tools('getdata', clbvsa, MC_Setup.unweighted_study_data, MC_Setup.volInfo);
% [prop_by_condition,se,num_by_condition,n] = Meta_cluster_tools('count_by_condition',studybyroi, R.Xi, R.w, 1, {'Diff'}, MC_Setup.Xinms, {[0 0 1] [1 1 0]});
%
% cluster_orthviews(clavsb, {[.1 .1 1]}, 'solid');
% cluster_orthviews(clbvsa, {[1 1 0]}, 'solid', 'add');
% 
% Then you can use the point slice plotter to visualize actual peaks:
% Meta_interactive_point_slice_plot(MC_Setup, DB);
%
% tor wager
%
%
% Example:
% ------------------------------------------------------
% Run an analysis with Meta_Chisq_new, and then use these tools to make bar
% plots of regions. The lines below run the entire analysis.
% mask = which('Activation_FWE_all.img')
% R = Meta_Chisq_new('compute', MC_Setup, 'mask', mask);
% R = Meta_Chisq_new('write', R);
% [studybyroi, studybyset] = Meta_cluster_tools('getdata', cl(1), R.dat, R.volInfo);
% [prop_by_condition, se, num_by_condition, n] = Meta_cluster_tools('count_by_condition', studybyset, R.Xi, R.w, 1);
% set(gca, 'XTickLabel', MC_Setup.Xinms)
% ------------------------------------------------------
%
% Example 2: start from the beginning, all steps
% DB = Meta_Setup(DB);
% mask = DB.maskname;
% MC_Setup = Meta_Activation_FWE('setup', DB);
% % ***note: You must 'add a contrast' and select the studies you want to
% % do chi-sq analysis on, and some arbitrary contrast over those studies.
% % this function will use the MC_Setup.Xi field, which has indicators for
% % the different categories of study you're analyzing with chi-sq.
% R = Meta_Chisq_new('compute', MC_Setup, 'mask', mask); save Meta_Chisq_R R

% Programmers' notes
% cosmetic edits/minor edits, tor wager, july 2011

function varargout = Meta_Chisq_new(meth, varargin)
    R = [];
    maskimg = 'Activation_FWE_all.img';

    spm_defaults;
    
    for i = 1:length(varargin)
        if isstruct(varargin{i})
            switch meth
                case 'compute'
                    MC_Setup = varargin{i};
                case 'write'
                    R = varargin{i};
                otherwise
                    error('Unknown method.  Could create problems!!')
            end
        elseif ischar(varargin{i})
            switch varargin{i}

                % functional commands
                case 'mask', maskimg = varargin{i+1}; varargin{i+1} = [];
                case 'data', dat = varargin{i+1};

                case 'MC_Setup', MC_Setup = varargin{i + 1}; varargin{i + 1} = [];
                    
                otherwise, warning(['Unrecognized input string option:' varargin{i}]);
            end
        end
    end

    switch meth
        case 'compute'
            % do first just to make sure no errs after running!
            Rtmp = [];
            Rtmp.contrastname = input('Enter contrast name for this analysis (no spaces or special chars) ', 's');


            disp('Preparing data.')
            if ~exist('dat', 'var'), dat = MC_Setup.unweighted_study_data; end

            % Get only in-mask voxels
            volInfo = iimg_read_img(maskimg, 2);
            is_significant = volInfo.image_indx(MC_Setup.volInfo.wh_inmask);

            dat = full(double(dat(is_significant,:)));

            fprintf(' %3.0f valid voxels.\n', size(dat, 1));
            
            disp('Preparing task indicators and weights.')
            % task and weights
            Xi = MC_Setup.Xi;

            w = MC_Setup.wts;
            w = w ./ mean(w);

            [nvox, nstudies] = size(dat);
            ntasks = size(Xi, 2);

            disp('Getting proportion activations for each condition.')
            % get proportions for the entire matrix
            [icon, ctxtxi, betas, num_by_condition, prop_by_condition] = meta_apply_contrast(dat, Xi, w, MC_Setup.contrasts);

            disp('Running chi-square with nonparametric options for each in-mask voxel.')
            % run nonparametric chi2 test
            [chi2, df, p, sig, warn, tab, expected, isnonpar] = meta_analyze_data(dat', 'X', Xi, 'chi2', 'w', w);

            R = struct('volInfo', volInfo, 'in_mask', is_significant, 'contrastname', Rtmp.contrastname, ...
                'dat', sparse(logical(dat)), 'Xi', Xi, 'w', w, 'prop_by_condition', prop_by_condition, ...
                'num_by_condition', num_by_condition, 'chi2', chi2, 'df', df, 'p', p, 'sig', sig, 'warn', warn, ...
                'isnonpar', isnonpar);

            varargout{1} = R;

        case 'write'
            disp('Meta_Chisq: Getting and saving results');
            pvals = [.05 .01 .005 .001 Inf];
            str = cell(length(pvals), 1);
            sig_vector = zeros(length(R.p), length(pvals));
            for i = 1:length(pvals)
                [nsig, str{i}, sig_vector(:,i)] = get_sig_vector(pvals(i), R.p);

                fprintf(1, '%3.0f\tP < %3.4f %s\t%3.0f\n', i, pvals(i), str{i}, nsig);
            end
            fprintf(1, '\n');

            % tor added july 2011 to fix bug with out-of-mask voxels and sizes
            sig_vector = zeroinsert(~R.in_mask, sig_vector);
            R.chi2 = zeroinsert(~R.in_mask, R.chi2')';
            R.prop_by_condition = zeroinsert(~R.in_mask, R.prop_by_condition);
            
            wh = input('Enter index of map to write: ');

            % p-value image
            name = [R.contrastname '_' num2str(pvals(wh)) '_' str{wh} '_chi2p'];
            name(name == '.') = '_';
            name = [name '.img'];
            iimg_reconstruct_vols(sig_vector(:,wh), R.volInfo, 'outname', name, 'descrip', 'Created by Meta_Chisq_new.m');
            fprintf(1, 'Written: %s\n', name);

            % chi2 image
            name = [R.contrastname '_' num2str(pvals(wh)) '_' str{wh} '_chi2'];
            name(name == '.') = '_';
            name = [name '.img'];
            chi2dat = R.chi2' .* (sig_vector(:,wh) > 0);
            iimg_reconstruct_vols(chi2dat, R.volInfo, 'outname', name, 'descrip', 'Created by Meta_Chisq_new.m');
            fprintf(1, 'Written: %s\n', name);

            % save chi-square values in clusters
            cl = mask2clusters(name);
            varargout{1} = cl;

            if size(R.prop_by_condition, 2) == 2
                % two conditions. Write pos and neg results
                avsb = -1 * diff(R.prop_by_condition')';

                % A > B
                chi2dat = R.chi2' .* (sig_vector(:,wh) > 0) .* (avsb > 0);
                name = [R.contrastname '_AgreaterthanB_' num2str(pvals(wh)) '_' str{wh} '_chi2'];
                name(name == '.') = '_';
                name = [name '.img'];

                iimg_reconstruct_vols(chi2dat, R.volInfo, 'outname', name, 'descrip', 'Created by Meta_Chisq_new.m');
                clavsb = mask2clusters(name);
                varargout{2} = clavsb;
                fprintf(1, 'Written: %s\n', name);

                % B > A
                chi2dat = R.chi2' .* (sig_vector(:,wh) > 0) .* (avsb < 0);
                name = [R.contrastname '_BgreaterthanA_' num2str(pvals(wh)) '_' str{wh} '_chi2'];
                name(name == '.') = '_';
                name = [name '.img'];

                iimg_reconstruct_vols(chi2dat, R.volInfo, 'outname', name, 'descrip', 'Created by Meta_Chisq_new.m');
                clbvsa = mask2clusters(name);
                varargout{3} = clbvsa;
                fprintf(1, 'Written: %s\n', name);

                % Extent threshold
                if ~isempty(clbvsa), clbvsa(cat(1, clbvsa.numVox) < 3) = []; end
                if ~isempty(clavsb), clavsb(cat(1, clavsb.numVox) < 3) = []; end

                % Make plots
                if exist('MC_Setup', 'var')
                elseif exist(fullfile(pwd, 'MC_Info.mat'), 'file')
                    load MC_Info MC_Setup
                end
                
                if exist('MC_Setup', 'var')
                    [studybyroi,studybyset] = Meta_cluster_tools('getdata', clavsb, MC_Setup.unweighted_study_data, MC_Setup.volInfo);
                    [prop_by_condition,se,num_by_condition,n] = Meta_cluster_tools('count_by_condition',studybyroi, R.Xi, R.w, 1, {'Diff'}, MC_Setup.Xinms, {[0 0 1] [1 1 0]});

                    [studybyroi,studybyset] = Meta_cluster_tools('getdata', clbvsa, MC_Setup.unweighted_study_data, MC_Setup.volInfo);
                    [prop_by_condition,se,num_by_condition,n] = Meta_cluster_tools('count_by_condition',studybyroi, R.Xi, R.w, 1, {'Diff'}, MC_Setup.Xinms, {[0 0 1] [1 1 0]});

                    cluster_orthviews(clavsb, {[.1 .1 1]}, 'solid');
                    cluster_orthviews(clbvsa, {[1 1 0]}, 'solid', 'add');
                else
                    disp('MC_Info.mat file not found in the current directory. Skipping plots.');
                    cluster_orthviews(clavsb, {[.1 .1 1]}, 'solid');
                    cluster_orthviews(clbvsa, {[1 1 0]}, 'solid', 'add');
                end

            else
                varargout{2} = [];
                varargout{3} = [];

                disp('More than two conditions: Omitting writing of results images containing directional effects.');

            end
                
        otherwise
            error('Unknown method.')
    end
    
end % Main function


% sub-functions

function [nsig, str, sig_vector] = get_sig_vector(p, chi2p)
    if isinf(p)
        str = 'FDR';
        p = FDR(chi2p, .05);
        if isempty(p), p = -Inf; end
    else
        str = 'Unc';
    end
    sig_vector = (chi2p .* (chi2p < p))';
    nsig = sum(sig_vector > 0);
end

