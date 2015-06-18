function plot(metaobj, plotmethod)
% plot(metaobj, [plotmethod])
%
% Plot methods:
% ----------------------------------------
% Plot data matrix
% plot(fmri_data_object)
%
% Plot means by condition
% plot(fmri_data_object, 'means_for_unique_Y')
%
%
if nargin < 2
    plotmethod = 'data';
end


switch plotmethod
    %  ==============================================================
    case 'data'
        %  ==============================================================
        
        if isempty(metaobj.dat)
            warning('No data in .dat field.');
            return
        end
        
        create_figure('fmri data matrix', 2, 2);
        imagesc(metaobj.dat');
        colorbar; colormap gray
        axis tight; set(gca, 'YDir', 'Reverse')
        title('meta-analysis .dat Data matrix');
        xlabel('Voxels'); ylabel('Images');
        drawnow
        
        if ~isempty(metaobj.classes)
            p = get(gca, 'Position'); ystart = p(2); ylen = p(4);
            
            axh = axes('Position', [.05 ystart .03 ylen]);
            imagesc(metaobj.classes);
            title('Classes');
            axis tight;
        end
        drawnow
        
        % Centering
        % datc = scale(scale(metaobj.dat, 1)', 1)';
        % datc = scale(scale(metaobj.dat)')';
        
        
        % ---------------------------------------------------------------
        % Histogram
        % ---------------------------------------------------------------
        dattmp = metaobj.dat(:);
        subplot(2, 2, 2);
        [h, x] = hist(dattmp, 100);
        han = bar(x, h);
        set(han, 'FaceColor', [.3 .3 .3], 'EdgeColor', 'none');
        axis tight;
        xlabel('Values'); ylabel('Frequency');
        title('Histogram of values');
        drawnow
        
        clear dattmp
        
        
        % ---------------------------------------------------------------
        % Covariance
        % ---------------------------------------------------------------
        covmtx = cov(metaobj.dat);
        subplot(2, 2, 3);
        imagesc(covmtx);
        axis tight; set(gca, 'YDir', 'Reverse')
        title('cov(dat''), cov(rows of dat)');
        colorbar
        drawnow
        
        if ~isempty(metaobj.classes)
            p = get(gca, 'Position'); ystart = p(2); ylen = p(4);
            
            axh = axes('Position', [.05 ystart .03 ylen]);
            imagesc(metaobj.classes);
            title('Y');
            axis tight;
            
        end
        drawnow
        
        subplot(2, 2, 4);
        
        globalmean = nanmean(metaobj.dat);  % global mean of each obs
        globalstd = nanstd(metaobj.dat);  % global mean of each obs
        
        nobs = length(globalmean);
        
        Y = metaobj.classes;
        Yname = 'Y values in fmri data obj';
        if isempty(Y)
            Y = 1:nobs;
            Yname = 'Case number';
        end
        
        sz = rescale_range(globalstd, [4 14]); % marker size related to global std
        sz(sz < .5) = .5;
        
        for i = 1:nobs
            plot(Y(i), globalmean(i), 'ko', 'MarkerSize', sz(i), 'LineWidth', 1);
        end
        ylabel('Global mean');
        xlabel(Yname);
        title('Globals for each case (size = spatial std)')
        drawnow
        
        % [coeff, score, latent] = princomp(metaobj.dat, 'econ');
        % %d2 = mahal(score, score);
        % plot(latent)
        
        
        % ---------------------------------------------------------------
        % Orthviews
        % ---------------------------------------------------------------
        % check to be sure:
        metaobj.dat(isnan(metaobj.dat)) = 0;
        
        m = mean(metaobj.dat')';
        s = std(metaobj.dat')';
        d = m./s;  
        d(m == 0 | s == 0) = 0;
        vecs_to_reconstruct = [m s d];
 
        if isempty(metaobj.volInfo)
            disp('.volInfo is empty. Skipping orthviews and other brain plots.');
        else
            create_orthviews(vecs_to_reconstruct, metaobj);
            spm_orthviews_name_axis('Mean data', 1);
            spm_orthviews_name_axis('STD of data', 2);
            spm_orthviews_name_axis('Mean / STD', 3);
            set(gcf, 'Name', 'Orthviews_fmri_data_mean_and_std');
        end
        
        %  ==============================================================
    case 'means_for_unique_Y'
        %  ==============================================================
        
        u = unique(metaobj.classes);
        
        [v, n] = size(metaobj.dat);
        nu = length(u);
        
        if nu > 20
            error('More than 20 unique values of Y.  For means_by_condition, Y should be discrete integer-valued.');
        end
        
        [means, stds] = deal(zeros(nu, v));
        
        for i = 1:nu
            means(i, :) = nanmean(metaobj.dat(:, metaobj.classes == u(i))');
            stds(i, :) = nanstd(metaobj.dat(:, metaobj.classes == u(i))');
        end
        
        create_figure('means by condition (unique Y values)', 2, 1);
        imagesc(means);
        colorbar
        axis tight; set(gca, 'YDir', 'Reverse')
        title('Means by condition');
        xlabel('Voxels');
        if iscell(metaobj.classes_names) && ~isempty(metaobj.classes_names)
            set(gca, 'YTick', u, 'YTickLabel', metaobj.classes_names);
        else
            ylabel('Unique Y values');
        end
        
        drawnow
        
        subplot(2, 1, 2)
        imagesc(stds);
        colorbar
        axis tight; set(gca, 'YDir', 'Reverse')
        title('Standard deviations by condition');
        xlabel('Voxels');
        if iscell(metaobj.classes_names) && ~isempty(metaobj.classes_names)
            set(gca, 'YTick', u, 'YTickLabel', metaobj.classes_names);
        else
            ylabel('Unique Y values');
        end
        drawnow
        
        % ---------------------------------------------------------------
        % Orthviews
        % ---------------------------------------------------------------
        if ~isempty(metaobj.volInfo)
            vecs_to_reconstruct = means';
            create_orthviews(vecs_to_reconstruct, metaobj);
            n = size(vecs_to_reconstruct, 2);
            for i = 1:n
                
                if iscell(metaobj.classes_names) && ~isempty(metaobj.classes_names)
                    spm_orthviews_name_axis(metaobj.classes_names{i}, i);
                end
                
            end
            set(gcf, 'Name', 'Orthviews_means_by_unique_Y');
            
            
            % ---------------------------------------------------------------
            % Montage: variance across conditions
            % ---------------------------------------------------------------
            vecs_to_reconstruct = std(means)';
            fig_handle = create_montage(vecs_to_reconstruct, metaobj);
            set(fig_handle, 'Name', 'Montage_cariability_across_conditions')
        end
        
        
    otherwise
        error('Unknown plot method');
end

end


function create_orthviews(vecs_to_reconstruct, metaobj)

n = size(vecs_to_reconstruct, 2);
overlay = which('SPM8_colin27T1_seg.img');
spm_check_registration(repmat(overlay, n, 1));

for i = 1:n
    
    cl{i} = iimg_indx2clusters(vecs_to_reconstruct(:, i), metaobj.volInfo);
    cluster_orthviews(cl{i}, 'add', 'handle', i);
    
    spm_orthviews_change_colormap([0 0 1], [1 1 0], [0 1 1], [.5 .5 .5], [1 .5 0]);
    
end

end


function fig_handle = create_montage(vecs_to_reconstruct, metaobj)

n = size(vecs_to_reconstruct, 2);
overlay = which('SPM8_colin27T1_seg.img');

for i = 1:n
    
    dat = vecs_to_reconstruct(:, i);
    % top and bottom 10%
    dat(dat > prctile(dat, 10) & dat < prctile(dat, 90)) = 0;
    
    cl{i} = iimg_indx2clusters(dat, metaobj.volInfo);
    
    fig_handle(i) = montage_clusters(overlay, cl{i}, [2 2]);
    
    set(fig_handle, 'Name', sprintf('Montage %3.0f', i), 'Tag', sprintf('Montage %3.0f', i))
    
end

end


function rx = rescale_range(x, y)
% re-scale x to range of y
m = range(y)./range(x);

if isinf(m)
    % no range/do not rescale
    rx = x;
else
    b = y(1) - m * x(1);
    rx = m*x + b;
end
end


