%Meta_Activation_FWE(meth)
%
% Meta_Activation_FWE('all', DB, iterations)   [default]
% > Run all of the commands below in sequence
%
% MC_Setup = Meta_Activation_FWE('setup', DB)
% > Create activation map .img and MC_Setup structure; save MC_Info
%
% Meta_Activation_FWE('mc', iterations)
% > Create or add null hypothesis iterations to MC_Info
%
% Meta_Activation_FWE('results')
% > Get and plot results
% Examples:
% For contrasts, positive con values.  The 1 below indicates to make plots.
% cl = Meta_Activation_FWE('results', 1, 'poscon')
% Negative con values
% cl = Meta_Activation_FWE('results', 1, 'negcon')
% Which thresholds to use:
% Meta_Activation_FWE('results', 1, 'height', 'stringent', 'medium', 'lenient') [default]
% Meta_Activation_FWE('results', 1, 'height', 'stringent') <- only use height and most stringent extent
% Meta_Activation_FWE('results', 1, 'negcon', 'height', 'stringent') <- only use height and most stringent extent and neg con values
% Which contrast to use:
% Meta_Activation_FWE('results', 1, 'negcon', 'contrast', 3) <- negative con values for contrast 3
%
% Enter colors:
% Meta_Activation_FWE('results','colors',[1 1 0;1 .5 0;1 .3 .3;.8 .4 .4])
% Meta_Activation_FWE('results','grayscale');
%
% Plot surfaces as well:
% Meta_Activation_FWE('results','surface');
%
% Save slices showing each cluster center:
% Meta_Activation_FWE('results','slices');
% Meta_Activation_FWE('results', 1, 'negcon','contrast',8,'slices');
%
% Add different colors to an existing plot
% cl_posneg = Meta_Activation_FWE('results', 1, 'height', 'stringent', 'medium', 'poscon'); 
% cl_negpos = Meta_Activation_FWE('results', 1, 'height', 'stringent', 'medium', 'negcon', 'add', 'colors', { [0 0 1] [0 .5 1] [.3 .3 1] });
%
% by Tor Wager, May 2006
% edited by Matthew Davidson June 2006
% modified: tor, Aug 2006
% tor: Sept, 2009: add cl output to "results" mode (minor)
% tor: march 2010, to add bivalent colors 'add' option

% Programmers' notes:
% 4.7.2013  There was a bug in Meta_Activation_FWE that made it not robust to using 
% non-ascending numerical order of contrasts in DB.Contrast.  
% The function meta_prob_activation was reconstructing the maps in ascending order of 
% contrasts, returning maps in ascending sorted contrast order in 
% MC_Setup.unweighted_study_data.  MC_Setup.Xi, the indicators for task type, are 
% sorted in the order of contrasts entered in the database, which respects the order 
% in DB.pointind and DB.connumbers from Meta_Setup.  Tor changed the sort order in 
% meta_prob_activation to use DB.connumbers.  
% 
% This affects MKDA difference analyses (contrasts) when contrast numbers are not 
% entered in ascending numerical order in your database. It also affects classification 
% with Meta_NBC_from_mkda, but not Meta_SVM_from_mkda.

function varargout = Meta_Activation_FWE(meth, varargin)
    varargout = {};

    switch lower(meth)
        
        case 'all'
            MC_Setup = Meta_Activation_FWE('setup', varargin{1});
            Meta_Activation_FWE('mc', varargin{2})
            Meta_Activation_FWE('results')
            
        case 'setup'
            % ------------------------------------------------
            % Setup and write activation probability map
            % ------------------------------------------------
            DB = meta_check_dbsetup(varargin{1});

            docons = input('Compute contrasts across conditions also? (1/0)');
            if docons
                [X,connames, DB, Xi, Xinms, condf, testfield, contrasts] = Meta_Logistic_Design(DB);
                [MC_Setup, activation_proportions] = meta_prob_activation(DB, Xi, contrasts, connames, Xinms);
            else
                [MC_Setup, activation_proportions] = meta_prob_activation(DB);
                disp('Activation map written: Activation_proportion.img');
            end

            % THE CODE BELOW IS REQUIRED FOR BLOB-SHUFFLING ONLY
            % ------------------------------------------------
            % set up blobs for all studies
            % ------------------------------------------------
            str = sprintf('Getting blobs for all studies'); fprintf(1, str);
            s = size(MC_Setup.unweighted_study_data, 2);
            MC_Setup.cl = cell(1,s);
            for i = 1:s
                MC_Setup.cl{i} = iimg_indx2contiguousxyz(MC_Setup.unweighted_study_data(:,i),MC_Setup.volInfo,1);
            end
            erase_string(str);

            % ------------------------------------------------
            % Save setup info
            % ------------------------------------------------
            if exist(fullfile(pwd, 'MC_Info.mat'),'file')
                save MC_Info -append MC_Setup activation_proportions
                disp('Appended info to existing MC_Info.mat.');
            else
                save MC_Info MC_Setup activation_proportions
                disp('Created new MC_Info.');
            end
            
            if exist(fullfile(pwd, 'SETUP.mat'),'file')
                save SETUP -append DB
                disp('Appended DB to existing SETUP.mat.');
            else
                save SETUP DB
                disp('Created new SETUP.mat file and saved DB in it.');
            end
            
            varargout{1} = MC_Setup;

        case 'mc'
            % ------------------------------------------------
            % Monte Carlo iterations (create or add to existing)
            % ------------------------------------------------
            tic
            iter = varargin{1};     % iterations to create or add

            [MC_Setup, maxprop, uncor_prop, maxcsize, startat] = setup_mc_vars(iter);

            last = length(maxprop.act);

            fprintf(1,'Iteration ');

            for i = startat:last

                fprintf(1,'%05d ',i);

                %[maxprop(i),uncor_prop(i,:),maxcsize(i,:)] = meta_stochastic_activation(MC_Setup);

                [mp,up,cs] = meta_stochastic_activation_blobs(MC_Setup);

                maxprop.act(i) = mp.act;
                uncor_prop.act(i,:) = up.act;
                maxcsize.act(i,:) = cs.act;
                if isfield(mp,'poscon')
                    maxprop.poscon(i,:) = mp.poscon;
                    maxprop.negcon(i,:) = mp.negcon;

                    nconds = size(MC_Setup.ctxtxi,1);
                    for j = 1:nconds

                        uncor_prop.poscon{j}(i,:) = up.poscon(j,:);
                        maxcsize.poscon{j}(i,:) = cs.poscon(j,:);

                        uncor_prop.negcon{j}(i,:) = up.negcon(j,:);
                        maxcsize.negcon{j}(i,:) = cs.negcon(j,:);
                    end
                end

                if mod(i,10) == 0 || i == last
                    str = sprintf('Saving results in MC_Info'); fprintf(1,str);
                    save MC_Info -append maxprop uncor_prop maxcsize
                    erase_string(str);
                end

                fprintf(1,'\b\b\b\b\b\b');

            end

            fprintf(1,'Done!\n');
            toc

        case 'results'
            % ------------------------------------------------
            % Get results, threshold, and display
            % ------------------------------------------------
            maptype = 'act';
            conidx = [];

            doplots = 1;
            useheight = 1;
            usestringent = 1;
            usemedium = 1;
            uselenient = 1;
            colors = [1 1 0;1 .5 0;1 .3 .3;.8 .4 .4];   % for activation maps [default]
            dosurface = 0;
            doslices = 0;
            using_defaults = 1;
            dotables = 1;
            addstr = 'noadd';

            read_results_inputs(); % sets up defaults for surfaces, slices, etc.

            % load Setup and prepare thresholds and data vectorized image
            [MC_Setup, statmap, imgprefix, maxprop, uncor_prop, maxcsize, conidx] = setup_results_vars(maptype, conidx);

            % get and plot thresholds
            [maxthr,uncor_thr,maxcthr] = get_thresholds(maxprop, uncor_prop, maxcsize, maptype, doplots, conidx);

            % thresholds to save, based on input
            whomsave = logical([usestringent usemedium uselenient]);

            
            % plot proportions or contrast in statmap
            if doplots
                subplot(1,3,3);
                plot_act_prop(statmap,maxthr,uncor_thr,whomsave);
            end

            % get rid of thresholds we're not interested in seeing in
            % figures
            uncor_thr = uncor_thr(whomsave);
            maxcthr = maxcthr(whomsave);
            
            % Height threshold -- whole brain FWE corrected
            if(useheight)
                i1 = statmap >= maxthr;
                wh = find(i1);
                
                fprintf(1, '\np <= .05 Corrected: Thresh = %3.4f, %3.0f voxels\n\n', maxthr, length(wh));

                % Write results image
                indic2mask(i1, statmap, MC_Setup, [imgprefix '_FWE_height.img']);
            else
                i1 = 0;
            end

            % get clusters to save; old plotting functions
            %cl = mask2clusters([imgprefix '_FWE_height.img']);
            %if doplots, cluster_orthviews(cl); end

            if(usestringent || usemedium || uselenient)
                e1 = 0;
                e2 = 0;
                e3 = 0;

                % Extent threshold --  corrected based on cluster size
                if(usestringent)
                    e1 = meta_cluster_extent_threshold(MC_Setup.volInfo.xyzlist', statmap, uncor_thr(1), maxcthr(1));
                end
                if(usemedium)
                    e2 = meta_cluster_extent_threshold(MC_Setup.volInfo.xyzlist', statmap, uncor_thr(2), maxcthr(2));
                end
                if(uselenient)
                    e3 = meta_cluster_extent_threshold(MC_Setup.volInfo.xyzlist', statmap, uncor_thr(3), maxcthr(3));
                end

                e = e1 | e2 | e3;
                % Write results image
                indic2mask(e, statmap, MC_Setup, [imgprefix '_FWE_extent.img']);
            else
                e = 0;
            end
            % any of the above
            indic = i1 | e;

            % Write results image
            indic2mask(indic, statmap, MC_Setup, [imgprefix '_FWE_all.img']);

            % plot if asked for
            %cl = mask2clusters('Activation_FWE_all.img');
            if doplots

                % multi-threshold view
                [cl,dat] = iimg_multi_threshold(statmap, ...
                    'thresh', [maxthr sort(uncor_thr,'descend')], ...
                    'size', [1 sort(maxcthr)], ...
                    'volInfo', MC_Setup.volInfo,'colors',colors, addstr);

                disp(['Saving clusters as cl variable in ' imgprefix '_clusters.mat']);
                eval(['save ' imgprefix '_clusters cl dat']);

                %cluster_orthviews(cl);
            end

            if exist('cl', 'var'), varargout{1} = cl; end
            
            if doslices
                
                %overlay = which('scalped_single_subj_T1.img');
                overlay = which('spm2_single_subj_T1_scalped.img');
                
                if ~isempty(cl{1})
                    cluster_orthviews_showcenters(cl{1},'coronal',overlay,0,1);
                    saveas(gcf,[imgprefix '_act_coronal_slices'],'fig');
                    saveas(gcf,[imgprefix '_act_coronal_slices'],'png');
                    cluster_orthviews_showcenters(cl{1},'sagittal',overlay,0,1);
                    saveas(gcf,[imgprefix '_act_sag_slices'],'fig');
                    saveas(gcf,[imgprefix '_act_sag_slices'],'png');
                    cluster_orthviews_showcenters(cl{1},'axial',overlay,0,1);
                    saveas(gcf,[imgprefix '_act_axial_slices'],'fig');
                    saveas(gcf,[imgprefix '_act_axial_slices'],'png');
                end

                if length(cl) > 1 && ~isempty(cl{2})
                    cluster_orthviews_showcenters(cl{2},'coronal',overlay,0,1);
                    saveas(gcf,[imgprefix '_extent1_coronal_slices'],'fig');
                    saveas(gcf,[imgprefix '_extent1_coronal_slices'],'png');
                    cluster_orthviews_showcenters(cl{2},'sagittal',overlay,0,1);
                    saveas(gcf,[imgprefix '_extent1_sag_slices'],'fig');
                    saveas(gcf,[imgprefix '_extent1_sag_slices'],'png');
                    cluster_orthviews_showcenters(cl{2},'axial',overlay,0,1);
                    saveas(gcf,[imgprefix '_extent1_axial_slices'],'fig');
                    saveas(gcf,[imgprefix '_extent1_axial_slices'],'png');
                end
                
            end
            
            if dosurface

% %                 % lateral surfaces in orthviews
% %                 create_figure('Lateral surface');
% %                 surfhan = make_surface([],cl{1},cl(2:end),3,colors);
% %                 %
% % 
% %                 saveas(gcf,[imgprefix '_surf'],'png');
% % 
% %                 [hh1,hh2,hh3,hl,a1,a2,a3] = make_figure_into_orthviews;
% %                 view(180,-90); [az,el]=view; h = lightangle(az,el);
% %                 scn_export_papersetup(500)
% % 
% %                 saveas(gcf,[imgprefix '_orth_surf'],'fig');
% %                 saveas(gcf,[imgprefix '_orth_surf'],'png');

                
% %                 % left and right medial
% %                 create_figure('Medial surfaces', 1, 2);
% %                 surfhan = make_surface('hires left',cl{1},cl(2:end),3,colors);
% %                 subplot(1,2,2);
% %                 surfhan = make_surface('hires right',cl{1},cl(2:end),3,colors);
% %                 
                % limbic surfaces
                create_figure('Limbic surface');
                surfhan = make_surface('brainstem',cl{1},cl(2:end),3,colors);
                surfhan = make_surface('cerebellum',cl{1},cl(2:end),3,colors);
                surfhan = make_surface('amygdala',cl{1},cl(2:end),3,colors);
                surfhan = make_surface('nucleus accumbens',cl{1},cl(2:end),3,colors);
                surfhan = make_surface('hippocampus',cl{1},cl(2:end),3,colors);
                surfhan = make_surface('left',cl{1},cl(2:end),3,colors);
                scn_export_papersetup(500)

                saveas(gcf,[imgprefix '_limbic_surf'],'fig');
                saveas(gcf,[imgprefix '_limbic_surf'],'png');
                
                % Note: to replace colored backgrounds with gray, try
                % this:
                %f = findobj(gcf,'Type', 'patch'); set(f, 'FaceColor', 'interp')
                %for i = 1:length(f), set(f(i), 'FaceVertexCData', repmat([.5 .5 .5], size(get(f(i), 'FaceVertexCData'), 1), 1)); cluster_surf(f(i),cl{1},10,{[1 1 0]}); end
                

            end
            
            if dotables
                cluster_table_successive_threshold(cl,10);
            end
            
            
        otherwise
            error('Unknown method.');
    end % switch Method


    % --------------------------------------------------------------------
    %
    % IN-line functions: These have access to all vars in main function
    % Needs matlab 2006a or later
    %
    % --------------------------------------------------------------------

    function read_results_inputs()
        if(length(varargin) > 0), doplots = varargin{1}; end
        if(length(varargin) > 1)
            args = varargin(1:end);  % edited 9/1/2009; do not force first input to be doplots only
            for i=1:length(args)
                if ischar(args{i})  % fixed bug, 9/1/2009
                    switch(args{i})
                        case {'negcon', 'poscon'}
                            maptype = args{i};

                        case {'height', 'stringent', 'medium', 'lenient'}
                            if(using_defaults)
                                useheight = 0;
                                usestringent = 0;
                                usemedium = 0;
                                uselenient = 0;
                                using_defaults = 0;
                            end

                            switch(args{i})
                                case 'height'
                                    useheight = 1;
                                case 'stringent'
                                    usestringent = 1;
                                case 'medium'
                                    usemedium = 1;
                                case 'lenient'
                                    uselenient = 1;
                            end
                        case 'contrast'
                            conidx = args{i+1};

                        case 'surface'
                            dosurface = 1;

                        case 'colors'
                            colors = args{i+1};

                        case 'grayscale'
                            colors = [0 0 0; .25 .25 .25; .4 .4 .4; .7 .7 .7];

                        case 'slices'
                            doslices = 1;
                            
                        case 'add'
                            addstr = 'add';
                            
                    end
                end
            end
        end
    end
end



% ====================================================================
% --------------------------------------------------------------------
%
% Monte Carlo subfunctions
%
% --------------------------------------------------------------------
% ====================================================================

function [MC_Setup, maxprop, uncor_prop, maxcsize,startat] = setup_mc_vars(iter)

    % These should be loaded in file
    MC_Setup = []; maxprop = [];  uncor_prop = []; maxcsize = [];

    if ~exist('MC_Info.mat','file')
        error('You must have an MC_Info file in the current dir that contains MC_Setup.');
    else
        load MC_Info
    end

    if isempty(maxprop) % no existing results; initialize
        fprintf(1,'No MC iterations found yet.  Creating maxprop, uncor_prop, maxcsize\n');

        maxprop.act = zeros(iter,1);
        uncor_prop.act = zeros(iter,3);
        maxcsize.act = zeros(iter,3);

        if isfield(MC_Setup,'ctxtxi')  % do by-conditions also
            nconds = size(MC_Setup.ctxtxi,1);

            maxprop.poscon = zeros(iter, nconds);
            maxprop.negcon = zeros(iter, nconds);

            for i = 1:nconds
                uncor_prop.poscon{i} = zeros(iter,3);
                maxcsize.poscon{i} = zeros(iter,3);

                uncor_prop.negcon{i} = zeros(iter,3);
                maxcsize.negcon{i} = zeros(iter,3);
            end
        end
    else
        % existing iterations; append
        % check for empty (incomplete) iterations and remove them
        [maxprop, uncor_prop, maxcsize] = remove_incomplete_iterations(maxprop, uncor_prop, maxcsize);

        fprintf(1,'MC iterations found: %3.0f.  Appending new: %3.0f iterations\n', length(maxprop.act), iter);

        maxprop.act = [maxprop.act; zeros(iter,1)];
        uncor_prop.act = [uncor_prop.act; zeros(iter,3)];
        maxcsize.act = [maxcsize.act; zeros(iter,3)];

        if isfield(MC_Setup,'ctxtxi')  % do by-conditions also

            nconds = size(MC_Setup.ctxtxi,1);

            maxprop.poscon = [maxprop.poscon; zeros(iter, nconds)];
            maxprop.negcon = [maxprop.negcon; zeros(iter, nconds)];

            for i = 1:nconds
                uncor_prop.poscon{i} = [uncor_prop.poscon{i}; zeros(iter,3)];
                maxcsize.poscon{i} = [maxcsize.poscon{i}; zeros(iter,3)];

                uncor_prop.negcon{i} = [uncor_prop.negcon{i}; zeros(iter,3)];
                maxcsize.negcon{i} = [maxcsize.negcon{i}; zeros(iter,3)];

            end
        end
    end

    startat = find(maxprop.act == 0); startat = startat(1);
end





% ====================================================================
% --------------------------------------------------------------------
%
% Results subfunctions
%
% --------------------------------------------------------------------
% ====================================================================


% SUBFUNCTION:
% --------------------------------------------------------------------
% Load some key variables for getting results
% --------------------------------------------------------------------
function [MC_Setup, statmap, imgprefix, maxprop, uncor_prop, maxcsize, conidx] = setup_results_vars(maptype, conidx)

    % should be loaded in file
    MC_Setup = []; maxprop = [];  uncor_prop = []; maxcsize = [];

    % load file
    % ------------
    if ~exist('MC_Info.mat','file')
        error('You must have an MC_Info file in the current dir that contains MC_Setup.');
    else
        load MC_Info
    end

    % check mc iterations
    % ------------
    if isempty(maxprop)
        % no existing results; error
        fprintf(1,'No MC iterations found in MC_Info. You must create using Meta_Activation_FWE(''mc'',num_iterations)\n');
    end

    % check for empty (incomplete) iterations and remove them
    [maxprop, uncor_prop, maxcsize] = remove_incomplete_iterations(maxprop, uncor_prop, maxcsize);


    % Set up statistic image data to be thresholded
    % ------------
    switch maptype
        case 'act'
            % overall activation
            % use existing props, or set up empty to read from file
            if ~(exist('activation_proportions', 'var')), activation_proportions = []; end
            % get activation proportions indicator or just check it
            statmap = setup_activation_proportions(activation_proportions, MC_Setup);

            imgprefix = 'Activation';
        otherwise
            % contrast across conditions
            % statmap = MC_Setup.con_data;

            if(~exist('conidx', 'var') || isempty(conidx))
                num_contrasts = size(MC_Setup.con_data, 2);
                if num_contrasts > 1
                    conidx = input(['Choose contrast number (1 thru ' num2str(num_contrasts) '): ']);
                else
                    conidx = 1;
                end
            end
            statmap = MC_Setup.con_data(:,conidx);

            % check for correct field names
            switch maptype
                case 'poscon'
                    imgprefix = [MC_Setup.connames{conidx} '_Pos'];
                case 'negcon'
                    % flip sign so thresholding will work the same way
                    statmap = -statmap;

                    imgprefix = [MC_Setup.connames{conidx} '_Neg'];
                otherwise
                    error('Unknown maptype: %s - must be ''act'', ''poscon'', or ''negcon''', maptype);
            end
    end
end

% SUBFUNCTION:
% --------------------------------------------------------------------
% Check that we have activations loaded, or load them from image
% --------------------------------------------------------------------

function activation_proportions = setup_activation_proportions(activation_proportions, MC_Setup)
    % check/load true meta activations
    % ------------
    if isempty(activation_proportions)
        name = 'Activation_proportion.img';
        str = sprintf('Loading %s', name);fprintf(1,str);
        if ~exist(name,'file')
            error('You must have an Activation_proportion.img file in the current dir, or run ''all'' option to create.');
        end

        mask = MC_Setup.volInfo.fname;
        if ~exist(mask,'file')
            error(['Cannot find mask image: ' mask]);
        end

        %         switch maptype
        %             case 'poscon'
        %                 imgprefix = [MC_Setup.connames(conidx) '_Pos'];
        %             case 'negcon'
        %                 % flip sign so thresholding will work the same way
        %                 statmap = -statmap; %*****
        %
        %                 imgprefix = [MC_Setup.connames(conidx) '_Neg'];
        %             otherwise
        %                 error('map type (3rd arg.) must be ''act'', ''poscon'', or ''negcon''');
        %         end
        activation_proportions = meta_read_image(name,mask);
        erase_string(str);
    end

    % check length
    len = size(MC_Setup.volInfo.xyzlist, 1);
    if len ~= length(activation_proportions)
        error('activation image and results have different numbers of in-analysis voxels.  Has the mask image changed??');
    end
end

% SUBFUNCTION:
% --------------------------------------------------------------------
% Get significance thresholds and plot histograms
% --------------------------------------------------------------------

function [maxthr,uncor_thr,maxcthr] = get_thresholds(maxprop, uncor_prop, maxcsize, maptype, doplots, conidx)

    % flip signs of thresholds for negative contrast so code will work generically
    % cell only if maptype is poscon or negcon (not act)
    switch(maptype)
        case 'act'
            current_maxprop = maxprop.(maptype);
            current_uncor_prop = uncor_prop.(maptype);
            current_maxcsize = maxcsize.(maptype);
        case 'poscon'
            current_maxprop = maxprop.(maptype)(:,conidx);
            current_uncor_prop = uncor_prop.(maptype){conidx};
            current_maxcsize = maxcsize.(maptype){conidx};
        case 'negcon'
            current_maxprop = -maxprop.(maptype)(:,conidx);
            current_uncor_prop = -uncor_prop.(maptype){conidx};
            current_maxcsize = maxcsize.(maptype){conidx};
    end

    maxthr = prctile(current_maxprop, 95);
    uncor_thr = mean(current_uncor_prop);
    maxcthr = prctile(current_maxcsize, 95);

    if doplots
        tor_fig(1,3);

        [h,x] = prepare_histogram(current_maxprop);

        plot(x,h,'k','LineWidth',2);
        title('Null hypothesis activation proportion');
        yl = get(gca,'YLim');
        plot([maxthr maxthr],yl,'Color',[.5 .5 .5],'LineWidth',2);
        sym = {'--' '.-' ':'};

        for i = 1:3
            plot([uncor_thr(i) uncor_thr(i)],yl,sym{i},'Color',[.5 .5 .5],'LineWidth',2);
        end
        legend({'Max whole-brain density' 'p < .05 corrected' 'p < .05' 'p < .01' 'p < .001'});
        drawnow

        subplot(1,3,2);
        [h,x] = prepare_histogram(current_maxcsize);

        plot(x,h,'LineWidth',2);
        colors = {'b' 'g' 'r'};
        yl = get(gca,'YLim');
        title('Null hypothesis cluster sizes');

        for i = 1:3
            plot([maxcthr(i) maxcthr(i)],yl,sym{i},'Color',colors{i},'LineWidth',2);
        end
        legend({'p < .05' 'p < .01' 'p < .001'});
        drawnow
    end
end

% --------------------------------------------------------------------
% Plot activation proportions with thresholds
% --------------------------------------------------------------------
function plot_act_prop(activation_proportions,maxthr,uncor_thr,whomsave)

    if nargin < 4, whomsave = [1 1 1]; end  % all thresholds
    [h,x] = prepare_histogram(activation_proportions);

    plot(x,h,'k','LineWidth',2);
    title('Observed activation proportions');
    yl = get(gca,'YLim');
    plot([maxthr maxthr],yl,'Color',[.5 .5 .5],'LineWidth',2);
    sym = {'--' '.-' ':'};

    for i = 1:3
        if whomsave(i)
            plot([uncor_thr(i) uncor_thr(i)],yl,sym{i},'Color',[.5 .5 .5],'LineWidth',2);
        end
    end
    
    legstr = {'Observed proportions' 'p < .05 corrected' 'p < .05' 'p < .01' 'p < .001'};
    legstr = legstr(find([1 1 whomsave]));
    legend(legstr);
    drawnow

end


% SUBFUNCTION:
% --------------------------------------------------------------------
% Get significant regions based on extent
% --------------------------------------------------------------------

function  indic = meta_cluster_extent_threshold(xyzlist,data,thr1,sizethr)


    norig = size(xyzlist,2);
    indic = zeros(norig,1);

    meet_primary = data > thr1;
    wh1 = find(meet_primary);        % voxels that do meet primary


    % consider only vox that meet primary
    wh = find(~meet_primary);        % voxels that do not meet primary threshold
    xyzlist(:,wh) = [];
    clear wh

    % check size OK
    n = size(xyzlist,2);

    if n > 50000
        disp('Too many voxels meet primary threshold. Cannot perform cluster-based thresholding.');
        return
    elseif n == 0
        disp('Warning!  No voxels meet primary threshold; this is unusual...');
        return
    else

        % cluster
        clid = spm_clusters(xyzlist);
    end

    csize = zeros(norig,1);     % in Original (full voxel) index

    u = unique(clid);
    for i = 1:length(u)
        wh = find(clid == i);   % voxels in this cluster
        csize(wh1(wh)) = length(wh); % put cluster size there; translate to full index
    end

    indic = csize > sizethr;

    % results display stuff
    sigvox = find(indic);

    fprintf(1,'\nThresh = %3.4f, size threshold = %3.0f: %3.0f sig. voxels.\n',thr1,sizethr,length(sigvox));

    if isempty(sigvox)
        return
    else
        rg = [min(csize(sigvox)) max(csize(sigvox))];
    end
    fprintf(1,'Cluster sizes range from %3.0f to %3.0f\n\n', rg(1), rg(2));


end

% SUBFUNCTION:
% --------------------------------------------------------------------
% Quick way to write a mask image for significant data
% --------------------------------------------------------------------

function indic2mask(indic, activation_proportions, MC_Setup, name)
    disp(['Writing image: ' name])
    wh = find(indic);

    data = zeros(size(activation_proportions));
    data(wh) = activation_proportions(wh);
    meta_reconstruct_mask(data, MC_Setup.volInfo.xyzlist, MC_Setup.volInfo.dim(1:3), 1, MC_Setup.volInfo,name);
end


% SUBFUNCTION:
% --------------------------------------------------------------------
% Surface plots at multiple thresholds
% --------------------------------------------------------------------

function surfhan = make_surface(typestr,cl,morecl,mydistance,color)
    % NOTE: this is the same as in cluster_surf_batch, except for input
    % type for color

    if isempty(typestr), typestr = mydistance; end  % kludgy fix to avoid error on empty input

    surfhan = [];   %cluster_surf(cl,mydistance,color{1});
    for i = length(morecl):-1:1
        if i == length(morecl)
            surfhan = cluster_surf(surfhan,morecl{i},mydistance,{color(i+1,:)},typestr);
        else
            cluster_surf(surfhan,morecl{i},mydistance,{color(i+1,:)});
        end
        %surfhan = [surfhan sh];
    end

    surfhan = cluster_surf(surfhan,cl,mydistance,{color(1,:)}, typestr);
    %surfhan = [surfhan sh];

end


% ====================================================================
% --------------------------------------------------------------------
%
% Other utility functions
%
% --------------------------------------------------------------------
% ====================================================================

function erase_string(str)
    fprintf(1,repmat('\b',1,length(str))); % erase string
    %fprintf(1,'\n');
end


function h = makeColumn(h)
    if length(h) > size(h,1) && size(h,1)==1, h = h'; end
end


% SUBFUNCTION:
function [h,x,nbins] = prepare_histogram(data)

    nbins = min(length(data)./5,100);

    [h,x] = hist(data,nbins);

    % pad to make nice display
    if length(x)>1, firstx = 2*x(1) - x(2); else firstx = 0; end
    lastx = 2*x(end) - x(end-1);

    % make column
    h = makeColumn(h);
    x = makeColumn(x);

    z = zeros(1,size(h,2));

    h = [z;h;z]; x = [firstx;x;lastx];

    h = h ./ repmat(sum(h),size(h,1),1);

end



% SUBFUNCTION:
function DB = meta_check_dbsetup(DB)
    % runs Meta_Setup if necessary

    % needs to set up design:
    %DB.(fields)   % lists of fields containing task conditions for each coordinate point
    %DB.pointind   % indices of which coord points are in which unique
    %                contrast
    %DB.x          % x coordinates, just to get # of points
    % DB.maskname
    % DB.radius
    % DB.Contrast      (study ID number; contrast within study)
    % DB.x, DB.y DB.z
    % DB.studyweight

    if isempty(DB), error('No DB structure! Create by reading from text file using read_database.m or create manually.'); end

    DBfields = fieldnames(DB);
    reqfields = {'maskname' 'radius' 'Contrast' 'x' 'y' 'z' 'studyweight' 'pointind'};

    for i = 1:length(reqfields)
        field_exists(i) = any(~(cellfun(@isempty, strfind(DBfields, 'maskname'))));
        if ~field_exists(i), fprintf(1,'DB.%s does not exist yet.\n', reqfields{i}); end
    end

    if ~all(field_exists)
        fprintf(1,'Running Meta_Setup to set up DB structure.\n');
        DB = Meta_Setup(DB);
    end

end


% SUBFUNCTION:
function [maxprop, uncor_prop, maxcsize] = remove_incomplete_iterations(maxprop, uncor_prop, maxcsize)
    %
    % check for empty (incomplete) iterations and remove them
    if any(maxprop.act == 0)
        wh = find(sum(maxprop.act,2) == 0); % empty iterations
        fprintf(1,'Warning! %3.0f Empty iteration results: Removing from all iteration variables.\n',length(wh));
        maxprop.act(wh) = [];
        uncor_prop.act(wh,:) = [];
        maxcsize.act(wh,:) = [];

        if isfield(uncor_prop,'poscon')  % do by-conditions also

            nconds = size(maxprop.poscon,2);

            maxprop.poscon(wh,:) = [];
            maxprop.negcon(wh,:) = [];

            for i = 1:nconds
                wh = find(sum(uncor_prop.poscon{i},2) == 0); % empty iterations
                % do separately because of potential for mismatch between
                % poscon and act
                uncor_prop.poscon{i}(wh,:) = [];
                maxcsize.poscon{i}(wh,:) = [];

                wh = find(sum(uncor_prop.negcon{i},2) == 0); % empty iterations
                uncor_prop.negcon{i}(wh,:) = [];
                maxcsize.negcon{i}(wh,:) = [];
            end
        end
    end

end
