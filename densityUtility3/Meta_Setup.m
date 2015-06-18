% DB = Meta_Setup(DB, [radius_mm], [con_dens_images])
% Set up Meta-analysis dataset
%
% use after read_database.m
% See the Manual for more complete information.
%
% Special Fields
%
% Subjects or N     : sample size
% FixedRandom       : fixed or random effects
% Subjective Weights : weighting vector based on FixedRandom and whatever
%   else you want to weight by; e.g., study reporting threshold
% x, y, z           : coordinates
% study or Study    : name of study
% Contrast          : unique indices (e.g., 1:k) for each independent
%                     contrast

function DB = Meta_Setup(DB, varargin)

    diary ANALYSIS_INFORMATION.txt

    write_contrast_image_flag = 0;  % flag to write contrast images;
    % Contrast images are not needed for Meta_Activation_FWE

    %P = ['brain_avg152T1.img'];     % 2 mm voxels
    P = 'scalped_avg152T1_graymatter_smoothed.img';

    if length(varargin) > 0
        radius_mm = varargin{1};
    else
        radius_mm = 10;
    end
    disp(['Radius is ' num2str(radius_mm) ' mm.']);

    %mask_file = P;
    t1 = clock;
    fprintf(1,'Setup. ')

    % -----------------------------------------------------
    % * load standard brain
    % -----------------------------------------------------

    switch spm('Ver')
        case 'SPM2'
            % spm_defaults is a script
            disp('WARNING: spm defaults not set for spm2. Make sure your defaults are set correctly');

        case 'SPM5'
            % spm_defaults is a function
            spm_defaults()

        case 'SPM8'
            % spm_defaults is a function
            spm_defaults()

        otherwise
            % unknown SPM
            disp('Unknown version of SPM!');
            spm_defaults()
    end
    
    defaults.analyze.flip = 0;

    P = which(P);
    if isempty(P), disp('Error: Cannot find mask image on path: scalped_avg152T1_graymatter_smoothed.img. Nothing done.'); return; end

    V = spm_vol(P);
    mask = zeros(V.dim(1:3));

    voxsize = diag(V(1).mat)';
    voxsize = voxsize(1:3);
    radius = radius_mm ./ mean(abs(voxsize));
    sphere_vol = 4 * pi * radius_mm ^ 3 / 3;

    DB.radius_mm = radius_mm;
    DB.radius = radius;

    % mask image -- gray matter
    DB.maskname = P;
    DB.maskV = V;
    DB.sphere_vol = sphere_vol;
    DB.voxsize = voxsize;
    DB.maskzeros = mask;


    % -----------------------------------------------------
    % identify field info to save with each contrast image in DB struct
    % later saved in SETUP.mat
    % -----------------------------------------------------
    N = fieldnames(DB); fieldn = {};
    for i = 1:length(N)
        str = ['tmp = DB.' N{i} ';'];
        eval(str)
        if length(tmp) == length(DB.x)  % if this field has all info for each peak
            fieldn(end+1) = N(i);
        end
    end




    % -----------------------------------------------------
    % Make sure we have some things that are often mis-labeled
    % -----------------------------------------------------

    % Special Fields
    %
    % Subjects or N     : sample size
    % FixedRandom       : fixed or random effects
    % Subjective Weights : weighting vector based on FixedRandom and whatever
    %   else you want to weight by; e.g., study reporting threshold
    % x, y, z           : coordinates
    % study or Study    : name of study
    % Contrast          : unique indices (e.g., 1:k) for each independent
    %                     contrast

    if ~isfield(DB,'Subjects') && isfield(DB,'N')
        disp('Cannot find Subjects field.  Copying from N.')
        DB.Subjects = DB.N;
    end


    if ~isfield(DB,'study') && isfield(DB,'Study')
        disp('Cannot find study field.  Copying from Study.')
        DB.study = DB.Study;
    end

    if ~isfield(DB,'x') && isfield(DB,'X')
        disp('Cannot find x field.  Copying from X.')
        DB.x = DB.X;
    end

    if ~isfield(DB,'y') && isfield(DB,'Y')
        disp('Cannot find y field.  Copying from Y.')
        DB.y = DB.Y;
    end

    if ~isfield(DB,'z') && isfield(DB,'Z')
        disp('Cannot find z field.  Copying from Z.')
        DB.z = DB.Z;
    end

    if ~isfield(DB,'xyz'),
        disp('No xyz field in DB.  Creating from x, y, z.');
        DB.xyz = [DB.x DB.y DB.z];
    end

    if ~isfield(DB,'Subjects'), error('Enter Subjects or N field.'); end

    % -----------------------------------------------------
    % identify independent contrasts and get sqrt(N)
    % -----------------------------------------------------


    % kludgy fix if you haven't entered indices of independent contrasts in
    % Contrast.  You should never use this with the new meta-analyses.
    % -----------------------------------------------------
    if ~isfield(DB,'Contrast')
        disp('No DB.Contrast field.  You can create one using unique levels of a variable,')
        disp('or create one in the text file and re-run.')

        disp('Run: get_contrast_indicator_improved')
        disp('Store results: DB.Contrast = contrasts; DB.pointind = first_peak_in_each_con;');

        % [contrasts, first_peak_in_each_con] = get_contrast_indicator_improved(DB, {'Memtype' 'Control' 'Design' 'Feature' 'Material'});
        % DB.Contrast = contrasts;
        % DB.pointind = first_peak_in_each_con';
        %
        %
        %
        %         testfield = input('Enter field name to use: ','s');
        %
        %         OUT = get_contrast_indicator(DB,testfield);
        %
        %         pointind = OUT.conindex';
        %         DB.pointind = pointind;
        %         for i = 1:length(pointind),
        %             ind = pointind; ind(ind <= pointind(i)) = []; nextlargest = min(ind);
        %             DB.Contrast(pointind(i):nextlargest-1) = i;
        %         end
        %
        %         [DB.connumbers,DB.pointind] = unique(DB.Contrast);
        %         DB.rootn = OUT.rootn';
        %
        %         disp('Warning!!! This part probably has bugs!!!')
        %         DB.connumbers(DB.connumbers == 0) = max(DB.connumbers)+1;
        %         DB.Contrast(end+1) = DB.Contrast(end);


    else
        % we have unique contrasts
        % you have entered them in the Contrast variable
        % this is the right thing to do.
        % -----------------------------------------------------

        DB.rootn = sqrt(DB.Subjects);
        [DB.connumbers,DB.pointind] = unique(DB.Contrast);
        DB.rootn = DB.rootn(DB.pointind);

    end

    % -----------------------------------------------------
    % Make sure contrast weights are sorted
    % Needed later for Meta_Select_Contrasts and Log. Design
    % Note: this should not actually affect order of contrasts, which
    % should still be in order of entry in dataset. Other functions should
    % respect this order, however.  See programmers' notes in
    % Meta_Activation_FWE
    % -----------------------------------------------------
    [DB.pointind,ii] = sort(DB.pointind);
    DB.rootn = DB.rootn(ii);
    DB.connumbers = DB.connumbers(ii);

    % -----------------------------------------------------
    % get Final contrast weights
    % -----------------------------------------------------

    if isfield(DB,'SubjectiveWeights')
        disp('Scaling by SubjectiveWeights');
        w = DB.rootn .* DB.SubjectiveWeights(DB.pointind);
    else
        w = DB.rootn;
    end

    % make sure no NaNs or zeros; impute mean
    whbad = find(isnan(w) | w == 0);
    if ~isempty(whbad)
        disp(['Warning! ' num2str(length(whbad)) ' contrasts have bad or missing weights. Imputing mean. ']);
        wok = w; wok(whbad)=[];
        w(whbad) = mean(wok);
    end

    % these must sum to 1 !
    DB.studyweight = w ./ sum(w);


    % -----------------------------------------------------
    % plot studies with weights
    % -----------------------------------------------------
    DB.connames = DB.study(DB.pointind);

    figure('Color', 'w');
    set(gca, 'FontSize', 14);
    plot(DB.rootn, 'k:.', 'LineWidth', 2);
    hold on;
    plot(zscore(DB.studyweight)*2+max(DB.rootn)*1.5, 'k.-','LineWidth', 2);
    set(gca, 'XTick', 1:length(DB.pointind), 'XTickLabel', DB.connames);
    try
        legend({'Sqrt(N)' 'Study weights (scaled for display)'})
    catch
        disp('Failed to make legend...Matlab 2010b bug?');
    end
        
    drawnow


    % save input variables in this directory
    disp('DB structure saved in SETUP.mat');
    save SETUP DB
    diary off


    % -----------------------------------------------------
    % for each contrast, make a density image
    % -----------------------------------------------------
    PP = [];     % image names
    if length(varargin) > 1, PP = varargin{2}; end

    if write_contrast_image_flag
        DB = write_contrast_images(DB, t1, PP);
    else
        disp('Contrast image writing is OFF.  Use Meta_Activation_FWE for analysis.');
    end

    disp('DB struct output of Meta_Setup saved in SETUP.mat');
    disp('Use info in this file to run next step');
    save SETUP DB
end


% ====================================================================
% --------------------------------------------------------------------
%
% subfunctions
%
% --------------------------------------------------------------------
% ====================================================================

function str = remove_special_chars(str)
    wh = (str == '''' | str=='"' | str == '&' | str == '~' | str == '/' | str == '\');
    wh = find(wh);
    str(wh) = [];
end



function DB = write_contrast_images(DB,t1,varargin)

    PP = [];     % image names
    if length(varargin) > 0, PP = varargin{1}; end            % load input filenames here!!


    if isempty(PP)

        % normalize density maps so that 1 activation = value of 1 in map

        % Try to find images first, and if they don't exist in current dir,
        % create.
        fprintf(1,'Creating images.')

        DB.XYZmm = {};

        for i = 1:length(DB.pointind)

            % peaks in this study, this contrast
            wh = find(DB.Contrast == DB.Contrast(DB.pointind(i)));

            studyname = DB.study{wh(1)};

            studyname = remove_special_chars(studyname);

            XYZmm = DB.xyz(wh,:);

            DB.XYZmm{i} = XYZmm;

            str = [studyname '_contrast_' num2str(DB.connumbers(i)) '.img'];    % name of image
            if exist(str,'file'), disp(['Found existing: ' str])

            else
                % doesn't exist yet; create it
                % NOTE: CHANGED 1/3/05 TO NOT WEIGHT BY STUDYWEIGHT HERE.
                % Weights are added later, in Meta_Activation and Meta_Logistic
                % This ensures that if you select a subset of studies with
                % Meta_Select_Contrasts, the weights will be recomputed and
                % reapplied correctly.
                % -------------------------------------------------
                %xyz2density(XYZmm,mask,DB.maskV,str,DB.radius,DB.studyweight(i));
                xyz2density(XYZmm,DB.maskzeros,DB.maskV,str,DB.radius,1);

            end % if find existing

            PP = str2mat(PP,str);

        end % loop thru contrasts
        PP = PP(2:end,:);
    else
        % we have filenames already created
    end

    fprintf(1, 'Done %3.0f images in %3.0f s\n', length(DB.pointind), etime(clock, t1))

    DB.PP = PP;
    DB.PP = [repmat([pwd filesep], size(DB.PP, 1), 1) DB.PP];
end
