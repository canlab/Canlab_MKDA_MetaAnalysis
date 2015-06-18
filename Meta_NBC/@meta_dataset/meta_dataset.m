classdef meta_dataset
    % Define a meta-analysis dataset object for classification analyses
    %
    % Properties:
    % ---------------------------------------------------------------
    %     dat               Data; sparse logical
    %     connames          cell array of contrast names
    %     volInfo           volInfo structure for 3-D reconstruction; see iimg_get_data.m
    %     classes           Contrast classes; integers
    %     classnames        Names for classes; cell of strings
    %     params            Parameter structure
    %     terms             Terms to classify
    %     include_cases     Logical vector of contrasts to include
    %     include_vox       Logical vector of voxels to include (features)
    %     train_pred_fcn    Function handle to create predicted classes from binary data
    %     test__pred_fcn    Function handle to create predicted classes
    %     loss_fcn          Function handle of loss (error) function to evaluate
    %   
    %
    % Class constructor: Creating a meta_dataset object
    % ---------------------------------------------------------------
    % obj = meta_dataset(varargin)
    % varargin can be any property, value pair
    %
    % Examples:
    % data = meta_dataset;
    % data = meta_dataset('dat', MC_Setup.unweighted_study_data);
    %
    % data = meta_dataset('MC_Setup', MC_Setup); % special to convert MC_Setup structure
    %
    % Full classification example for pain vs. touch (from Neurosynth):
    % data = meta_dataset('MC_Setup', MC_Setup);
    % data = setup_data(data, {'pain', 'touch'});
    % nbc = train(meta_nbc, data);
    % nbc = test(nbc, data.dat);
    % nbc = get_error(nbc, data.classes);
    % report(nbc);
    %
    % Full classification example 2:
    %   ---------------------------------------------------------
    % DB = Meta_Setup(DB, 10);
    % MC_Setup = Meta_Activation_FWE('setup', DB); % Here, choose Valence / positive negative
    %            - you must enter a contrast (any contrast) across the
    %            levels you wish to classify
    % data = meta_dataset('MC_Setup', MC_Setup);
    % data.params.min_vox = 50; % do not exclude sparse maps, just those < 50 vox
    % data = setup_data(data, {'positive', 'negative'});
    % % Cross-validate and get cv error rate
    % nbc = cv(nbc, data);
    % nbc = report(nbc);

    
    properties
        
        dat
        connames
        volInfo
        classes
        classnames
        params
        terms
        include_cases
        include_vox
        train_pred_fcn
        test_pred_fcn
        loss_fcn
        
        
    end
    
    methods
        
        % Class constructor
        function obj = meta_dataset(varargin)
            %
            % obj = meta_dataset(varargin)
            
            % ---------------------------------
            % Create empty fmri_data object, and return if no additional
            % arguments
            % ---------------------------------
            
            obj.dat = sparse(logical(false));
            obj.connames = {};
            obj.volInfo = [];
            obj.classes  = single(0);
            obj.classnames = {};
            obj.terms = {};
            obj.include_cases = logical(false);
            obj.include_vox = logical(false);
            
            obj.train_pred_fcn = [];
            obj.test_pred_fcn = [];
            obj.loss_fcn = [];
            
            obj.params.min_vox = 5000;  % Studies must report at least this many sig voxels (after kernel conv)
            obj.params.term_freq_cutoff = 0.001; % Minimum frequency of word use for study to be counted as involving term
            obj.params.min_prop_active = 0.03;  % Voxels must have at least this proportion of active studies
            
            % Varying the threshold has a substantial effect. A
            % cut-off around .001 seems to work well for different
            % datasets.
            obj.params.predict_pval_thresh = 0.05;
            
            obj.params.k = 10;  % number of cross-validation folds
            obj.params.s_param = 2; % smoothing param for Naive Bayes m-estimator; Caffo says 2
            
            % The code below can be generic to any class definition
            % It parses 'fieldname', value pairs of inputs
            % and returns a warning if unexpected strings are found.
            
            if nargin == 0
                return
            end
            
            % all valid fieldnames
            valid_names = fieldnames(obj);
            valid_names = [valid_names; {'MC_Setup', 'MC_setup', 'mc_setup'}'];
            
            for i = 1:length(varargin)
                if ischar(varargin{i})
                    
                    % Look for a field (attribute) with the input name
                    wh = strmatch(varargin{i}, valid_names, 'exact');
                    
                    % behaviors for valid fields
                    if ~isempty(wh)

                        % special methods for specific fields
                        switch varargin{i}
                            
                            % Re-cast vars as correct data types
                            case {'include_cases' 'include_vox'}
                                obj.(varargin{i}) = logical(varargin{i+1});
                                
                            case 'dat'
                                obj.(varargin{i}) = sparse(logical(varargin{i+1}));
                            
                            case {'MC_Setup', 'MC_setup', 'mc_setup'}
                                M = varargin{i+1};
                                M.unweighted_study_data = sparse(logical(M.unweighted_study_data));
                                obj.dat = M.unweighted_study_data;
                                M = rmfield(M, 'unweighted_study_data'); % save memory space
                                obj.volInfo = M.volInfo;
                                
                                %This worked for Neurosynth...below is mod
                                %for MC_Setup from Meta_Activ_FWE
                                %obj.connames = M.connames;
                                
                                obj.connames = cell(size(obj.dat, 2), 1);
                                if isfield(M, 'Xi')
                                    for j = 1:size(M.Xi, 2)
                                        obj.connames(logical(M.Xi(:, j))) = M.Xinms(j);
                                    end
                                end
                                
                            otherwise
                                obj.(varargin{i}) = varargin{i + 1};
                        end
                        
                        % eliminate strings to prevent warnings on char
                        % inputs
                        if ischar(varargin{i + 1})
                            varargin{i + 1} = [];
                        end

                    else
                        warning('inputargs:BadInput', 'Unknown field: %s', varargin{i});
                        
                    end
                end % string input
            end % process inputs
            
        end % constructor function
        
    end % methods
    
end % classdef

