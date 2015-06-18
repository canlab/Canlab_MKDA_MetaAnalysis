function obj = train(obj, meta_data_obj)
% nbc_obj = train(nbc_obj, meta_data_obj)
%
% Training function for Naive Bayes classifier
% Tor Wager
% - Modified by Tal Yarkoni 10/27/2010
%
% Parameters used:
% -----------------------------------------------------------------------
% s_param, default = 2; % smoothing param for shrinking near-zero probabilities; Caffo says 2
%
% Inputs:
% -----------------------------------------------------------------------
% obj, a meta_nbc class object
% meta_data_obj, a meta_dataset class object
%
% Outputs:
% -----------------------------------------------------------------------
% meta_nbc class object with assigned values for these fields:
% n_classes
% p_act
% prior_probs
% train_pprobs_act
% train_pprobs_noact
% class_map
% 

% assign all parameters in object to their respective variable names
N = fieldnames(obj.params);

for i = 1:length(N)
    eval([N{i} ' = obj.params.' N{i} ';'])
end

if isempty(obj.trainIdx)
    fprintf('trainIdx is empty. Using all contrasts for training.\n')
    obj.trainIdx = true(size(meta_data_obj.dat, 2), 1);
end

% Select training data from dataset object
% -----------------------------------------------------------------------

% Matlab's crossval function uses nominal vars, which are incompatible with
% some commands below. So recast, and select samples:
train_class = single(meta_data_obj.classes(obj.trainIdx));

obj.n_classes = length(unique(train_class(train_class > 0)));

% was obs x voxels (for cv compat). make vox x obs
train_maps = meta_data_obj.dat(:, obj.trainIdx); 

% -----------------------------------------------------------------------
% Activation probabilities
% -----------------------------------------------------------------------
n_train = size(train_maps, 2);

for ii = 1:obj.n_classes
    
    n_term(ii) = sum(train_class == ii);
    
    % p(activation = yes | class)
    obj.p_a_given_c(:, ii) = (sum(train_maps(:, train_class == ii), 2)+s_param) ./ (n_term(ii)+2*s_param);
   
    % p(activation = no | class)
    obj.p_noa_given_c(:, ii) = (sum(~train_maps(:, train_class == ii), 2)+s_param) ./ (n_term(ii)+2*s_param);
    
    % p(activation = yes | no class)
    obj.p_a_given_no_c(:, ii) = (sum(train_maps(:, train_class ~= ii), 2)+s_param) ./ (n_train-n_term(ii)+2*s_param);
    
    % p(activation = no | no class)
    obj.p_noa_given_no_c(:, ii) = (sum(~train_maps(:, train_class ~= ii), 2)+s_param) ./ (n_train-n_term(ii)+2*s_param);
    
end

% -----------------------------------------------------------------------
% Base rates
% -----------------------------------------------------------------------

% Now calculate
% p(class | activation)
% act_term * p(class) / p(activation)
% balanced error rate: p(class) = 1/obj.n_classes
% need to recalculate p(activation) based on balanced error rate

obj.prior_probs = ones(obj.n_classes, 1) ./ obj.n_classes;
obj.p_act = obj.p_a_given_c * obj.prior_probs; %+ obj.p_a_given_no_c * (1-obj.prior_probs);

%p_noact = p_noa_given_c * obj.prior_probs;
p_noact = 1 - obj.p_act;

% -----------------------------------------------------------------------
% Posterior probabilities of class given activation or lack of activation
% at each voxel
% -----------------------------------------------------------------------

% p(class | activation = yes) ... this is a matrix of vox x classes
obj.train_pprobs_act = obj.p_a_given_c * diag(obj.prior_probs) ./ obj.p_act(:, ones(1, obj.n_classes));

% p(class | activation = no) ... this is a matrix of vox x classes
obj.train_pprobs_noact = obj.p_noa_given_c * diag(obj.prior_probs) ./ p_noact(:, ones(1, obj.n_classes));

[discard, obj.class_map] = max(obj.train_pprobs_act');
obj.class_map = obj.class_map';

% these values are scaled more appropriately than act_term values.
% there are not a large number of zero values, for example, because
% the base rate is taken into consideration.

% save volInfo and excluded vox for reconstruction, 
% so we can insert missing values back into class_map
% -----------------------------------------------------------------------
obj.volInfo = meta_data_obj.volInfo;
obj.include_vox = meta_data_obj.include_vox;

obj.terms = meta_data_obj.terms;  % save for reporting later

end % function
