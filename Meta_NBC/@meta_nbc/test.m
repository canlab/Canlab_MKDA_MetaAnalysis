function obj = test(obj, test_data, varargin)
% obj = test(obj, test_data, varargin)
%
% test_data should be voxels x maps matrix of test data
% can be either binary (i.e., meta-analysis peaks) 
% or continuous (e.g., single-subject t- or contrast maps)
% Different methods are used for each data type.
%
% 

dobinarize = 0;
if nargin > 2, dobinarize = varargin{1}; end

if isempty(obj.train_pprobs_act)
    fprintf('You must train the nbc object first: train(nbc_obj, data_obj)\n')
    error('Exiting');
end

% -----------------------------------------------------------------------
% Sum of log posterior probabilities of each study 
% predictions for the test set
% -----------------------------------------------------------------------

if islogical(test_data) || length(unique(test_data(:))) < 3
    % we have a binary map
    fprintf('Predicting binary test data\n')
    [class_pred, posterior_probs] = predict_binary_images(logical(test_data), obj);
    
else
    % assume continuous-valued images unless 'binarize' is set, in
    % which case, use a threshold of p < .05 and a minimum cut-off of 1%
    % of voxels
    if dobinarize
        thresh = 1.65;
        min_sig = 0.01;
        nvox = size(test_data,1);
        for i=1:size(test_data,2)
            active = zeros(nvox, 1);
            sig_vox = (test_data(:,i) >= thresh);
            nsig = sum(sig_vox)/nvox;
            [sor, sindex] = sort(test_data(:,i), 1, 'descend');
            if nsig < min_sig, nsig = min_sig; end
            active(sindex(1:(round(nsig*nvox))))=1;
            test_data(:,i) = active;

        end
        fprintf('Binarizing continuous-valued test data\n') 
        [class_pred, posterior_probs] = predict_binary_images(logical(test_data), obj);
    else
        fprintf('Predicting continuous-valued test data\n')   
        [class_pred, posterior_probs] = predict_continuous_images(test_data, obj);
    end
    
end

obj.class_pred = class_pred;
obj.posterior_probs = posterior_probs;

fprintf('Updated class_pred, posterior_probs in meta_nbc object.\n');


end % function





function [class_pred, posterior_probs] = predict_binary_images(test_data, obj)

[n_vox, n_maps] = size(test_data);
n_classes = obj.n_classes;
 
if isempty(n_classes), error('obj.n_classes is empty.'); end

posterior_probs = zeros(n_maps, n_classes);

fprintf('Map %3.0f', 0);
for ii = 1:n_maps
    fprintf('\b\b\b%3.0f', ii);
    
    % construct weights for voxels, based on two parameters
%     nbc.params.w_act_noact:    weight for activation vs. lack of activation, 0-1  
%     if overall prob(activation) varies among classes, we can
%     get strange results (always classifying as the least activating
%     class) unless we downweight no-act.  
%     nbc.params.w_prior_vs_like: weight for priors vs. likelihood; 0 is
%     standard NBC classifier, 1 divides likelihood by # voxels, so weights
%     priors and likelihood equally

    w = ones(n_vox, 1, 'single');
    w(test_data(:, ii)) = 2 * obj.params.w_act_noact; % scale so that param .5 gives standard sum log like (w all ones)
    w(~test_data(:, ii)) = 2 * (1 - obj.params.w_act_noact);
    
    % scale so that weighted sum gives something between standard log like
    % and normalization by number of voxels
    w = w ./ (1 + (obj.params.w_prior_vs_like * (n_vox - 1)));
    
    % post prob map: voxels x classes matrix of posterior probabilities
    % given either activation or lack thereof
    [p_c_map, p_no_c_map] = deal(zeros(size(obj.p_a_given_c), 'single'));

    p_c_map(test_data(:, ii), :) = obj.p_a_given_c(test_data(:, ii), :);
    p_c_map(~test_data(:, ii), :) = obj.p_noa_given_c(~test_data(:, ii), :);
    
    p_no_c_map(test_data(:, ii), :) = obj.p_a_given_no_c(test_data(:, ii), :);
    p_no_c_map(~test_data(:, ii), :) = obj.p_noa_given_no_c(~test_data(:, ii), :);

    %posterior_probs(ii, :) = sum(log(p_c_map./p_no_c_map)); % divide in log space for numerical stability
    pp = log(p_c_map) - log(p_no_c_map);
    posterior_probs(ii, :) = w' * pp;   % weighted sum
end
fprintf('\n');

[discard, class_pred] = max(posterior_probs');
class_pred = class_pred';

% Now, the thing is that posterior probs are not scaled nicely because of
% the assumption of independence across voxels. I'm not sure how to correct
% for this -- but we do want to do that if we want interpretable
% probability numbers.  This also might be useful for defining a
% posterior probability-based error metric, which should be better than
% the misclassification rate for model/parameter selection.

n_vox = size(test_data, 1);
pp = exp(posterior_probs ./ n_vox);
pp = pp ./ repmat(sum(pp, 2), 1, size(pp, 2)); % sum to 1 across classes (this may mainly fix rounding error??)
posterior_probs = pp;  % this produces the same class prob estimates as before.
% these all tend to be very close to the base rate, which may be a product
% of the fact that the base rate goes into the post. prob calc at each
% voxel. may not matter for picking most likely class?

end



function [class_pred, posterior_probs] = predict_continuous_images(test_data, obj)

% now we have to make predictions about new images.
% using z-scores based on p(class | activation) provides a natural
% weighting.  images that we're classifying might be z-scores from
% individual subjects, or they might be something else.
% so, rather than binarizing the test image (and thus ignoring information
% about de-activations, and imposing an arbitrary threshold, etc.)
% let's just take a weighted sum over the cross-products of probabilities
% and normed test images
% 
% using this norm: (sum(test_data .^ 2)) .^ .5 == 1
for ii = 1:size(test_data, 2)
    test_data(:, ii) = test_data(:, ii) ./ norm(test_data(:, ii));
end

% test images x classes matrix of class posterior_probs -
% this is not scaled correctly, but should have relative rank info...
% ideas about how to improve??
posterior_probs = (obj.train_pprobs_act' * test_data)';

[discard, class_pred] = max(posterior_probs');
class_pred = class_pred';

end
