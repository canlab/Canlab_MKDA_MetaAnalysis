classdef meta_nbc
% meta_nbc: Naive Bayesian classifier object for meta-analysis
%
% An empty object is created by calling the function with no arguments.
% Properties below are assigned using the methods train, test, and cv
% along with meta_dataset objects.
%
% % test properties
% class_pred
% posterior_probs
% 
% % error properties
% p_err
% t_err
% b_err
%
% % data selection properties
% trainIdx 
% (note: select test dataset before running test method;
% testIdx not stored here)
% terms
% 
% % training properties
% n_classes
% p_act
% prior_probs
% train_pprobs_act
% train_pprobs_noact
% class_map
%
% % reconstruction/map properties
% volInfo
% include_vox
%
% parameters:
% k, num of cross-val folds, default = 10
% s_param, regularization for near-zero probabilities, default = 2
% w_act_noact, [0 - 1] weight for posteriors on activation vs.
% no-activation, default = .5 (standard naive bayes on both act and no-act)
% w_prior_vs_like, [0 - 1] weight for posteriors on priors vs.
% likelihood, default = 0 (standard naive bayes).  higher values will give
% more weight to priors, as though there were fewer independent measures, which may be more appropriate if there is positive
% dependence across voxels. Ideally, this parameter might be 1 / effective
% number of independent tests
%
% Full classification example for pain vs. touch:
% data = meta_dataset('MC_Setup', MC_Setup);
% data = setup_data(data, {'pain', 'touch'});
% nbc = train(meta_nbc, data);
% nbc = test(nbc, data.dat);
% nbc = get_error(nbc, data.classes);
% report(nbc);

    properties
        
        params = struct('k', 10, 's_param', 2, 'predict_pval_thresh', .05, 'w_act_noact', .5, 'w_prior_vs_like', 0);

        % test properties
        class_pred
        posterior_probs
        test_class
        
        % error properties
        p_err
        t_err
        b_err
        
        % data selection properties
        trainIdx
        terms
        
        % training properties
        n_classes
        p_act
        prior_probs
        train_pprobs_act
        train_pprobs_noact
        class_map
        p_a_given_c
        p_noa_given_c
        p_a_given_no_c
        p_noa_given_no_c
        
        
        % reconstruction/map properties
        volInfo
        include_vox
        
        % not needed
        %train_pred_fcn = @(trainX, trainclass, testX) nbc_predfun(trainX,
        %trainclass, testX);
        %test_pred_fcn  = @(trainX, trainclass, testX) nbc_predfun(trainX, trainclass, testX);
        
        loss_fcn = @(train_maps, train_class, test_maps, test_class) nbc_errfun(train_maps, train_class, test_maps, test_class)
        
    end
    
end
