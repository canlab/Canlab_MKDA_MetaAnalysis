function [p_ste_avg,iterations] = meta_simulate_nonparamchi2(respfreq,numconds,N)
% [p_ste_avg,iterations] = meta_simulate_nonparamchi2(respfreq,numconds,N)
%
% Perform tests on made-up data to determine what the variability in
% nonparametric chi2 p-value estimates is as a function of the number
% of iterations in the nonparam_chi2 test.
%
% tor wager
%
% e.g., to set up data
%
%respfreq = .2;      % frequency of 'hits' (yesses)
%numconds = 2;       % number of conditions (nominal)
%N = 50;             % number of observations


[y,condf] = generate_data(respfreq,numconds,N);


% repeat to get STEs

iterations = [50 100 200 400 800 1600 3200 6400 12800]; %[500 1000 5000 10000]; %50000];

reps_within_dataset = 30;
dataset_reps = 30;

for i = 1:length(iterations)

    iter = iterations(i);

    e1 = clock;
    fprintf(1,'\nIterations: %3.0f  rep %03d',iter,0);

    dataindex = 1;  % keeps track of indices of reps that belong to the same dataset

    for dataindex = 1:dataset_reps

        for j = 1:reps_within_dataset


            [realchi2,chi2p(j,1)] = nonparam_chi2(y,condf,iter);

            fprintf(1,'\b\b\b%03d',reps_within_dataset.*(dataindex-1) + j );
        end

        [y,condf] = generate_data(respfreq,numconds,N);
        pste(dataindex,i) = std(chi2p);



    end

    p_ste_avg = mean(pste);

    p_ste_ste = ste(pste);  % std. error across datasets (realizations)
    
    e2 = etime(clock,e1);

    fprintf(1,' Time: %3.0f s  \n',e2);

end

% plot data

%tor_fig; plot(iterations,p_ste_avg,'ko-','LineWidth',2); hold on;

tor_fig;
tor_fill_steplot(pste,'k');

return





function [y,condf] = generate_data(respfreq,numconds,N)

num_ones = round(N.*respfreq);
y = zeros(N,1);     % data

indx = randperm(N)';
y(indx(1:num_ones)) = 1;

% equal frequencies in each condition
npercond = ceil(N./numconds);
condf = [];
for i=1:numconds

    condf = [condf; i*ones(npercond,1)];

end

return

