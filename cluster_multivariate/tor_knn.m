function OUT = tor_knn(data,clas,k)
% OUT = tor_knn(data,clas,k)
%
% tor wager
%
% k nearest neighbor classification of T. Denoeux
% with cross-validation
%
%
% Example:
% Create a fake dataset and test
% x = randn(100,3);  y = repmat((1:4)',25,1);
%x(find(y==1),1) = x(find(y==1),1) + 2;  % add some discriminability to first 2 classes
%x(find(y==2),2) = x(find(y==2),2) + 2;
%tor_knn(x,y,3);

%[gamm,alpha,err] = knndsfit(data,clas,k);
%[m,classest] = knndsval(data,clas,k,gamm,alpha,0,data);   % apparent misclassification on whole set

% PROBLEM: THIS ALGORITHM DEPENDS ON ORDER OF ENTRIES IN DATA AND CLASS!
%
% References:
% 
% T. Denoeux. A k-nearest neighbor classification rule based on 
%  Dempster-Shafer theory. IEEE Transactions on Systems, Man
%  and Cybernetics, 25(05):804-813, 1995.
%
% L. M. Zouhal and T. Denoeux. An evidence-theoretic k-NN rule with 
% parameter optimization. IEEE Transactions on Systems, Man and 
% Cybernetics - Part C, 28(2):263-271,1998.

% simple one, for now:
classest = simple_knn(data,clas,k);


[m,dprime,corr,far,misclass] = confusion_matrix(clas,classest);

fprintf(1,'Apparent misclassification\n')
m
dprime
corr
far

% cross-validate

p = .85;    % percentage of data to use as training set in cross-validation

for j = 1:max(clas)
    dat{j} = [data(clas == j,:) clas(clas == j)];
    tot(j) = size(dat{j},1);
    num(j) = round(size(dat{j},1) .* p);
end

if any(num==0) | any(num==tot), warning('Some classes have no test data or no trainin data!');, keyboard, end

iter = 200; mc = []; dprimec = []; corrc = []; farc = [];

for i = 1:iter

    % get random sample of data for train and test
    % class labels are appended in last column to keep track
    
    train = []; test = [];
    for j = 1:max(clas)
        wh = randperm(tot(j));
        train = [train; dat{j}(find(wh <= num(j)),:)];
        test = [test; dat{j}(find(wh > num(j)),:)];
    end
    
    %[tmp,classestc] = knndsval(train(:,1:end-1),train(:,end),k,gamm,alpha,0,test(:,1:end-1));
    [tmp,classestc] = simple_knn(train(:,1:end-1),train(:,end),k,test(:,1:end-1));
    
    % cross-validated stats for this sample
    [mc(:,:,i),dprimec(i,:),corrc(i,:),farc(i,:)] = confusion_matrix(test(:,end),classestc);
end


% output table

OUT.hit_rate_xval = nanmean(corrc);
OUT.hit_rate_xval_sd = nanstd(corrc);
OUT.far_xval = nanmean(farc);
OUT.far_xval_sd = nanstd(farc);
OUT.dprime_xval = nanmean(dprimec);
OUT.dprime_xval_sd = nanstd(dprimec);

OUT.confusion_xval_mean = mean(mc,3);
OUT.confusion_xval = mc;

OUT.misclass = misclass;
OUT.class = classest;

fprintf(1,'Cross-validated classification\n')
fprintf(1,'Class\tHit rate\tFalse Alarms\tA-prime\t\n')
for i = 1:max(clas)
    fprintf(1,'%3.0f\t%3.2f (%3.2f)\t%3.2f (%3.2f)\t%3.2f (%3.2f)\t\n',i,OUT.hit_rate_xval(i),OUT.hit_rate_xval_sd(i), ...
        OUT.far_xval(i),OUT.far_xval_sd(i),OUT.dprime_xval(i),OUT.dprime_xval_sd(i))
end

disp('Cross-validated confusion matrix (mean)')
OUT.confusion_xval_mean


return
    


return