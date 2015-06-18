function [chi2, df, p, sig, warn, tab, isnonpar] = Meta_Parcel_chisq(MC_Setup, parcels, parcel_stats) 
%
% Documentation here.
%
% [chi2, df, p, sig, warn, tab, isnonpar] = Meta_Parcel_chisq(MC_Setup, parcels, parcel_stats) 
%
% First run Meta_Activation_FWE('setup'....) which creates MC_Setup (stored in MC_Info.mat)
% Then run Meta_Parcel, which gives you parcels and parcel_stats
% Then run this.

%% Chi-square
Xi = MC_Setup.Xi;

w = MC_Setup.wts;
w = w ./ mean(w);

studybyparcel = parcel_stats.studybyparcel;

[chi2,df,p,sig,warn,tab,expected,isnonpar] = meta_analyze_data(studybyparcel,'X',Xi,'chi2','w',w, 'table');

%%
diary Meta_parcel_chisq_output.txt
nparcels = length(parcels);
spm_defaults

for i = 1:nparcels
    
    parcel_number = i;

sig_chi2 = p < .05;

create_figure('Barplot', 1, 1); 

% % create_figure('Barplot', 1, 3); 
cnt = Xi' * studybyparcel(:, parcel_number); 
% % bar(cnt, 'FaceColor', [.5 .5 .5]);
% % set(gca, 'XTick', 1:length(cnt), 'XTickLabel', MC_Setup.Xinms);
% % title('Unweighted counts');
% % 
% % wcnt = tab{parcel_number}(2, :);
% % subplot(1, 3, 2);
% % bar(wcnt, 'FaceColor', [.5 .5 .5]);
% % set(gca, 'XTick', 1:length(cnt), 'XTickLabel', MC_Setup.Xinms);
% % title('Weighted');

propact = tab{parcel_number}(2, :) ./ sum(tab{parcel_number});
% % subplot(1, 3, 3);
bar(propact, 'FaceColor', [.5 .5 .5]);
set(gca, 'XTick', 1:length(cnt), 'XTickLabel', MC_Setup.Xinms);
title('Weighted Proportion')

saveas(gcf, ['parcel_' num2str(parcel_number)], 'png');


create_figure('Brainslice'); montage_clusters_maxslice([], parcels(parcel_number), {'r'});
str = sprintf('Parcel %3.0f, chi2 = %3.2f, p = %3.6f, nonpar = %01d', parcel_number, chi2(parcel_number), p(parcel_number), isnonpar(parcel_number));
title(str, 'FontSize', 14);

saveas(gcf, ['parcel_region' num2str(parcel_number)], 'png');


disp(str);

 
end
diary off

end
