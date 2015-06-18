function plot_sens_spec(obj)
% plot_sens_spec(obj)
% 
% Barplots of sensitivity and specificity for meta_nbc object.
% You should have tested the nbc already: e.g., obj = test(obj); 

nterms = length(obj.terms);

colors = {'y' 'b' 'r' 'g' [1 .5 0] [1 0 .5] [0 1 .5] [.5 1 0] [0 .5 1] [.5 0 1]};
while length(colors) < nterms, colors = [colors colors]; end

names = obj.terms;


[cm, dprime, hr, far, misclass]  = confusion_matrix(obj.test_class, obj.class_pred);  
disp(cm)

for i = 1:nterms, h(i) = bar(i, hr(i)); set(h(i), 'FaceColor', colors{i}); end
for i = 1:nterms, h(i) = bar(1+i+nterms, 1 - far(i)); set(h(i), 'FaceColor', colors{i}); end

set(gca, 'FontSize', 24)

set(gca, 'XTick', 1:2*nterms+1, 'XTickLabel', [names {' '} names]);
ylabel('Probability');
title('Prediction sensitivity and specificity', 'FontSize', 28);

dat = [hr 0 1 - far];
for i = 1:nterms, n(i) = sum(obj.test_class == i); end

errorbars = sqrt(dat .* (1 - dat) ./ [n Inf n]);
tor_bar_steplot(dat, errorbars, {'k'});

han = plot_horizontal_line(1./nterms); set(han, 'LineStyle', '--');

end % function 