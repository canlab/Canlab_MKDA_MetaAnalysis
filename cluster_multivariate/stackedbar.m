function stackedbar(v,colors,labels)
% stackedbar(v,colors,labels)
% plots a stacked bar plot in order of magnitude, so that you can see all
% values for all variables

figure('Color','w'); set(gca,'FontSize',18); hold on;

if isempty(colors), colors = {'b' 'g' 'r' 'y' 'c' 'm'};, end
v(isnan(v)) = 0;

for k = 1:size(v,2)
    
for i = 1:size(v,1)

    w = find(v(i,:) == max(v(i,:)));
    
    if length(w) > 1
        % two values have same max
        % could change offset here 
        w = w(1);
        alph = .5;
        
    else
        alph = 1;
    end
    
    tmp = zeros(size(v,1),1);
    tmp(i) = max(v(i,:));
    
    h = barh(1:size(v,1) + .1*w,tmp,colors{w});
    %set(h,'FaceAlpha',alph);
    v(i,w) = -Inf;
    
end

end

set(gca,'YTick',1:size(v,1),'YTickLabel',labels,'YLim',[.5 size(v,1) + .5])

return