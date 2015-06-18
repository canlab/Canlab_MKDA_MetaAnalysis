function SOMResults_view_connections(studyByCluster,cl,wh)

% function SOMResults_view_connections(studyByCluster,cl,wh)

% a) find clusters significantly connected to this one
wh_conx = studyByCluster.sigbonf(wh,:); wh_conx(wh) = 0;
wh_conx = find(wh_conx);

cluster_orthviews(cl{wh},{[0 0 0]}); %studyByCluster.colors(wh));

if isempty(wh_conx)
    
    disp('No significantly connected clusters.')
    
else
    
    fprintf(1,'Found: %3.0f correlated clusters.\n',length(wh_conx));

    colors = studyByCluster.colors(wh_conx);
  
    for i = 1:length(wh_conx)
    
        cluster_orthviews(cl{wh_conx(i)},colors(i),'add');
    
    end
    
end

return

