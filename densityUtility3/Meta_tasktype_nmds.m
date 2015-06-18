function Meta_tasktype_nmds(Xi,dat)
    
ntasks = size(Xi,2);
[nstudies,nvox] = size(dat);


% get average activation of each task
taskvec = zeros(ntasks,nvox);

    for i = 1:ntasks
        
        wh = find(Xi(:,i));
        taskvec(i,:) = mean(dat(wh,:));
        
    end
    
    taskvec = taskvec';
    
    % NMDS
    coord = 
    
    
    
    return
    
    