function Meta_Chisq(DB,whfield,levels)
% Meta_Chisq(DB,whfield)
%
%DB.PP       % list of image names

%DB.XYZmm    % mm coords for each study
%DB.maskV    % mask mapped volume
%DB.radius   % in voxels


% -----------------------------------------------------
% across all density images, compute prob(activation), p(a)
% -----------------------------------------------------
v = spm_read_vols(spm_vol(DB.PP));
pa = sum(v,4);

% -----------------------------------------------------
% for each type
% -----------------------------------------------------
myvar = DB.(whfield);
myvar = myvar(DB.pointind);

% select contrasts with levels of interest

% include = ones(size(DB.Contrast));
% for i = 1:length(DB.pointind), 
%     test = strmatch(myvar{i},levels); 
%     if isempty(test), include(DB.Contrast==i) = 0;,end
% end
  
include = ones(size(DB.pointind));
for i = 1:length(DB.pointind), 
    test = strmatch(myvar{i},levels); 
    if isempty(test), include(i) = 0;,end
end

v = v(:,:,:,find(include));
myvar = myvar(find(include));
 
for i = 1:length(levels)
    
    wh = find(strcmp(myvar,levels{i}));
    
    palev{i} = sum(v(:,:,:,wh),4);
    
    V = DB.maskV; V.fname = ['Activation_' levels{i} '.img']; spm_write_vol(V,palev{i});

    
end
  
% try Fisher's Exact test...
    
% -----------------------------------------------------
% get expected null-hypothesis probability
% -----------------------------------------------------
[p0,pa_h0] = p0_monte_carlo_FWE(DB,10);


pa(pa <= p0) = 0;
V = DB.maskV; V.fname = 'Activation_thresholded.img'; spm_write_vol(V,pa);

