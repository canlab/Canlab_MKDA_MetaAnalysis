function weighted_sum_image(P,Pout,w)
%weighted_sum_image(input file names,output file name);
%
% Creates a weighted sum of images based on weights for each image
%
% weighted_sum_image(DB.PP,'Activation.img',DB.studyweight);

switch spm('Ver')
    case 'SPM2'
        % spm_defaults is a script
        disp('WARNING: spm defaults not set for spm2. Make sure your defaults are set correctly');

    case 'SPM5'
        % spm_defaults is a function
        spm_defaults()

    case 'SPM8'
        % spm_defaults is a function
        spm_defaults()

    otherwise
        % unknown SPM
        disp('Unknown version of SPM!');
        spm_defaults()
end
defaults.analyze.flip = 0;

% initialize output
V = spm_vol(deblank(P(1,:)));
V.fname = Pout;

dat = zeros(V.dim(1:3));
datu = dat;                 % unweighted counts

% read and weight

for i = 1:size(P,1)
    v = spm_read_vols(spm_vol(deblank(P(i,:))));
    
    dat = dat + v .* w(i);
    
    % unweighted counts
    datu = datu + double(v > 0);
    
end

% write output
spm_write_vol(V,dat);

V.fname = [Pout(1:end-4) '_counts.img'];
spm_write_vol(V,datu);

return