function [odm,clusters] =  i_density_diff(XYZ,XYZ2,radius,niterations,varargin)
% function [odm,clusters] =  i_density_diff(XYZ,XYZ2,radius,niterations,[resume at iteration])
%
% iterated density algorithm, to mask out significant areas from
% density analysis and re-run on remaining points, until convergence
% this is a step-down test.
%
% XYZ in mm coordinates
% radius of smoothing sphere in mm
%
% odm is cumulative (overall) thresholded density mask
% this is the union of all previous density masks.
%
% if this quits in the middle, load the last saved working.mat file
% input the XYZ list saved there, which is the one for the iteration
% that quit.  
% take the union of the results you get after restarting and the 
% saved odm, and re-make clusters using mask2clusters 
% 
% by Tor Wager
%
% Example:
% XYZ = [randn(5,3).* 20;rand(3,3).*20+20;rand(3,3).*20 - 20];
% XYZ = [XYZ; randn(80,3) - repmat([40 40 0],80,1)];
% XYZ = [XYZ; randn(12,3) - repmat([30 -30 20],12,1)];
%XYZ2 = [randn(5,3).* 20;rand(3,3).*20+20;rand(3,3).*20 - 20];
%XYZ2 = [XYZ2; randn(12,3) + repmat([30 30 20],12,1)];
%XYZ2 = [XYZ2; randn(80,3) + repmat([40 -40 0],80,1)];
% figure('Color','w'); xyz = XYZ;plot3(xyz(:,1),xyz(:,2),xyz(:,3),'ro'); hold on; addbrain
% plot3(XYZ2(:,1),XYZ2(:,2),XYZ2(:,3),'bs');
% view(0,90),camzoom(1.3),drawnow
% [tmp,tmpcl] = i_density_diff(XYZ,XYZ2,10,5);    % 10 mm radius, 5 iterations




n = size(XYZ,1); n2 = size(XYZ2,1);
converged = 0;
indx = 1;

P = ['brain_avg152T1.img'];     % 2 mm voxels; don't change this filename w/o changing density.m as well.
P = which(P); V = spm_vol(P);

if length(varargin) > 0 & ~ischar(varargin{1})
    % -----------------------------------------------------------------------
	% Code to resume at a specified iteration
	% -----------------------------------------------------------------------
    indx = varargin{1};

    % load mask/XYZ
    load(['i_dens_working_step' num2str(indx-1)])
    disp(['loaded: i_dens_working_step' num2str(indx-1)])
    
    % -----------------------------------------------------------------------
	% smooth odm and save smoothed overall mask (odms) in tmp_brainmask
	% -----------------------------------------------------------------------
    vmask = spm_read_vols(V);	% read the FULL brain mask 
    V = spm_vol('tmp_brainmask.img'); V.M = V.mat;
    vox_radius = radius ./ mean(diag(V.mat(1:3,1:3)));
            
    tmp = zeros(size(odm)); tmp(round(size(odm,1)./2),round(size(odm,2)./2),round(size(odm,3)./2)) = 1;
    tmp = tor_conv_sphere(tmp,vox_radius,'sphere'); thresh = 1./max(tmp(:));
    odms = tor_conv_sphere(odm,vox_radius,'sphere',thresh);
    odms(odms > 0) = 1; % threshold to get mask of all results, smoothed
    vmask = vmask & ~odms;		% mask out the smoothed density mask    

    V.fname = 'tmp_brainmask.img'; spm_write_vol(V,vmask);
    P = 'tmp_brainmask.img';
    
    % -----------------------------------------------------------------------
	% get rid of points in XYZ in mm list that are in smoothed odms
	% -----------------------------------------------------------------------
    
    XYZout = mm2voxel(XYZ,V,1);
    XYZind = sub2ind(size(odms),XYZout(:,1),XYZout(:,2),XYZout(:,3));
    whrm = odms(XYZind); whrm(whrm > 0) = 1;
    XYZ(find(whrm),:) = [];
    
    disp(['Starting at iteration ' num2str(indx) ': original XYZ is ' num2str(size(XYZout,1)) ' new is ' num2str(size(XYZ,1))])
    disp([num2str(sum(odm(:) > 0)) ' Significant voxels so far, ' num2str(sum(vmask(:) > 0)) ' left in brain mask'])
    clear vmask
    n = size(XYZ,1);
    skip1st = 0;
    
elseif length(varargin) > 0 & ischar(varargin{1})
    % resume by loading mat file from previous density_diff analysis
    odm = zeros(V.dim(1:3));
    disp(['Loading ' varargin{1}])
    load(varargin{1})
    if ~exist('OUT')==1, error('Variable OUT not found in loaded mat file.'), end
    skip1st = 1;
    
else
    % no resume, normal function
	odm = zeros(V.dim(1:3));
    skip1st = 0;
end



diary i_density_log.txt

while ~converged
	
	disp(' '); disp(['Starting iteration ' num2str(indx)])
	disp('_________________________________________________')

	% -----------------------------------------------------------------------
	% compute density threshold and apply to the XYZ points
	% -----------------------------------------------------------------------

    if ~skip1st
        [u,OUT] = density_diff_pdf(n,n2,radius,niterations,['i_density_rad' num2str(radius) '_step' num2str(indx)],P); close all;
        u = OUT.du; u2 = OUT.du2;    
    else
        u = OUT.du; u2 = OUT.du2;    % we already have it
        skip1st = 0;
        disp('Skipping to iteration 2')
    end
    
    v = 4*pi*radius^3 / 3; p = u .* v;
    disp(['Threshold is at least ' num2str(p) ' points difference within r = ' num2str(radius) ' mm with a volume of ' num2str(v) ' mm3'])
    
    diary off, vv = spm_vol(P);vol = spm_read_vols(vv);, diary on
    fprintf(1,'\nMask contains %5.0f voxels, XYZ contains %3.0f points, doing %5.0f iterations this step', sum(vol(:)>0),n,niterations)
    clear vol
    
    %[u,maxd] = fast_max_density(n,niterations,radius,P);
    eval(['save i_density_rad' num2str(radius) '_step' num2str(indx) ' u OUT n n2 niterations radius indx']) 
    
    [cl1,cl2,dmt,dmt2,d1,d2,mask,mask2] = density_diff(XYZ,OUT.du,XYZ2,OUT.du2,radius,['id_a_rad' num2str(radius) '_step' num2str(indx)],['id_b_rad' num2str(radius) '_step' num2str(indx)]);
    %[dmt,clusters,dm,mask] =  density(XYZ,u,radius,['i_density_rad' num2str(radius) '_step' num2str(indx)],P);,close
	odm = odm | dmt;

	fprintf(  1,'\nThresh is %3.6f, %3.0f clusters, %5.0f voxels significant, cumulative %5.0f sig voxels.\n',u,length(cl1),sum(dmt(:) > 0),sum(odm(:) > 0)  )

    v = 4*pi*radius^3 / 3;
    p = u .* v;

    disp(['Threshold is at least ' num2str(p) ' points difference within r = ' num2str(radius) ' mm with a volume of ' num2str(v) ' mm3'])

    
	% mask is the mask of points
	% dmt is thresholded density mask

    % -----------------------------------------------------------------------
    % make a "smoothed" or convolved version of odm, so that 
    % points w/i a radius of significant voxels are removed from both
    % mask and point list.  If you don't do this, the significant
    % voxels may be between points, and so the points that contributed
    % to the significant region will not be removed.  
    % -----------------------------------------------------------------------
    
    vox_radius = radius ./ mean(diag(V.mat(1:3,1:3)));
    
    tmp = zeros(size(odm)); tmp(round(size(odm,1)./2),round(size(odm,2)./2),round(size(odm,3)./2)) = 1;
    tmp = tor_conv_sphere(tmp,vox_radius,'sphere'); thresh = 1./max(tmp(:));
    odms = tor_conv_sphere(odm,vox_radius,'sphere',thresh);
    odms(odms > 0) = 1; % threshold to get mask of all results, smoothed
    
	% -----------------------------------------------------------------------
	% save a new brainmask file with smoothed odm removed
	% -----------------------------------------------------------------------
	disp('Updating brain mask')
	v = spm_read_vols(V);
	v(odms>0) = 0; V.fname = 'tmp_brainmask.img';
	spm_write_vol(V,v);
	P = V.fname;	% next time around it will load the new mask


	% -----------------------------------------------------------------------
	% get rid of points and convert  back to XYZ in mm
	% -----------------------------------------------------------------------
    
	mask(find(odms)) = 0;       % get rid of any points in overall smoothed dens mask

	X = []; Y = []; Z = [];
	for i = 1:max(mask(:))
		[x y z] = ind2sub(size(mask),find(mask>=i));
		X = [X;x]; Y = [Y;y]; Z = [Z;z];
	end
		
	if isempty(X)
		newmm = [];
	else
		newmm = voxel2mm([X Y Z]',V.mat)'; 
	end

    mask2(find(odms)) = 0;       % get rid of any points in overall smoothed dens mask

	X = []; Y = []; Z = [];
	for i = 1:max(mask2(:))
		[x y z] = ind2sub(size(mask2),find(mask2>=i));
		X = [X;x]; Y = [Y;y]; Z = [Z;z];
	end
		
	if isempty(X)
		newmm2 = [];
	else
		newmm2 = voxel2mm([X Y Z]',V.mat)'; 
	end
    
	%mostmoved = max(max(sort(newmm) - sort(XYZ)));
	%disp(['New mm coordinates found.  Moved in transformation at most by ' num2str(mostmoved) ' mm'])
	
	% -----------------------------------------------------------------------
	% check for convergence and save stuff if yes.
	% -----------------------------------------------------------------------

	% if no points masked out, then we converge.
	if length(cl1) == 0 | (isempty(newmm) & isempty(newmm2))
		disp(['Converged at iteration ' num2str(indx)])
	
		converged = 1;
		VV = V;  VV.fname = 'i_density_res_mask';
		spm_write_vol(VV,odm);
        
        VV.fname = 'i_density_smoothed_res_mask';
		spm_write_vol(VV,odms);
        
		if sum(odm(:)) > 0, 
            clusters = mask2clusters('i_density_res_mask');
            clusters_smoothed = mask2clusters('i_density_smoothed_res_mask');
            save i_density_clusters clusters clusters_smoothed
        else
            clusters = [];
		    disp('No significant results!')
        end

	else
	
		XYZ = newmm; XYZ2 = newmm2;
        n = size(XYZ,1);
        n2 = size(XYZ2,1);
		eval(['save i_dens_working_step' num2str(indx) ' XYZ XYZ2 n n2 odm u'])
        
		indx = indx + 1;
		
	end

end 	% end while loop


return
