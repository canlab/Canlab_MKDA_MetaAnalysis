function [odm,clusters] =  i_density(XYZ,radius,niterations,varargin)
% function [odm,clusters] =  i_density(XYZ,radius,niterations,[resume at iteration])
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
% figure('Color','w'); xyz = XYZ;plot3(xyz(:,1),xyz(:,2),xyz(:,3),'ro'); hold on; addbrain
% view(0,90),camzoom(1.3),drawnow
%
% add z-scores
% XYZ(:,4) = rand(size(XYZ,1),1)
% [tmp,tmpcl] = i_density(XYZ,10,5);    % 10 mm radius, 5 iterations
%
% FOR PLOTTING, CHANGE HARD-CODED FLAG IN THIS FUNCTION

doplot = 1;

n = size(XYZ,1);

% format n for z-scores, if they are entered
if size(XYZ,2) > 3
    disp('Z-scores entered in XYZ column 4.  Using z-score weighting. WILL NOT WORK FOR RESUME RIGHT NOW')
    disp('ALSO, PROBLEMS COULD ARISE IF SOME POINTS HAVE THE EXACT SAME COORDINATES...CHECK THIS FUNCTION IF SO.')
    n = XYZ(:,4);
    XYZ = XYZ(:,1:3);
end

converged = 0;
indx = 1;

P = ['brain_avg152T1.img'];     % 2 mm voxels; don't change this filename w/o changing density.m as well.
P = 'scalped_avg152T1_graymatter_smoothed.img';
P = which(P); V = spm_vol(P);

if length(varargin) > 0
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
            
    odm = double(odm);
    tmp = zeros(size(odm)); tmp(round(size(odm,1)./2),round(size(odm,2)./2),round(size(odm,3)./2)) = 1;
    tmp = tor_conv_sphere(tmp,vox_radius,'sphere'); thresh = 1./max(tmp(:));
    odms = tor_conv_sphere(odm,vox_radius,'sphere',thresh);
    odms(odms > 0) = 1; % threshold to get mask of all results, smoothed
    vmask = double(vmask & ~odms);		% mask out the smoothed density mask    

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
    
else
    % no resume, normal function
	odm = zeros(V.dim(1:3));
end

   if doplot, 
      tor_fig(1,2); subplot(1,2,1); hold on; 
      tmpv = spm_read_vols(spm_vol(P));
      tmp = cat(1,sum(sum(tmpv))); % find max slice
      whz = find(tmp == max(tmp)); whz = whz(1);
      imagesc(flipud(rot90(tmpv(:,:,whz)))); 
      title('Search mask'), colormap gray,axis off; axis image,drawnow
      
      subplot(1,2,2); hold on; 
      plot3(XYZ(:,1),XYZ(:,2),XYZ(:,3)+100,'ko','MarkerFaceColor','k');
      addbrain; view(0,90)
      title('Peaks MIP'),drawnow
   end

diary i_density_log.txt

while ~converged
	
	disp(' '); disp(['Starting iteration ' num2str(indx)])
	disp('_________________________________________________')

	% -----------------------------------------------------------------------
	% compute density threshold and apply to the XYZ points
	% -----------------------------------------------------------------------

    % old, slow way.
	%[u,OUT] = density_pdf(n,niterations,radius,['i_density_rad' num2str(radius) '_step' num2str(indx)],P); close all;
	%u
    %v = 4*pi*radius^3 / 3;
    %p = u .* v;
    %disp(['OLD way: Threshold is at least ' num2str(p) ' points within r = ' num2str(radius) ' mm with a volume of ' num2str(v) ' mm3'])
    %figure;hist(OUT.maxd)
    % testing
    %spm_image('init',P),
    
    
    diary off, vv = spm_vol(P);vol = spm_read_vols(vv);, diary on
    npts = n; if length(npts)>1, npts = length(npts);,end
    fprintf(1,'\nMask contains %5.0f voxels, XYZ contains %3.0f points, doing %5.0f iterations this step', sum(vol(:)>0),npts,niterations)
    clear vol

    % -----------------------------------------------------------------------
	% compute density threshold 
	% -----------------------------------------------------------------------
    
    [u,maxd] = fast_max_density(n,niterations,radius,P);    % n is vector of z-scores if entered
    eval(['save i_density_rad' num2str(radius) '_step' num2str(indx) ' u maxd n niterations radius indx']) 
    
    if doplot, 
        tor_fig;
        hist(maxd,round(length(maxd)./10)); hold on; plot([u u],[-1 1],'k','LineWidth',3),drawnow
    end
    
    % -----------------------------------------------------------------------
	% compute density map of image 
	% -----------------------------------------------------------------------
    if length(n) > 1    % we have z-scores
        [dmt,clusters,dm,mask] =  density(XYZ,u,radius,['i_density_rad' num2str(radius) '_step' num2str(indx)],n);,close
    else    % no z-scores
        [dmt,clusters,dm,mask] =  density(XYZ,u,radius,['i_density_rad' num2str(radius) '_step' num2str(indx)]);,close
    end
    
    if doplot, 
      tor_fig(1,2); subplot(1,2,1); hold on; 
      whz = flipud(rot90(max(dm,[],3)));
      %tmp = cat(1,sum(dmt,3); % find max slice
      %whz = find(tmp == max(tmp)); whz = whz(1);
      %imagesc(dmt(:,:,whz)); colorbar;
      imagesc(whz); 
       title('Density MIP'),axis off; axis image,drawnow,colorbar;

      subplot(1,2,2); hold on; 
      whz = flipud(rot90(max(dmt,[],3)));
      %tmp = cat(1,sum(dmt,3); % find max slice
      %whz = find(tmp == max(tmp)); whz = whz(1);
      %imagesc(dmt(:,:,whz)); colorbar;
      imagesc(whz);
       title('Results MIP'),axis off; axis image,drawnow, colorbar;
    end
    
    odm = double(odm | dmt);    % overall density mask, saves latest thresholded values total

	fprintf(  1,'\nThresh is %3.6f, %3.0f clusters, %5.0f voxels significant, cumulative %5.0f sig voxels.\n',u,length(clusters),sum(dmt(:) > 0),sum(odm(:) > 0)  )

    v = 4*pi*radius^3 / 3;
    p = u .* v;

    disp(['Threshold is at least ' num2str(p) ' points within r = ' num2str(radius) ' mm with a volume of ' num2str(v) ' mm3'])

    
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
    % get rid of z-scores too -- if entered
	% -----------------------------------------------------------------------
    
	mask(find(odms)) = 0;       % get rid of any points in overall smoothed dens mask

	X = []; Y = []; Z = [];
	
    if length(n) > 1
        % if we're doing z-scores, multiple points get added together into
        % one summed z-score.  this may cause inaccuracies in subsequent
        % steps if 2 + points fall on the same exact coordinate, because we
        % can't recover which point was which.
        [X Y Z] = ind2sub(size(mask),find(mask));
        n = mask(find(mask));
    else
        
        for i = 1:max(mask(:))   
            % this is to preserve duplicate coordinates
            [x y z] = ind2sub(size(mask),find(mask>=i));
            X = [X;x]; Y = [Y;y]; Z = [Z;z];
        end
    end
    
	if isempty(X)
		newmm = [];
	else
		newmm = voxel2mm([X Y Z]',V.mat)'; 
	end

	%mostmoved = max(max(sort(newmm) - sort(XYZ)));
	%disp(['New mm coordinates found.  Moved in transformation at most by ' num2str(mostmoved) ' mm'])
	
	% -----------------------------------------------------------------------
	% check for convergence and save stuff if yes.
	% -----------------------------------------------------------------------

	% if no points masked out, then we converge.
	if length(clusters) == 0 | isempty(newmm)
		disp(['Converged at iteration ' num2str(indx)])
	
		converged = 1;
		VV = V;  VV.fname = 'i_density_res_mask.img';
		spm_write_vol(VV,odm);
        
        VV.fname = 'i_density_smoothed_res_mask.img';
		spm_write_vol(VV,odms);
        
		if sum(odm(:)) > 0, 
            clusters = mask2clusters('i_density_res_mask.img');
            clusters_smoothed = mask2clusters('i_density_smoothed_res_mask.img');
            save i_density_clusters clusters clusters_smoothed
        else
            clusters = [];
		    disp('No significant results!')
        end

	else
	
		XYZ = newmm;
        
		eval(['save i_dens_working_step' num2str(indx) ' XYZ odm u'])
        n = size(XYZ,1);
		indx = indx + 1;
		
	end

end 	% end while loop


return
