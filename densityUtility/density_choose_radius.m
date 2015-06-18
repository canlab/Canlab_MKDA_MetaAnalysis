function density_choose_radius(XYZ,radii,niterations)
% function density_choose_radius(XYZ,radii,niterations)
%
% XYZ is 3-column vector of coordinates
% radii is a vector of radius values in mm for density analysis
% niterations is number of iterations per Monte Carlo simulation
%
% This function takes an XYZ list of brain activation coordinates
% and runs a series of nonparametric thresholds with smoothing kernels
% of different radii (see density.m and density_pdf.m)
% in order to choose an appropriately sensitive radius for the
% density analysis.
%
% it iteratively runs density_pdf.m

n = size(XYZ,1);

for radius = 1:length(radii)

	disp(['Starting simulation for radius ' num2str(radii(radius))])

	%[u(radius),OUT] = density_pdf(n,niterations,radii(radius),['tmp_dens_rad' num2str(radii(radius))]);
    [u(radius),maxdd] = fast_max_density(n,niterations,radii(radius));
    eval(['save dens_choose_rad_' num2str(radius) ' u maxdd n niterations radius radii']) 
    
	[dmt,clusters,dm] =  density(XYZ,u(radius),radii(radius),['tmp_pts_rad' num2str(radii(radius))]);
	eval(['save dens_rad' num2str(radii(radius)) '_clusters clusters'])
	nclust(radius) = length(clusters);
	maxd(radius) = max(dm(:));
	close all

	P = ['dens_tmp_pts_rad' num2str(radii(radius)) '_filtered.img'];
	try
		clusters = smooth_dens_and_mask(P);
		nclusts(radius) = length(clusters);
	catch
		nclusts(radius) = 0;
	end

end

try
figure('Color','w'); subplot 121; 
plot(radii,u,'ro-','LineWidth',2), hold on; plot(radii,maxd,'bs-','LineWidth',2)
title('Height threshold and max density by radius')

subplot 122; plot(radii,nclust,'ko-','LineWidth',2),
plot(radii,nclusts,'ro-','LineWidth',2), legend({'sig clusters' 'sig clusters 6 mm apart'})
title('Significant clusters by radius')
saveas(gcf,'choose_radius_output','tif')
save choose_radius_output nclust nclusts maxd u radii
catch
    disp('Cannot make figure')
end

return


