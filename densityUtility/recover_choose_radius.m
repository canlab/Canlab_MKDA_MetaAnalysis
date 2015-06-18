function [numcl,rad,nclusts] = recover_choose_radius(varargin)
% [numcl,rad,nclusts] = recover_choose_radius([plot color])
%  
% input of an argument signals to plot.  color = input argument
%
% tor wager
% 
%

d = dir('dens*clusters.mat'); d = str2mat(d.name);
for i = 1:size(d,1),

	load(deblank(d(i,:)))
	
	numcl(i) = length(clusters);
	rad(i) = nums_from_text(d(i,:));
	
end

% Sort the output into the correct order, if necessary
tmp = [numcl; rad];
tmp2 = sortrows(tmp',2);
numcl = tmp2(:,1)';
rad = tmp2(:,2)';

if length(varargin) > 0

	plot(rad,numcl,varargin{1},'LineWidth',2)

end



d = dir('dens_tmp_pts_rad*filtered.img');
for i = 1:size(d,1),

    P = d(i).name;
	try
		clusters = smooth_dens_and_mask(P);
		nclusts(i) = length(clusters);
    catch
        % does this if empty mask
		nclusts(i) = 0;
    end
    
    aa = nums_from_text(d(i).name);
    rad(i) = aa(1);
end

% Sort the output into the correct order, if necessary
tmp = [nclusts; rad];
tmp2 = sortrows(tmp',2);
nclusts = tmp2(:,1)';
rad = tmp2(:,2)';

if length(varargin) > 0
    hold on; 
    plot(rad,nclusts,varargin{1},'LineWidth',2,'LineStyle','--')
end

return
