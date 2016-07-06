function [axishan, pointhan, surfhan] = plot_points_near_regions_on_surface(DB, r, wh_points, color)
% [axishan, pointhan, surfhan] = plot_points_near_regions_on_surface(DB, r, wh_points, color)
%
% DB: Database of meta-analysis coordinates, from Meta_SETUP.m
% r:  region object with clusters of regions
% wh_points: 1/0 logical vector of which coordinates in DB.xyz to include. 
%           Should be length of DB.x
% color: rgb triplet for color


DB.xyz = [DB.x DB.y DB.z];

contrastnums = DB.Contrast(wh_points);
xyz = DB.xyz(wh_points, :);
condf = ones(size(xyz, 1), 1);

% Select coordinates near regions - for display
% (Eliminate background noise)
% ----------------------------------------------
[xyz, indx] = select_coordinates_near_regions(r, xyz, 10);
condf = condf(indx);
contrastnums = contrastnums(indx);
%descrip = descrip(indx);
study = DB.study(indx);

% average nearby coordinates within 12 mm that come from the same contrast
% ----------------------------------------------
[xyz,codesall, indx] = average_nearby_xyz(xyz, 12, contrastnums);
condf = condf(indx);
contrastnums = contrastnums(indx);
%descrip = descrip(indx);
study = study(indx);


%%
[axishan, pointhan, surfhan] = plot_points_on_surface2(xyz, color, condf);
material(cat(2, pointhan{:}), 'dull');
subplot(3, 2, 6);
% camzoom(1.1)
% camzoom(1.05)
for i = 1:length(pointhan)
    set(pointhan{i}, 'SpecularColorReflectance', .7)
end

axes(axishan(1))
lightRestoreSingle;
axes(axishan(6))
camzoom(1.1);

drawnow
snapnow

end % function