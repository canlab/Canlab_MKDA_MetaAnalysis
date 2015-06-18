function [parcel_stats, all_surf_handles] = Meta_Parcel_results(parcels, parcel_stats)
% [parcel_stats, all_surf_handles] = Meta_Parcel_results(parcels, parcel_stats)
%
% Get results and images for parcels created in Meta_Parcel

    nparcels = length(parcels);
    fprintf('There are %3.0f parcels\n', nparcels);

    for i = 1:nparcels
        colors{i} = rand(1, 3);
    end

    % orthviews
    % ------------------
    cluster_orthviews(parcels(1), colors(1), 'solid');
    for i = 2:nparcels
        cluster_orthviews(parcels(i), colors(i), 'solid', 'add');
    end

    % NMDS map
    % ------------------
    % scale correlations (or beta weights) so that they are a similarity matrix between 0 and 1.
    parcel_stats.D = (1 - parcel_stats.taumtx) ./ 2;
    [parcel_stats.GroupSpace,parcel_stats.obs,parcel_stats.implied_dissim] = shepardplot(parcel_stats.D,[]);

    create_figure('nmdsfig');
    f1 = nmdsfig(parcel_stats.GroupSpace,'classes',(1:nparcels)','sig',parcel_stats.sig,'legend',{'Pos' 'Neg'}, 'sizes', 18, 'colors', colors);
    axis image
    axis equal
    axis off

    % surfaces
    % ------------------
    all_surf_handles = mediation_brain_surface_figs([], []);

    for i = 1:nparcels

        cluster_surf(parcels(i), 3, colors(i), all_surf_handles);

    end

end