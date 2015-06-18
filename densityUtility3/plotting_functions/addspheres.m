function sphere_handles = addspheres(surface_handles, DB, varargin)

% DB                    database struct
%                   - must have xyz, contrast, condf
% optional args:
% 'contrast_num'        vector of contrast numbers, for averaging coords
% 'average_nearby'      average nearby clusters from same contrast
% 'max_dist'            max distance from surface coordinate for plotting
%                       - could be a problem for interior points
%                       - may need to add clusters
% 'reflect_medial'      reflect medial points...
% 'color'               cell

addtext = 0;
reflect_medial = 0;
max_dist = 8;   % max dist in mm

shapes = {'sphere' 'sphere'}; % one for each type, unique vals in condf
colors = {[1 0 0] [1 0 1]};

% do the same within-contrast averaging as is done for other plots
% --------------------------------------------------------------------
[xyzall,codesall,orderall] = average_nearby_xyz(DB.xyz,12,DB.Contrast);
types = DB.condf(orderall);

connums = DB.Contrast(orderall); % for text labels only

[xyzall,codesall,orderall] = average_nearby_xyz(xyzall,8,codesall);
types = types(orderall);
connums = connums(orderall);

% --------------------------------------------------------------------

xyz_orig = xyzall;
labels_orig = codesall;

indx = 1;
sphere_handles = {};

for t = unique(types)'
    %     xyzt = xyz_orig(types==t,:);
    
    xyzt = [xyz_orig connums];
    xyzt = xyzt(types==t,:);
    
    if reflect_medial
        % keep only the medial ones, within 18 mm
        wh = abs(xyzt(:, 1)) < 18;
        xyzt(~wh, :) = [];
        
        % make all para-sagg for visibility
        xyzt(:, 1) = 6;
    end
    
    % for each surface
    for s = 1:length(surface_handles)
        
        d = dist_to_surface(xyzt(:,1:3), get(surface_handles(s), 'Vertices'));
        wh = d <= max_dist; % <- close enough
        
        xyzts = xyzt(wh, :); %<- s for in-surface
        
        if ~isempty(xyzts)
            switch shapes{t}
                case 'sphere'
                    sphere_handles{indx} = cluster_image_shape(xyzts(:,1:3), 'colors', colors{t}, 'radius', 4, 'sphere');
                case 'cube'
                    sphere_handles{indx} = cluster_image_shape(xyzts(:,1:3), 'colors', colors{t}, 'radius', 2, 'cube');
            end
            material dull
        end
        
        indx = indx + 1;
       
        if addtext
            % For text
            
            wh = [];
            for j = 1:size(xyzts,1)
                wh(j) = length(num2str(xyzt(j,4)))>1;
            end;
            wh = wh';
            
            xyzts(:,1) = xyzts(:,1) + 9;
            xyzts(wh==1,2) = xyzts(wh==1,2) + -2;
            
            for p = 1:size(xyzts,1)
                % text labels
                text(xyzts(p,1),xyzts(p,2),xyzts(p,3),...
                    num2str(xyzts(p,4)),'Color','w','FontSize',12,'FontWeight','bold');
            end
        end
        
    end % surfaces
    
end % types

end % function




% % select coords that are close to surface vertices5
%     % coords, cscale{1}, and alphascale{1} all have same indices
%     % ---------------------------------------------------------------------
%     cmax = max(smallv, [], 1);
%     cmin = min(smallv, [], 1);
%     fprintf('%3.0f coords.  selecting: ', size(coords, 1));
%
%     % omit wh and wh2 - outside scope of this coord set
%     wh = any(coords - repmat(cmax, size(coords, 1), 1) > mind, 2);
%     wh2 = any(repmat(cmin, size(coords, 1), 1) - coords > mind, 2);
%
%     % if cscale is matrix, must select these values of cscale as well!
%     if length(cscale) > 0 && size(cscale{1}, 1) == size(coords, 1)
%         cscale{1}(wh | wh2,:) = [];
%     end
%
%     if length(alphascale) > 0 && size(alphascale{1}, 1) == size(coords, 1)
%         alphascale{1}(wh | wh2,:) = [];
%     end
%
%     coords(wh | wh2,:) = [];
%
%     if isempty(coords), return, end
%
%     nc = size(coords, 1);
%     fprintf('%3.0f\n', nc);

%
%     function [vertex_indices, d] = find_in_radius(xyz2, setno, i, smallv, mind, whverts)
%     % output: indices of vertices in BIG list, and distances
%
%     % get vertices v within box -- fast method
%     wh = find(all(abs(bsxfun(@minus, xyz2{setno}(i,:), smallv)) <= mind, 2));
%     d = dist_tmp(smallv(wh,:), xyz2{setno}(i,:)');
%
%     % convert back to big list
%     vertex_indices = whverts(wh(d < mind));
%     end
%


%% may need thalamus, PAG only -- with contrast numbers
%
% DB.x = DB.xyz(:, 1);
% DB.y = DB.xyz(:, 2);
% DB.z = DB.xyz(:, 3);
%
% p = '/Users/tor/Documents/matlab_code_external/3DheadUtility/SPM2_brains';
% addpath(p);
% cl = mask2clusters(which('spm2_thal.img'));
%
% cl = database2clusters(DB, cl, 4);

% %% Thalamus
%
% create_figure('thal')
%
% condf = cl(1).condf(cl(1).wh_points);
% for i = 1:size(cl(1).XYZmm, 2)
%
%     hpatch = cluster_image_shape(cl(1).XYZmm(:, i)', 'colors', colors{condf(i)}, 'radius', 3, 'sphere');
%     material dull
% end
%
% h2 = addbrain('thalamus');
% set(h2, 'FaceColor', [.4 .4 .4]);
% camlight right
% axis image; axis tight
% view(133, 14);
% % h3 = addbrain('brainstem');
% % set(h3, 'FaceColor', [.5 .5 .5]);
% axis off
%
% %% Amy
%
% cl = mask2clusters(which('spm2_amy.img'));
%
% cl = database2clusters(DB, cl, 4);
%
% create_figure('amy')
%
% condf = cl(1).condf(cl(1).wh_points);
% for i = 1:size(cl(1).XYZmm, 2)
%
%     hpatch = cluster_image_shape(cl(1).XYZmm(:, i)', 'colors', colors{condf(i)}, 'radius', 3, 'sphere');
%     material dull
% end
%
% % add other side - reflect to right
% condf = cl(2).condf(cl(2).wh_points);
% cl(2).XYZmm(1, :) = -cl(2).XYZmm(1, :);
% for i = 1:size(cl(2).XYZmm, 2)
%
%     hpatch = cluster_image_shape(cl(2).XYZmm(:, i)', 'colors', colors{condf(i)}, 'radius', 3, 'sphere');
%     material dull
% end
%
% h2 = addbrain('amygdala');
% set(h2, 'FaceColor', [.4 .4 .4]);
% camlight right
% axis image; axis tight
% view(133, 14);

