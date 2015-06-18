function xyzvox = meta_get_voxel_coords(DB,contrastnumber,V,maskdata)
% Given a) x,y,z fields in DB, 2) 3-D mask data; and 3) a contrast number indexed in DB.Contrast,
% return voxel coordinates xyzvox of in-mask coordinate points.
%
% xyzvox = meta_get_voxel_coords(DB,contrastnumber,V,maskdata)
% See meta_prob_activation (parent function)
%
% Tor Wager, 6/30/06

% get study coordinates for contrast
xyzstudy = [DB.x DB.y DB.z]'; xyzstudy = xyzstudy(:,find(DB.Contrast==contrastnumber)); % 3-rows
xyzvox = mm2voxel(xyzstudy,V,1);   % input: 3 rows, output: 3 columns

sz = size(maskdata);

% limit to in-mask voxel coords
for i = 1:size(xyzvox,1)
    
    if any(xyzvox(i,:) > sz) || any(xyzvox(i,:) < 1)
        fprintf(1,'\nWarning! Coords %3.0f %3.0f %3.0f in Contrast %3.0f are outside mask.  Omitting this coord\n     ', ...
            xyzstudy(1,i),xyzstudy(2,i),xyzstudy(3,i),contrastnumber);
    else
        inmask(i,1) = maskdata(xyzvox(i,1),xyzvox(i,2),xyzvox(i,3));
    end
    
end

if ~exist('inmask', 'var')
    disp('Contrast has no valid in-mask coordinates.');
    xyzvox = [];
    
else
    xyzvox = xyzvox(find(inmask),:);
end

return