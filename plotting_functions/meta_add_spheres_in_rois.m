function meta_add_spheres_in_rois(DB, varargin)
% Add spheres to selected ROIs
%
% [p, mesh_struct] = brainstem_slices_3d;
% meta_add_spheres_in_rois(PLOTINFO{1}, 'brainstem');
% colormap gray
%
% DB must be a DB or PLOTINFO struct with these fields:
% - colors
% - xyz
% - condf
%
% (for table)
% - nums
% - descrip
% - study
%
% (to save in diary)
% - diaryname

clear cl

colors = DB.colors;

indx = 0;

for i = 1:length(varargin)
    switch varargin{i}
        
        case 'brainstem'
            indx = indx + 1;
            cl{indx} = mask2clusters(which('spm8_brainstem.img'));
            cl{indx} = database2clusters(DB, cl{indx}, 5);
            cl{indx}(1).title = varargin{i};
            
        case 'amygdala'
            indx = indx + 1;
            cl{indx} = mask2clusters(which('spm2_amy.img'));
            cl{indx} = database2clusters(DB, cl{indx}, 5);
            cl{indx}(1).title = varargin{i};
            
        case 'hippocampus'
            indx = indx + 1;
            cl{indx} = mask2clusters(which('spm2_hipp.img'));
            cl{indx} = database2clusters(DB, cl{indx}, 5);
            cl{indx}(1).title = varargin{i};
            
        otherwise
            % not entered yet
    end
end

% Compile XYZ list
% -----------------------------------------------------------------

% Print diary if required field entered
if isfield(DB, 'diaryname') && ~isempty(DB.diaryname)
    diary(DB.diaryname)
end

    
clear xyzcut condfcut
for i = 1:length(cl)

    fprintf('\n---------------------------------\n%s\n---------------------------------\n', cl{i}(1).title);
    
    print_table(cl{i}, DB);
    
    fprintf('\n');
    
    for j = 1:length(cl{i})
        %whin = cat(1, cl{i}.wh_points);
        whin = cl{i}(j).wh_points;
        
        xyzcut{i}{j} = cl{i}(j).xyz(whin, :);
        condfcut{i}{j} = cl{i}(j).condf(whin);
        
    end
    
    xyzcut{i} = cat(1, xyzcut{i}{:});
    condfcut{i} = cat(1, condfcut{i}{:});
    
end

xyzcutaway = cat(1, xyzcut{:});
condfcutaway = cat(1, condfcut{:});

if isfield(DB, 'diaryname') && ~isempty(DB.diaryname)
    diary off
end

% Plot
% -----------------------------------------------------------------

u = unique(condfcutaway);
u(u == 0) = [];             % remove out-of-scope points

for clas = 1:length(u)
    pthan{clas} = cluster_image_sphere(xyzcutaway(condfcutaway==u(clas), 1:3), 'color', colors{u(clas)}, 'radius', 3);
end

for i = 1:length(pthan)
    set(pthan{i}, 'SpecularColorReflectance', .7)
end

lighting gouraud  % re-set for points


end % function



% Table
function print_table(cl, PLOTINFO)
for i = 1:length(cl)
    
    mynums = cl(i).nums(cl(i).wh_points);
    mydescrip = cl(i).descrip(cl(i).wh_points);
    mystudy = cl(i).study(cl(i).wh_points);
    
    mycolors = PLOTINFO.colors(cl(i).condf(cl(i).wh_points));
    
    for j = 1:length(mynums)
        if ~ischar(mycolors{j}), mycolors{j} = num2str(mycolors{j}); end
        
        fprintf('%3.0f\t%s\t%s\t%s\n', mynums(j), mycolors{j}, mystudy{j}, mydescrip{j});
    end
    
end
end % function