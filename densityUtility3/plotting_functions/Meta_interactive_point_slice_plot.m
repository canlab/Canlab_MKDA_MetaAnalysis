function Meta_interactive_point_slice_plot(MC_Setup, DB)
    % Meta_interactive_point_slice_plot(MC_Setup, DB)
    %
    % Set up interactive point plotting on slice
    %
    % tor wager, nov 2007
%
% Simple example for checking points manually, rather than setting up the interactive plotter:
% load SETUP
% V = DB.maskV;
% DB.xyz = [DB.x DB.y DB.z];
% 
% dd = DB.direction(strcmp(DB.disorder, 'SP'), :);
% coords = DB.xyz(strcmp(DB.disorder, 'SP'), :);
% 
% pos = spm_orthviews('Pos')
% xyzvox = mm2voxel(pos, V.mat);
% z = round(xyzvox(3));
% x = round(xyzvox(1));
% create_figure('points'); plot_points_on_slice(coords(strcmp(dd, 'ptmore'), :), 'slice', z, 'close_enough', 10, 'markerfacecolor', [1 0 0]);
% plot_points_on_slice(coords(strcmp(dd, 'ctrmore'), :), 'slice', z, 'close_enough', 10, 'markerfacecolor', [0 0 1], 'nodraw');


    disp('Setting up Meta point plot')


    conditions = MC_Setup.Xinms;
    colors = {'b' 'r' 'g' 'y' 'm' 'k'};
    maxdistance = DB.radius_mm;

    wh = find(MC_Setup.connames{1} == '_'); 
    if ~isempty(wh)
        wh = wh(1);
        fieldname = MC_Setup.connames{1}(1:wh-1);
    else
         fieldname = MC_Setup.connames{1};
    end
    
    myxyz = [DB.x DB.y DB.z];
    V = DB.maskV;

    disp('Setting graphics callback')

    callback_handle = @(str, pos, reg, hReg) point_plot_callback_wrapper(str, pos, reg, hReg);

    hSpmFig = spm_figure('GetWin', 'Graphics');

    hReg = uicontrol(hSpmFig, 'Style', 'Text', 'String', 'InteractiveViewer hReg', ...
        'Position', [100 200 100 025], 'Visible', 'Off', ...
        'FontName', 'Times', 'FontSize', 14, 'FontWeight', 'Bold', ...
        'HorizontalAlignment', 'Center');
    hReg = spm_XYZreg('InitReg', hReg, V.mat, V.dim(1:3)');
    spm_XYZreg('Add2Reg', hReg, 0, callback_handle);
    spm_orthviews('Register', hReg);


    disp('Ready!')


    % inline

    function point_plot_callback_wrapper(str, pos, reg, hReg)

        switch str
            case 'SetCoords'
                plot_callback(pos, DB, fieldname, myxyz, conditions, colors, maxdistance, V)

            otherwise
                disp('Unknown callback command from spm_XYZreg');
        end

    end

end  % main function



function plot_callback(pos, DB, fieldname, myxyz, conditions, colors, maxdistance, V)


    xyzvox = mm2voxel(pos, V.mat);
    z = round(xyzvox(3));
    x = round(xyzvox(1));
    
    create_figure('Slice view', 1, 2);

    subplot(1, 2, 1);

    for i = 1:length(conditions)

        xyz_to_plot = myxyz(strcmp(DB.(fieldname), conditions{i}), :);
        
        [handles, wh_slice, my_z] = plot_points_on_slice(xyz_to_plot, 'slice', z, 'markerfacecolor', colors{i}, 'marker', 'o', 'close_enough', maxdistance);

        left_right_counts(i, 1) = sum(xyz_to_plot(:,1) < 0 & xyz_to_plot(:, 3) >= pos(3) - maxdistance & xyz_to_plot(:, 3) <= pos(3) + maxdistance);
        left_right_counts(i, 2) = sum(xyz_to_plot(:,1) > 0 & xyz_to_plot(:, 3) >= pos(3) - maxdistance & xyz_to_plot(:, 3) <= pos(3) + maxdistance);
        
    end


    subplot(1, 2, 2);

    for i = 1:length(conditions)

        [handles, wh_slice, my_x] = plot_points_on_slice(myxyz(strcmp(DB.(fieldname), conditions{i}), :), 'slice', x, 'markerfacecolor', colors{i}, 'marker', 'o', 'close_enough', maxdistance, 'sagg');


    end
    

    print_matrix(left_right_counts, {'Left' 'Right'}, conditions);
    assignin('base', 'left_right_counts', left_right_counts);
  
end  