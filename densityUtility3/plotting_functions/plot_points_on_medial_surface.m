function plot_points_on_medial_surface(coords,colors,varargin)
    % plot_points_on_medial_surface(coords,colors,varargin)
    % plot_points_on_medial_surface(coords,colors,[factor variable],[factor levels],[newfig],[mytextlabels])
    %
    % examples:
    % plot_points_on_medial_surface([EMDB.x EMDB.y EMDB.z],{'ro'});
    %
    % plot_points_on_medial_surface([EMDB.x EMDB.y EMDB.z],{'bs' 'yo'},EMDB.valence,{'pos' 'neg'});
    %
    % plot_points_on_medial_surface([INHDB.x INHDB.y INHDB.z], ...
    % {'gs' 'ro' 'b^' 'mo' 'yv' 'cs' 'c^'},INHDB.Task);
    %
    %LTMDB.Memtype(find(strcmp(LTMDB.Memtype,'Autobiog'))) = {'Autobiographical'};
    %LTMDB.Memtype(find(strcmp(LTMDB.Memtype,'encoding'))) = {'Encoding'};
    %LTMDB.y(abs(LTMDB.y) > 100) = NaN;
    % plot_points_on_medial_surface([LTMDB.x LTMDB.y LTMDB.z], ...
    %{'ro' 'ys' 'bv' 'cs' 'mo' 'g^'},{'Autobiographical' 'Encoding' 'Recall' 'Recognition' 'Retrieval' 'Source'});
    %
    % plot_points_on_medial_surface([PAINDB3.x PAINDB3.y PAINDB3.z] ...
    %,{'go' 'rs' 'bv' 'cs' 'yo' 'g^' 'mo'},{'Cutaneous' 'Electric_shock' 'Intramuscular' 'Laser' 'Thermal' 'Vestibular' 'Visceral'});
    %
    % plot_points_on_medial_surface([PAINDB3.x PAINDB3.y PAINDB3.z],{'ro' 'bs'},PAINDB3.Right_vs_Left,{'Left' 'Right'});
    %
    % plot_points_on_medial_surface(XYZ,{'bo' 'go' 'ro' 'rd' 'rs'},color, ...
    %{'bo' 'go' 'ro' 'rd' 'rs'},1,letter);
    % plot_points_on_medial_surface(XYZ,{'bo' 'go' 'gs' 'r^' 'ro' 'rd'
    % 'rs'},color,{'bo' 'go' 'gs' 'r^' 'ro' 'rd' 'rs'},1,letter);

    % see also plot_points_on_surface, plot_points_on_surface2 (best)
    % tor wager
    newmethod = 0;  % needs to be developed

    surfacecutoff = 16;                 % plot points this distance or less from surface, in mm.
    mymarkersize = 8;
    plottext = 0;
    mytext = [];    % would be cell array of text strings
    colordefi = [];
    levelnames = [];

    newfig = 1;


    if ~exist('coords'), coords = [0 0 0]; end
    if ~exist('colors'), colors = {'ko'}; end

    if length(varargin)>0, colordefi = varargin{1}; end
    if length(varargin)>1, levelnames = varargin{2}; end
    if length(varargin)>2, newfig = varargin{3}; end
    if length(varargin)>3, plottext = 1; mytext = varargin{4}; end

    %if length(colors) == size(coords,1)

    [colorout,coords,colornames,whomit]   = get_colors(colors,colordefi,levelnames,coords);


    % add 5 to Z value to adjust origin
    disp('WARNING!  ADDING 5 TO Z VALUES TO MATCH WHAT I THINK IS CAUSED BY SHIFT IN ORIGIN.')
    coords(:,3) = coords(:,3) + 6;


    coords(:,2) = -coords(:,2);
    coordscopy = coords;

    %if length(colors) < size(coords,1), colors = repmat(colors,size(coords,1),1); end

    [D,hdr,origin] = read_image_data;


    % average coords together, if study-unique text labels are entered
    if ~isempty(mytext)
        mytext(find(whomit)) = [];

        fprintf(1,'Averaging nearby coords at 12 mm ... starting with %3.0f coords. ',size(coords,1));
        % average nearby coordinates together!  12 mm
        [coords,mytext,order] = average_nearby_xyz(coords,12,mytext);
        colorout = colorout(order);
        coordscopy = coordscopy(order,:);

        % second pass, 8 mm
        [coords,mytext,order] = average_nearby_xyz(coords,8,mytext);
        colorout = colorout(order);
        coordscopy = coordscopy(order,:);

        fprintf(1,'finished with %3.0f.\n ',size(coords,1));
    else
        for i = 1:size(coords,1), mytext(i,1) = 'x'; end
    end


    % ===================================================================
    % * Left Medial *
    % ===================================================================
    disp('Plotting left medial coordinates.');drawnow

    coords = coordscopy;


    disp('  Adjusting coordiates from mm to brain voxels');drawnow
    % coordinates entered in mm should be scaled to voxels

    coords(:,1) = round(coords(:,1) / hdr.xsize);
    coords(:,2) = round(coords(:,2) / hdr.ysize);
    coords(:,3) = round(coords(:,3) / hdr.zsize);

    % express all coordinates in voxel space, converting from distance from the origin
    for i = 1:size(coords,1)
        coords(i,:) = coords(i,:) + origin;
    end


    % --- make figure ----------------------------------------------
    plothandle = create_figure('Medial point plot', 1, 2, ~newfig);
    
% %     %if newfig, plothandle = tor_fig(1,2); ,end
subplot(1,2,1); plothandle = gca; subplot(1,2,2); plothandle(2) = gca;
subplot(1,2,1);
    plotorigin(origin,D)

    % ---- select based on laterality ------------------------------
    wh = coordscopy(:,1) >= -(surfacecutoff) & coordscopy(:,1) <= 0;
    pcoords = coords(wh,:);
    pcolors = colorout(wh,:);		% svert is coordinates to plot


    plot_points(pcoords,pcolors,colors,plottext,mytext(wh,:),mymarkersize)

    % ---- plot medial surface ----
    if newmethod && newfig
        addbrain('left');
    else
        E = D;
        if newfig
            D(:,1:48,:) = [];
            D(:,:,77:end) = [];
            make_surface(D)
            D = E;
        end
    end

    % ===================================================================
    % * Right Medial *
    % ===================================================================
    disp('Plotting right medial coordinates.');drawnow

    disp('Right Medial')
    coords = coordscopy;

    origin = hdr.origin(1:3,1)';
    D = E;
    % adjust origin because you removed voxels
    %origin(1) = origin(1) + -48;
    % origin is spec from top, matlab plots from bottom.  reverse y
    origin(2) = size(E,1) - origin(2);

    % coordinates entered in mm should be scaled to voxels
    coords(:,1) = round(coords(:,1) / hdr.xsize);
    coords(:,2) = round(coords(:,2) / hdr.ysize);
    coords(:,3) = round(coords(:,3) / hdr.zsize);

    disp('  Adjusting coordiates from mm to brain voxels');drawnow
    % adjust coords because you removed voxels
    %    coords(:,1) = coords(:,1) + -48;
    % express all coordinates in voxel space, converting from distance from the origin
    for i = 1:size(coords,1)
        coords(i,:) = coords(i,:) + origin;
    end


    % --- make figure ----------------------------------------------
    if ~isempty(plothandle), axes(plothandle(2))
    else tor_fig;
    end
    hold on

    origin
    % plot origin
    plot3([origin(1)+5 origin(1)+5],[origin(2) origin(2)],[1 size(D,3)],'k')
    plot3([origin(1)+5 origin(1)+5],[1 size(D,1)],[origin(3) origin(3)],'k')
    hold on; drawnow

    % ---- select based on laterality ------------------------------
    wh = coordscopy(:,1) <= surfacecutoff & coordscopy(:,1) >= 0;
    pcoords = coords(wh,:);
    pcolors = colorout(wh,:);

    disp(['plotting ' num2str(size(pcoords,1)) ' points.'])

    legh = plot_points(pcoords,pcolors,colors,plottext,mytext(wh,:),mymarkersize,80);
    hold on;



    % ---- plot medial surface ----

    if newmethod && newfig
        addbrain('left');
    elseif newfig
        D = E;
        D(:,45:end,:) = [];
        D(:,:,77:end) = [];
        disp('  Plotting surface');drawnow
        p1 = patch(isosurface(D, 50),'FaceColor',[1,.75,.65], ...
            'EdgeColor','none');
        p2 = patch(isocaps(D, 50),'FaceColor','interp', ...
            'EdgeColor','none');
        %view(75,10);
        view(90,5)
        axis tight; axis image; axis off
        colormap(gray(100));camzoom(1.4);
        camlight; camlight right; lighting gouraud
        isonormals(D,p1)
        drawnow
    end

    h = findobj('Type','light'); delete(h)
    subplot(1,2,1); camlight left; subplot(1,2,2); camlight left;
    subplot(1,2,1); camlight right; subplot(1,2,2); camlight right;

    if length(varargin) > 0
        legend(legh,colornames)
    end

%     subplot(1,2,1)
%     camzoom(.8)
%     subplot(1,2,2)
%     camzoom(1.1)

    return





function [D,hdr,origin] = read_image_data

    basename = 'scalped_single_subj_T1';
    [array,hdr] = readim2(basename);

    clear D
    for i = 1:size(array,3)
        E(:,:,i) = rot90(array(:,:,i));
    end

    origin = hdr.origin(1:3,1)';
    D = E;
    % adjust origin because you removed voxels
    origin(1) = origin(1) - 48;
    % origin is spec from top, matlab plots from bottom.  reverse y
    origin(2) = size(E,1) - origin(2);

    return





function [colorout,coords,colornames,whomit]   = get_colors(colors,colordef,levelnames,coords)
    % process colors, if a varargin string is entered
    % make new colors string with one color per point, defined by colordef
    % cell string to define colors = colordef

    colornames = [];

    if isempty(levelnames), levelnames = unique(colordef);  end

    if isempty(colordef)
        colorout = repmat(colors,size(coords,1),1)';
    else
        [indic,colornames] = string2indicator(colordef,levelnames);

        for i =1:size(indic,2)
            colorout(find(indic(:,i))) = colors(i);
        end

    end

    % get rid of nonselected levels
    whomit = [];
    for i = 1:length(colorout)
        if isempty(colorout{i}), whomit(i) = 1;  else, whomit(i) = 0;  end
    end
    colorout(find(whomit)) = [];
    coords(find(whomit),:) = [];

    colorout = colorout';
    return




function plotorigin(origin,D)
    %origin
    % plot origin
    plot3([origin(1)-5 origin(1)-5],[origin(2) origin(2)],[1 size(D,3)],'k')
    plot3([origin(1)-5 origin(1)-5],[1 size(D,2)],[origin(3) origin(3)],'k')
    hold on; drawnow
    return






function legh = plot_points(pcoords,pcolors,colornames,plottext,mytext,mymarkersize,varargin)

    legh = [];
    
    if length(varargin), add2x = varargin{1}; else, add2x = 0; end

    hold on
    disp(['plotting ' num2str(size(pcoords,1)) ' points.'])

    if plottext
        legh = plot3(0,0,0,colornames{1},'MarkerFaceColor',colornames{1}(1),'MarkerSize',mymarkersize);
        set(legh,'Visible','off');

        for i = 1:length(mytext)
            text(0+add2x,pcoords(i,2),pcoords(i,3),...
                mytext{i},'Color',pcolors{i}(1),'FontSize',12,'FontWeight','bold');
        end
    else

        for i = 1:length(colornames)
            % ---- select based on color group -------------------------
            plotcolors = pcolors(strcmp(pcolors,colornames{i}) == 1,:);
            svert = pcoords(strcmp(pcolors,colornames{i}) == 1,:);

            % ---- plot the points -------------------------------------
            hold on

            if ~isempty(svert)
                legh(i) = plot3(zeros(size(svert,1),1)+add2x,svert(:,2),svert(:,3),colornames{i},'MarkerSize',mymarkersize,'MarkerFaceColor',colornames{i}(1));
            end
        end
    end
    drawnow
    return





function make_surface(D)
    disp('  Plotting surface');drawnow
    p1 = patch(isosurface(D, 50),'FaceColor',[1,.75,.65],...,...
        'EdgeColor','none');
    p2 = patch(isocaps(D, 50),'FaceColor','interp',...
        'EdgeColor','none');
    view(269,5); axis tight; axis image; axis off
    colormap(gray(100)); camzoom(1.35);
    camlight left; camlight; camlight left; lighting gouraud
    isonormals(D,p1)
    drawnow
    hold on
    return

