function [con,NLCON] = cluster_indiview(cl,EXPT)
% function [con,NLCON] = cluster_indiview(cl,EXPT or con img list or con image data) 
%
% tor wager


% used in button-up fcn callback

global VOL
global con
global f
global f2
global NLCON

    % -------------------------------------------------------------------
    % * essential stuff for the viewing
    % -------------------------------------------------------------------
    
% view clusters
cluster_orthviews(cl,'unique');
set(gcf,'WindowButtonUpFcn','indiview_get_data')


% get coordinate mapping matrix
VOL = struct('M',cl(1).M);

% prepare figure
f1 = figure('Color','w','Name','Contrast value view');
f = f1;     % button callback uses figure f


% prepare contrast images
if isstruct(EXPT), 
    % it's EXPT
    con = EXPT.SNPM.P;,
    
    % H, T, W stuff
    if isfield(EXPT,'NLCON')
        NLCON = EXPT.NLCON;
        if isfield(EXPT,'DX'), if isfield(EXPT.DX,'dxnames'),  NLCON.names = EXPT.DX.dxnames;,end,end
        f2 = figure('Color','w','Name','HTW view');
    end
else, 
    % it's data
    con = EXPT; clear EXPT;, 
end




    % -------------------------------------------------------------------
    % * some names and stuff for interpretation
    % -------------------------------------------------------------------
    
% show name of contrast, if we can
if isfield(cl,'P'), disp('CONTRAST: '); disp(cl(1).P(1,:));, disp('----------------------------'); end



% list names of images, if we can
if isstr(con{1}),
    for i = 1:length(con)
        fprintf(1,'%3.0f\t%s\n',i,con{i}(1,:));
    end
end





    % -------------------------------------------------------------------
    % * call for the first time (internal function)
    % -------------------------------------------------------------------

coord = spm_orthviews('Pos');
coord = mm2voxel(coord',struct('M',cl(1).M));
mm = round(spm_orthviews('Pos')');

%[coord,con,mm] = indiview_get_data(coord,con,mm,f1);

input('Press a key to exit.')



return







%function [coord,con,mm] = indiview_get_data(coord,con,mm,f1)
% the stand-alone version is used in button callback, and uses global vars

    % -------------------------------------------------------------------
    % * load images, if necessary, and get data
    % -------------------------------------------------------------------

    if isstr(con{1}), disp('Loading images (this may take a few minutes).');,end
   
    % if we already have volumes, P is image data; otherwise, names, and
    % vols are loaded
   
    for i = 1:length(con)
        [ts(i),con{i}] = timeseries3(coord,con{i});
        
        dat(:,i) = ts(i).indiv;
    end
 
 
    figure(f1); barplot_columns(dat);
    
    title(['mm: ' num2str(mm) ' vox: ' num2str(coord)],'FontSize',24);
    
    
%return