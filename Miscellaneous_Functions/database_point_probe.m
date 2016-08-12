function newcl = database_point_probe(DB,XYZlist,r,varargin)
% newcl = database_point_probe(DB,XYZlist,r,varargin)
%
% given a point list (n x 3 mm coordinates) and a radius r,
% gets a clusters structure with points within r mm of any input coordinate in the
% list
% prints a table of unique contrasts within studies that have points 
% within the specified regions.
%
% USED IN: 
% meta_analysis_table.m
%
% tor wager
%
% Example:
% 
% fl = {'Right_vs_Left' 'Right' 'Left'};
% cl = database_point_probe(PAINDB2,[-42 0 2],10,fl);
% dbcluster_point_table(cl)
%

% OLD functionality:
% given a point list (n x 3 mm coordinates) and a radius r,
% gets a clusters structure with points within r mm of point list
% each center is a cluster.
% prints a table of unique contrasts within studies that have points 
% within the specified regions.
%
% OLD (deprecated/commented out in code):
% An optional argument specifies the field list to define 'unique
% contrasts' and print results.  only for "Long" output format - hard-coded
% in this script.  If you want to change to Long, edit the script

if size(XYZlist, 2) ~= 3
    XYZlist = XYZlist';
end

if size(XYZlist, 2) ~= 3
    error('XYZlist must be n x 3 matrix of x, y, and z coordinates.')
end

% define clusters with centers
cl = [];
cl.XYZmm = XYZlist';  % transpose just for consistency with clusters format
cl.mm_center = nanmean(XYZlist, 1);
% OLD -- center-by-center list
% for i = 1:size(XYZlist,1)
%     cl(i).mm_center = XYZlist(i,:);
%     
% end


% Get Points in database

if isfield(DB,'dbname'), 
    fprintf(1,'\n%s\n',DB.dbname),
    DB = rmfield(DB,'dbname');
end

XYZ = [DB.x DB.y DB.z];
N = fieldnames(DB);
newcl = cl;

if ~isfield(DB, 'x'), error('Must have x, y, z fields with coordinates.'), end
npoints = length(DB.x);

% ---------------------------------------------
% Put variables in structure
% ---------------------------------------------

for i = 1:length(newcl)
    
    % find which points are in-cluster (w/i 3.464 mm)
    % max euclidean distance in same 2 x 2 x 2 mm voxel is sqrt(4+4+4) = 3.46 mm
    % use dominance metric: max dist on any dimension is 2 mm
    
    %disp(['Doing cluster ' num2str(i)])

    % list of points in XYZ database that are within r mm of cluster center
    %[xyz2,wh,d] = points_in_sphere(cl(i).mm_center,r,XYZ);
    
    % list of points in XYZ database that are within r mm of any coords in
    % input list
    [xyz2, wh, d] = points_in_sphere(cl(i).XYZmm', r, XYZ, 'any');
    
    newcl(i).distance = d;
    
    if isempty(wh)
        
    else
        for ind = wh'

            for k = 1:length(N)

                mydat = DB.(N{k});
                if length(size(mydat)) == 2 && size(mydat,1) ~= npoints, mydat = mydat'; end

                if size(mydat,1) == npoints
                    % define new field, if necessary
                    
                    if ind == wh(1)

                        newcl(i).(N{k}) = mydat(ind,:);


                        %eval(['newcl(i).' N{k} ' = DB.' N{k} '(ind,:);'])
                    else

                        newcl(i).(N{k})(end+1, :) = mydat(ind,:);

                        %eval(['newcl(i).' N{k} '(end+1,:) = DB.' N{k} '(ind,:);'])
                    end
                end

            end
        end

        newcl(i).XYZ = XYZ(wh,:)';
        %disp(['Found ' num2str(ind-1) ' points'])
    
    end     % empty
    
end

% ---------------------------------------------
% Print table
% ---------------------------------------------

% if length(varargin) > 0
%     fl = varargin{1};
% end
% 

% LONG format
%OUT = dbcluster_contrast_table(newcl,DB,fl);

% SHORT format
dbcluster_point_table(newcl);


return



