function [newcl,DB] = database2clusters(db_mat,clusters,mydist)
% function [newcl,DB] = database2clusters(db_mat,clusters,mydist)
%
% tor wager
%
% This function turns a database file (stored as variables in a .mat file)
% as created with read_database.m
% into a clusters structure, by finding the points
% that fall within mydist mm from any in-cluster voxel
%
% db_mat: the .mat file containing the database
% clusters: a clusters structure, see tor_extract_rois.m, that has the
%           clusters to group data points in the database at a specified distance
% mydist: distance in mm (euclidean)
%   * if empty, uses 2 mm in a dominance metric
%       i.e., point must be w/i 2 mm of cluster voxel on every
%       dimension (x,y,z)
%
% This function also adds an 'alltask' field to each cluster whose entries
% are the name of the database, minus extension.  This is useful for
% comparing databases with dbcluster_meta_meta.m
%
% Having a field called 'Z', or anything already in clusters, will cause
% this to crash.
%
% newcl.XYZ contains point coordinates
%newcl = database2clusters('switching_database.mat',clusters,5);
%montage_clusters([],clusters,{'r' 'b'},XYZ)
%montage_clusters([],clusters,{'r' 'b'},cat(1,newcl.XYZ))
%
% Revised aug 06 by tor wager
% Now takes input structure as well as a mat filename
% Example: image points in brainstem
% --------------------------------------------
% pname = which('spm2_brainstem.img')
%    cl = mask2clusters(pname)
%    cl = database2clusters(DB,cl,2);
% figure; plot_points_on_brain(cl(1).xyz,{'go'},[],0);
% addbrain('brainstem')
% axis image; axis off; axis vis3d; view(90,10)
%
% Example: from saved SETUP.mat in densityUtility3
% load SETUP DB
% img = 'Activation_FWE_all.img';
% cl = mask2clusters(img);
% [cl] = database2clusters(DB,cl, 10);

% ---------------------------------------------
% Put variables in structure
% ---------------------------------------------
if isstruct(db_mat)
    DB = db_mat;

else
    load(db_mat)

    [tmp1,vname] = fileparts(db_mat);   % save name of calling string in vname

    if ~(exist('Study')==1) && exist('study')==1, Study = study; end

    if (exist('Study')==1), num = length(Study); else num = length(x); end
    N = whos;

    for i = 1:length(N)

        if any(N(i).size == num)
            %eval(['whos ' N(i).name])
            str = (['if size(' N(i).name ',2) > size(' N(i).name ',1), ' N(i).name '=' N(i).name ''';, end']);
            eval(str)
            %eval(['whos ' N(i).name])
            eval(['DB.' N(i).name ' = ' N(i).name ';'])
        end

    end

    DB.alltask(1:size(DB.x,1),1) = {vname};

end

if isfield(DB, 'xyz')
    XYZ = DB.xyz;
else
    XYZ = [DB.x DB.y DB.z];
end

if ~isfield(DB, 'x') || isempty(DB.x)
    DB.x = DB.xyz(:, 1);
    DB.y = DB.xyz(:, 2);
    DB.z = DB.xyz(:, 3);
    
end

N = fieldnames(DB);
newcl = clusters;

dodom = 0;
if isempty(mydist), mydist = 2; dodom = 1; end

% ---------------------------------------------
% Put variables in structure
% ---------------------------------------------

for i = 1:length(clusters)

    % find which points are in-cluster (w/i 3.464 mm)
    % max euclidean distance in same 2 x 2 x 2 mm voxel is sqrt(4+4+4) = 3.46 mm
    % use dominance metric: max dist on any dimension is 2 mm

    disp(['Doing cluster ' num2str(i)])
    
    XYZc = clusters(i).XYZmm';
    ind = 1;
    newcl(i).XYZ = [];
    
    % index of values within r mm of cluster
    indx = iimg_xyz2spheres(XYZc,XYZ,mydist);
    
    wh = find(indx);
    newcl(i).wh_points = wh;
    
    for k = 1:length(N)
        
        eval(['myvar = DB.' N{k} ';']);
        if size(myvar,1) == size(XYZ,1)
            
            %             if isfield(newcl(i), N{k}) && ~strcmp(class(newcl(i).(N{k})), class(DB.(N{k}))) %isfield(newcl(i), N{k}) && ~isempty(newcl(i))
            %                 fprintf('Field %s has mismatched type in DB and cl. skipping.\n', N{k});
            %             else
            try
                
                eval(['newcl(i).' N{k} '(wh,:) = DB.' N{k} '(wh,:);'])
                
            catch
                fprintf('Error with %s : mismatched type/special meaning? skipping.\n', N{k});
            end
            
            
        end
    end

    newcl(i).XYZ = XYZ(wh,:)';
    disp(['Found ' num2str(length(wh)) ' points'])

    % % %
    % % %         for j = 1:size(XYZ,1)   % find all XYZ within cluster
    % % %
    % % %             xm = abs(XYZ(j,1) - XYZc(:,1));
    % % %             ym = abs(XYZ(j,2) - XYZc(:,2));
    % % %             zm = abs(XYZ(j,3) - XYZc(:,3));
    % % %
    % % %             if dodom
    % % %                 % dominance metric
    % % %                 d = max([xm ym zm]');
    % % %             else
    % % %                 % euclidean dist
    % % %                 d = (xm.^2 + ym .^2 + zm .^2) .^.5;
    % % %             end
    % % %
    % % %
    % % %             %mincl = XYZc(find(d==min(d)),:);
    % % %             %disp(['Point: ' num2str(XYZ(j,:)) '--cluster: ' num2str(mincl(1,:)) ': min d = ' num2str(min(d))])
    % % %
    % % %             if min(d) <= mydist,
    % % %
    % % %                 %mincl = XYZc(find(d==min(d)),:);
    % % %                 %disp(['Point: ' num2str(XYZ(j,:)) '--cluster: ' num2str(mincl(1,:)) ': min d = ' num2str(min(d))])
    % % %
    % % %
    % % %                 for k = 1:length(N)
    % % %
    % % %                     eval(['myvar = DB.' N{k} ';']);
    % % %                     if size(myvar,1) == size(XYZ,1)
    % % %
    % % %                         eval(['newcl(i).' N{k} '(ind,:) = DB.' N{k} '(j,:);'])
    % % %                     end
    % % %                 end
    % % %                 %disp('found one')
    % % %                 %keyboard
    % % %
    % % %                 ind = ind + 1;
    % % %
    % % %             end
    % % %
    % % %         end

    % % %         newcl(i).XYZ = newcl(i).XYZ';
    % % %         disp(['Found ' num2str(ind-1) ' points'])

    if ~isfield(newcl(i), 'x') & isfield(newcl(i), 'xyz')
        newcl(i).x = newcl(i).xyz(:, 1);
        newcl(i).y = newcl(i).xyz(:, 2);
        newcl(i).z = newcl(i).xyz(:, 3);
    end
    
    % add 'all' field with variable calling name
    if ~isempty(newcl(i).x) && exist('vname','var')
        newcl(i).alltask(1:size(newcl(i).x,1),1) = {vname};
    else
        newcl(i).alltask = {};
    end

end

% format in standard format for SPM orthviews
for i = 1:length(newcl)
    newcl(i).XYZmm = newcl(i).XYZ;
    newcl(i).XYZ = mm2voxel(newcl(i).XYZ,newcl(i).M,1)';
    newcl(i).Z = ones(1,size(newcl(i).XYZ,2));

end

end
