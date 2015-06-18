function OUT = meta_analysis_table(XYZlist,r,varargin)
% OUT = meta_analysis_table(XYZlist,r,[master file name])
% takes a list of XYZ mm points, and, for each point you specify, makes a table of
% other studies in the meta-analysis files that report peaks within r mm (radius) of
% the point.
% Looks for a file on the path called "meta_analysis_master_file.mat",
% unless you specify another file name as an optional argument.  The file
% should contain one or more meta-analysis database structures, named *DB*,
% which are created with database2clusters.
%
% This function calls other useful functions that can also be used
% stand-alone:
%   cl = database_point_probe(PAINDB2,[-42 0 2],10,fl);
%   dbcluster_point_table(cl)

if size(XYZlist, 2) ~= 3
    XYZlist = XYZlist';
end

if size(XYZlist, 2) ~= 3
    error('XYZlist must be n x 3 matrix of x, y, and z coordinates.')
end

% get file name and load
name = 'meta_analysis_master_file.mat'; 
if length(varargin) > 0, name = varargin{1}; end
name2 = which(name);
if isempty(name2), error(['Cannot find ' name]), end
name = name2; clear name2;
load(name, '*DB*');

% get DB structures
d = whos('*DB*');
fprintf(1,'Found %3.0f database structures:\n', length(d));

for i = 1:length(d), fprintf(1,'%s\t',d(i).name); end
fprintf(1,'\nRadius is %3.0f mm\n', r)

% run table for each

for i = 1:length(d)
    
    eval(['myDB = ' d(i).name ';'])
    
    if isstruct(myDB)
        
        disp(['Database: ' d(i).name])
        disp('_____________________________')
        
        cl = database_point_probe(myDB, XYZlist, r);
        
        % summary
        disp('Summary of first coordinate/cluster');
        nstudies = length(unique(myDB.Study));
        if ~isfield(cl(1), 'Study')
            nyes = 0; 
        else
            nyes = length(unique(cl(1).Study));
        end
        OUT.dbname{i} = d(i).name;
        OUT.nyes(i) = nyes;
        OUT.nstudies(i) = nstudies;
        OUT.percentage_by_database(i) = 100*(nyes./nstudies);
        fprintf('%3.0f of %3.0f total studies, %3.0f%%\n', nyes, nstudies, OUT.percentage_by_database(i));
        
        disp('_____________________________')
        disp(' ')
        
    end
    
end

return

