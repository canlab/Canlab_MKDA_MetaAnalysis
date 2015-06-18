% script read_database
% reads a neuroimaging coordinate database with fields
% and transforms T88 coordinates to MNI coordinates
% RULES
%
% 1     specify database name as dbname in workspace
% 2     database must be text file, tab delimited
% 3     1st row of database must contain number of columns as 1st entry
% 4     2nd row of database must be variable names
% 5     3rd - nth row of dbase contains data
% 6     Talairach coordinates are indicated by T88 in CoordSys variable
% 7     coordinate fields should be called x, y, and z
% 8     Do NOT use Z, or any other field name in clusters structure
%
% Note: On MAC, seems to read stuff saved as Windows Text in Excel
% 
% AN EXAMPLE FROM THE SADNESS DATABASE
% % 1 - save sheet as sheet name.txt, Windows format text
% e.g., sheet name is MastervNeutralOnly.txt
% clear variables
% clear
% define the name of the sheet to load
% dbname = 'MastervNeutralOnly.txt'
% read the database in, all columns
% read_database
% save a matlab .mat data file with these columns (variables) in it
% call it the same name as the sheet
% save MastervNeutralOnly
% to load the file again, just type load filename
% load MastervNeutralOnly
%
% Special Fields
% ------------ These field names have special meaning in the program ----
% 
% Subjects or N     : sample size
% FixedRandom       : fixed or random effects
% Subjective Weights : weighting vector based on FixedRandom and whatever
%   else you want to weight by; e.g., study reporting threshold
% x, y, z           : coordinates
% study or Study    : name of study
% Contrast          : unique indices (e.g., 1:k) for each independent
%                     contrast


clear a DBnames
global study

% not to be used in variable names
badstrings = '!@#$%^&*(){}[] ~`\|<>,.?/;:"''+=';

% -----------------------------------------------------------------------------
% * load database and read in first column - first row of first col is num columns
% -----------------------------------------------------------------------------

if ~exist('dbname') == 1, dbname = input('Enter name of input file: ','s');,end

a = textread(dbname,'%s%*[^\n]','delimiter','\t');
numc = str2num(a{1});

str = repmat('%s',1,numc);

% -----------------------------------------------------------------------------
% * define output string
% -----------------------------------------------------------------------------

outstr = ['[a1'];
for i = 2:numc
    outstr = [outstr ',a' num2str(i)];
end
outstr = [outstr ']'];

% -----------------------------------------------------------------------------
% * load full database as strings 
% -----------------------------------------------------------------------------

eval([outstr  ' = textread(dbname,str,''delimiter'',''\t'');'])

% -----------------------------------------------------------------------------
% * change strings to numbers, if possible
% -----------------------------------------------------------------------------
myl = length(a1) - 2;
for i = 1:numc
    
    eval(['vec = a' num2str(i) '(3:end);'])
    
    vecout = [];
    for j = 1:myl
        try
            vecout(j,1) = str2num(vec{j});
        catch
            break
        end
    end
    if isempty(vecout) | length(vecout) < length(vec), vecout = vec; end

    eval(['ColumnNames{i} = a' num2str(i) '{2};'])
    eval(['myname = a' num2str(i) '{2};'])
    
    % get rid of problematic characters
    myname = deblank(myname); 
    for j = 1:length(badstrings), myname(strfind(myname, badstrings(j))) = []; end
    
    DBnames{i} = myname;
    
    eval([myname ' = vecout;'])
    eval(['clear a' num2str(i)])
end
   
% -----------------------------------------------------------------------------
% * check x, y, z
% -----------------------------------------------------------------------------
if ~(exist('x') == 1 && exist('y') == 1 && exist('z') == 1)
    disp('Database must have columns named x, y, and z, which contain numeric coordinates.');
    error('One of these is missing.')
    
elseif iscell(x) || iscell(y) || iscell(z)
    
    disp('Database must have columns named x, y, and z, which contain numeric coordinates.');
    
    disp('Finding errors:')
    for i = 1:length(x)
        try
            xx(i, 1) = str2num(x{i});
            yy(i, 1) = str2num(y{i});
            zz(i, 1) = str2num(z{i});
        catch
            fprintf('Row %3.0f : Bad value for x, y, or z.  Values entered are : %s %s %s\n', i, x{i}, y{i}, z{i});
        end
    end
    
    error('One of these seems to have non-numeric values in at least one row (check for blank spaces/empty row at the end!)')
end

clear xt yt zt

% -----------------------------------------------------------------------------
% * convert from talairach to MNI, if we can
% -----------------------------------------------------------------------------
if exist('CoordSys') == 1 & exist('x') == 1 & exist('y') == 1 & exist('z') == 1
    fprintf(1,'Converting T88 coordinates to MNI (M. Brett transform): ')
    numt = 0;
    for i = 1:length(CoordSys)
        if strcmp(CoordSys{i},'T88')
            XYZt = tal2mni([x(i) y(i) z(i)]);
            x(i) = XYZt(1); y(i) = XYZt(2); z(i) = XYZt(3);
            numt = numt + 1;
        end
    end
    fprintf(1,'%3.0f Transformed\n',numt)
    
    fprintf(1,'Converting MNI coordinates to T88 for XYZtal array: ')
    numt = 0;
    for i = 1:length(CoordSys)
        if strcmp(CoordSys{i},'MNI') | strcmp(CoordSys{i},'spm96') | strcmp(CoordSys{i},'spm99')
            XYZt = mni2tal([x(i) y(i) z(i)]);
            xt(i) = XYZt(1); yt(i) = XYZt(2); zt(i) = XYZt(3);
            numt = numt + 1;
        else
            xt(i) = x(i); yt(i) = y(i); zt(i) = z(i);
        end
    end
    fprintf(1,'%3.0f Transformed\n',numt)

    XYZall = [x y z];
    XYZtal = [xt' yt' zt'];

end

% -----------------------------------------------------------------------------
% * convert T to Z, if all elements are present
% -----------------------------------------------------------------------------
if exist('Tscore') == 1 & exist('Zscore') == 1 & exist('N') == 1

    fprintf(1,'Converting Tscores to Zscores\n ')
    for i = 1:length(Tscore),
        if ~isnan(Tscore(i)) & isnan(Zscore(i))
            Zscore(i) = spm_t2z(Tscore(i),N(i)-1);
        end
    end

end


clear a, clear i, clear j, clear myl, clear myname, clear vec, clear vecout, clear numc
clear str, clear XYZt, clear numt
    

%DBnames = whos;
if ~exist('xyz') == 1 && exist('x', 'var') && exist('y', 'var') && exist('z', 'var'), xyz = [x y z]; end
for i = 1:length(DBnames)
    eval(['LEN = length(' DBnames{i} ');'])
    if exist('xyz', 'var')
        % Meta-database, has XYZ, include only if right size
        if LEN == size(xyz,1)
            eval(['DB.(DBnames{i}) = ' DBnames{i} ';'])
        end
    else
        % Other database; include all in struct.
        eval(['DB.(DBnames{i}) = ' DBnames{i} ';'])
    end
end

if ~exist('dbsave') == 1, dbsave = input('Enter name of file to save (without .mat extension): ','s'); end

eval(['save ' dbsave])



   