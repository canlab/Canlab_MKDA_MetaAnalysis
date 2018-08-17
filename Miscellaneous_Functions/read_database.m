% script read_database
% reads a neuroimaging coordinate database with fields
% and transforms T88 coordinates to MNI coordinates
% When you set up your own database, several rules apply. Formatting the spreadsheet
% so it loads correctly can sometimes be the hardest part of running a
% meta-analysis.
%
% See https://canlabweb.colorado.edu/wiki/doku.php/help/meta/meta_analysis_mkda
% 
% And https://canlabweb.colorado.edu/wiki/doku.php/help/meta/database
%
% Here are some rules for setting up the file:
%   1     The variable "dbname" in the workspace should specify the name of your coordinate database file
%   2     The database must be a text file, tab delimited
%   3     The 1st row of the database must contain the number of columns as its 1st and only entry
%   4     The 2nd row of database must contain variable names (text, no spaces or special characters)
%   5     The 3rd - nth rows of the database contains data
%   6     Talairach coordinates are indicated by T88 in CoordSys variable
%   7     coordinate fields should be called x, y, and z
%   8     Do NOT use Z, or any other field name in clusters structure
% 
% The second row of your database should contain names for each variable
% you have coded.  Some variables should be named with special keywords,
% because they are used in the meta-analysis code. Other variables can be named 
% anything, as long as there are *no spaces or special characters* in the
% name (e.g., !@#$%^&*(){}[] ~`\/|<>,.?/;:"''+=). 
% Anything that you could not name a variable in Matlab will also not work
% here. 
% Here are the variable names with special meaning. They are case-sensitive:
%
% Subjects          : Sample size of the study to which the coordinate belongs
% FixedRandom       : Study used fixed or random effects. 
%                     Values should be Fixed or Random.
%                     Fixed effects coordinates will be automatically
%                     downweighted
% SubjectiveWeights : A coordinate or contrast weighting vector based on FixedRandom 
%                     and whatever else you want to weight by; e.g., study reporting threshold
%                     The default is to use FixedRandom only if available
% x, y, z           : X, Y, and Z coordinates
% study             : name of study
% Contrast          : unique indices (e.g., 1:k) for each independent
%                     contrast. This is a required variable!
%                     All rows belonging to the same contrast should
%                     (almost) always have the same values for every
%                     variable other than x, y, and z.
% CoordSys          : Values should be MNI or T88, for MNI space or Talairach space
%                     Talairach coordinates will be converted to MNI using
%                     Matthew Brett's transform
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
% If a variable dbsave is created with a text string, the DB structure will
% be saved in a .mat file with the name specified by the string.

clear a DBnames
global study

% not to be used in variable names
badstrings = '!@#$%^&*(){}[] ~`\/|<>,.?/;:"''+=';

% -----------------------------------------------------------------------------
% * load database and read in first column - first row of first col is num columns
% -----------------------------------------------------------------------------

if ~exist('dbname') == 1, dbname = input('Enter name of input file: ','s'); end

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
    
    if iscell(x)
        xx = cellfun(@str2num, x, 'UniformOutput', 0);
        fprintf('x has %d bad values.\nRows:\n', sum(cellfun(@isempty, xx)));
        find(cellfun(@isempty, xx))
    end
    
    if iscell(y)
        yy = cellfun(@str2num, y, 'UniformOutput', 0);
        fprintf('x has %d bad values.\nRows:\n', sum(cellfun(@isempty, yy)));
        find(cellfun(@isempty, yy))
    end
    
    
    if iscell(z)
        zz = cellfun(@str2num, z, 'UniformOutput', 0);
        fprintf('z has %d bad values.\nRows:\n', sum(cellfun(@isempty, zz)));
        find(cellfun(@isempty, zz))
        z(find(cellfun(@isempty, zz)))
    end
    

%     for i = 1:length(x)
%         try
%             xx(i, 1) = str2num(x{i});
%             yy(i, 1) = str2num(y{i});
%             zz(i, 1) = str2num(z{i});
%         catch
%             fprintf('Row %3.0f : Bad value for x, y, or z.  Values entered are : %s %s %s\n', i, x{i}, y{i}, z{i});
%         end
%     end
    
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

% Save, if dbsave var is entered

if ~exist('dbsave', 'var')
    % do not save; do nothing
    
else
    if isempty(dbsave)
        dbsave = input('Enter name of file to save (without .mat extension): ','s'); 
    end
end

eval(['save ' dbsave])



   