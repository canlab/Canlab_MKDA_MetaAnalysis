function Meta_interactive_table
% Make table output when you click on a voxel in orthviews.
%

mm = round(spm_orthviews('Pos')');

% table output
% ----------------------------------------
global DB
global testfield

go = 1;
if isempty(DB), try load SETUP DB, catch, disp('Cannot load DB...no table. Load DB and declare as global first!'), 
        fname = spm_get(1,'*mat','Select mat file with DB database structure variable:');,
        load(fname);
        if isempty(DB), go = 0;, end
        
    end, 
end

if isempty(testfield), try load SETUP testfield, catch, testfield = input('Cannot load testfield from SETUP.  Type name of field in DB to display: ','s'), go = 1;, end, end

if isempty(testfield), testfield = input('Cannot load testfield from SETUP.  Type name of field in DB to display: ','s'), go = 1;, end,
    
while ~isfield(DB,testfield), disp(['NO field called ' testfield]);
    disp(DB), testfield = input('Type field name: ','s');
end
    
if go
    d = distance(mm,[DB.x DB.y DB.z]);
    wh = find(d < DB.radius_mm);

    if isfield(DB,'N'), N = DB.N;, elseif isfield(DB,'Subjects'), N = DB.Subjects; else, N = NaN*zeros(size(DB.x));, end
    if isfield(DB,'Contrast'), con = DB.Contrast;, else, con = NaN*zeros(size(DB.x));, end

    fprintf(1,'Studies activating within %3.0f mm of \t%3.0f\t%3.0f\t%3.0f\t\n',DB.radius_mm,mm(1),mm(2),mm(3));
        
    fprintf(1,'Study\tx\ty\tz\t%s\tN\tContrast #\tpoint index\tdist\t\n',testfield);
    for i = 1:length(wh)
            pt = wh(i);
            fprintf(1,'%s\t%3.0f\t%3.0f\t%3.0f\t%s\t%3.0f\t%3.0f\t%3.0f\t%3.2f\t\n', ...
               DB.study{pt},DB.x(pt),DB.y(pt),DB.z(pt),DB.(testfield){pt},N(pt),con(pt),pt,d(pt));
    end
    
    fprintf(1,'\n');
    
    % summary table
    values = DB.(testfield); values = values(wh);   % values within x mm
    levels = unique(values);
    
    fprintf(1,'Level\tCount\tTotalContrasts\t%%Activating\tWeighted\tWeighted%%\n');
    for i = 1:length(levels)
        cnt = sum(strcmp(values,levels{i}));
        total = meta_count_contrasts(DB,testfield,levels{i});

        % find contrasts and then study weights for these points
        islevel = strcmp(DB.(testfield),levels{i});
        wh = find(d < DB.radius_mm & islevel);
        clear connums
        for ii=1:length(wh)
            connums(ii) = find(DB.connumbers == DB.Contrast(wh(ii)));
        end
        weighted = sum(DB.studyweight(unique(connums)));
        
        wh2 = find(islevel);
        clear connums
        for ii=1:length(wh2)
            connums(ii) = find(DB.connumbers == DB.Contrast(wh2(ii)));
        end
        weighted_tot = sum(DB.studyweight(unique(connums))); 

        fprintf(1,'%s\t%3.0f\t%3.0f\t%3.2f\t%3.2f\t%3.2f\t\n',levels{i},cnt,total,100*cnt./total,weighted,100.*weighted./weighted_tot);
    end
    
end
        