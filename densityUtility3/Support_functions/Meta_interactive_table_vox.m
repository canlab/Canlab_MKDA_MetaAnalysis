function Meta_interactive_table_vox(compareflag,varargin)
% Meta_interactive_table_vox(compareflag,[data matrix, volInfo struct])
% Make table output when you click on a voxel in orthviews.
%
% This version uses fields in DB.PP and computes voxel-based distances
% for consistency with Meta_Setup and other meta-analysis voxel-based
% functions.
%
% the compareflag option
% compares the computation from DB against saved images in DB.PP and Xi
% use in a Meta_Logistic directory; for debugging purposes



global DB
global testfield

vargs = {};
if length(varargin) > 0, vargs = varargin; end

% get mm and voxel coordinates for current spm_orthviews position
[VOL,coord,mm] = get_current_coords(DB);



% Set up required variables for table
[N,con,testfield] = check_fields(DB,testfield);



% get data from activated contrasts
% USE DB.PP images, or if data and volInfo are entered above, use those
[dat2,conindex,whcontrasts] = get_data_at_coordinate(DB,coord,vargs);


% run comparison from saved images
% must be in Meta_Logistic results directory
if compareflag
    compare_from_saved(dat2,DB,coord,testfield);
end



% Table header

fprintf(1,'Studies activating within %3.0f voxels of \t%3.0f\t%3.0f\t%3.0f\t\n',DB.radius,mm(1),mm(2),mm(3));

fprintf(1,'Study\tx\ty\tz\t%s\tN\tContrast #\tpoint index\tvox. dist\tmm dist\n',testfield);

% for summary table
values = DB.(testfield);
levels = unique(values);
w = DB.studyweight;
w = w ./ mean(w);
levelcnt = zeros(1,length(levels)); conwts = 0;

for i=1:length(whcontrasts)

    whpts = find(con == whcontrasts(i));    % indices in DB.point lists for this contrast

    whcon = conindex(i);    % index in DB.contrast lists for this contrast

    xyzmm = [DB.x(whpts) DB.y(whpts) DB.z(whpts)];

    % convert to voxels for voxel distance
    % transpose so that 3-voxel lists are not reoriented (produces wrong
    % transform otherwise!)
    xyz = mm2voxel(xyzmm',VOL(1),1);

    % voxel distances
    d = distance(coord,xyz);
    whinradius = find(d <= DB.radius);  % indices relative to this pt list of points for this contrast w/i radius

    % overall DB indices of points for this contrast w/i radius
    whpts2 = whpts(whinradius);

    if isempty(whpts2), 
        
        disp('Warning! Size or database mismatch. Contrast is in-radius, but no points within contrast are in-radius. Using closest:');, 
        disp(['check: ' DB.PP(whcon,:)])
        whinradius = find(d==min(d)); whinradius = whinradius(1);
        whpts2 = whpts(whinradius);
        keyboard
    
    end


    % save summary information
    for j = 1:length(levels)
        levelcnt(i,j) = double(any(strcmp(levels{j},values(whpts2))));  % any values of this level for this contrast
        conwts(i,1) = w(whcon);
    end

    for j = 1:length(whpts2)
        pt = whpts2(j);
        fprintf(1,'%s\t%3.0f\t%3.0f\t%3.0f\t%s\t%3.0f\t%3.0f\t%3.0f\t%3.2f\t%3.2f\t\n', ...
            DB.study{pt},DB.x(pt),DB.y(pt),DB.z(pt),DB.(testfield){pt},N(pt),con(pt),pt,d(whinradius(j)),distance(mm',xyzmm(whinradius(j),:)));
    end


end

fprintf(1,'\n');




% summary table

for i = 1:length(levels)
    [total(1,i),indic,wtotal(1,i),conindic(:,i)] = meta_count_contrasts(DB,testfield,levels{i});
    % returns total and weighted total

end

n = size(levelcnt,1);
cnt = sum(levelcnt,1);            % contrast count
prc = 100.*sum(levelcnt)./total;  % percentage of contrasts in each condition

wcnt = (levelcnt' * conwts)'; %weighted contrast counts
%./ sum(conwts);    % weighted percentage of contrasts
%wcnt = wprc .* total;            % weighted contrast counts
wprc = 100.*wcnt ./ wtotal;  % weighted percentage of contrasts



tabledat = [cnt;total;prc;wcnt;wtotal; wprc];
tablenms = {'Count' 'Total contrasts' 'Percentage' 'Wtd. count' 'Wtd. total' 'Wtd. perc.'};

% table header
fprintf(1,'Measure\tLevel\t\n\t\t');
for i = 1:length(levels)
    fprintf(1,'%s\t',levels{i});
end
fprintf(1,'\n');

for i = 1:length(tablenms)
    fprintf(1,'%s\t',tablenms{i});
    for j = 1:length(levels)
        fprintf(1,'%3.2f\t',tabledat(i,j));
    end
    fprintf(1,'\n');
    if i == 3, fprintf(1,'\n');, end
end
fprintf(1,'\n');

disp('Chi-square analyses based on frequency tables')
dat = tabledat(1:2,:);
% nos must be total minus yesses
dat(2,:) = dat(2,:) - dat(1,:);
wstrings = {'' ', Warning: expected counts < 5: P-value may be inaccurate; mapwise p-values use nonparametric test.'};
[chi2,df,p,sig,warn,freq_tableu,e,warn2] = chi2test(dat,'table'); warnstr = wstrings{warn+1};
fprintf(1,'Unweighted: chi2(%3.0f) = %3.2f, p = %3.4f %s\n',df,chi2,p,warnstr);

dat = tabledat(4:5,:);
% nos must be total minus yesses
dat(2,:) = dat(2,:) - dat(1,:);
[chi2,df,p,sig,warn,freq_tablew,e,warn2] = chi2test(dat,'table'); warnstr = wstrings{warn+1};
fprintf(1,'Weighted: chi2(%3.0f) = %3.2f, p = %3.4f %s\n',df,chi2,p,warnstr);

% nonparametric 
if exist('dat2','var')
    condf = indic2condf(conindic);
    %condf = condf(DB.connumbers);
    [chi2,df,p,sig,warn,freq_table,e,nonpar] = chi2test(full([dat2 condf]),'obs',w,1); 
    fprintf(1,'Weighted: chi2(%3.0f) = %3.2f, p = %3.4f %s Did nonpar: %3.0f\n',df,chi2,p,'Weighted with nonpar option.',nonpar);
end

% not completed, and won't work for weighted
%p = fishers_exact(freq_tableu);
%fprintf(1,'\nUnweighted: Fisher''s exact test: p = %3.4f\n',p);

return






% sub-functions




function [VOL,coord,mm] = get_current_coords(DB)
% get mm and voxel coordinates for current spm_orthviews position
if ~isfield(DB,'maskV'),
    error('No DB.maskV field, which is required.  See Meta_Setup to create this field and set up analysis.');
end

VOL = DB.maskV; VOL.M = VOL.mat;
coord = spm_orthviews('Pos');
mm = coord;
coord = mm2voxel(coord',VOL(1));
return






function [N,con,testfield] = check_fields(DB,testfield)
if isfield(DB,'N'), N = DB.N;, elseif isfield(DB,'Subjects'), N = DB.Subjects; else, N = NaN*zeros(size(DB.x));, end
if isfield(DB,'Contrast'), con = DB.Contrast;,
else, con = NaN*zeros(size(DB.x));,
    disp('Warning!  You must have a field called DB.Contrasts for the table function to work properly.');
end

if ~isfield(DB,'connumbers'),
    error('No DB.connumbers field, which is required.  See Meta_Setup to create this field and set up analysis.');
end

% Define testfield (field to display in table)
if isempty(testfield), try load SETUP testfield, catch, testfield = input('Cannot load testfield from SETUP.  Type name of field in DB to display: ','s'), go = 1;, end, end

if isempty(testfield), testfield = input('Cannot load testfield from SETUP.  Type name of field in DB to display: ','s'), go = 1;, end,

while ~isfield(DB,testfield), disp(['NO field called ' testfield]);
    disp(DB), testfield = input('Type field name: ','s');
end
return



function [dat2,conindex,whcontrasts] = get_data_at_coordinate(DB,coord,vargs)

if ~isempty(vargs)
    dat = vargs{1};
    volInfo = vargs{2};
    tmp = sub2ind(volInfo.dim(1:3),coord(1),coord(2),coord(3));
    whcol = find(volInfo.wh_inmask == tmp);
    if isempty(whcol)
        disp(['Coordinate is not in mask in volInfo.'])
        dat2 = []; conindex = []; whcontrasts = [];
        return
    else
        dat2 = dat(whcol,:)';
    end
        
else
    P2 = check_valid_imagename(DB.PP);
    dat2 = spm_get_data(P2,coord');
end
conindex = find(dat2);
whcontrasts = DB.connumbers(conindex);

return




function compare_from_saved(dat2,DB,coord,testfield)
try
    load SETUP Xi Xinms X names 
catch
    disp('Warning: compare will not work right unless you are the dir with SETUP.mat, which contains Xi and Xinms variables');
    disp('See Meta_Logistic for creating Xi and Xinms variables');
end

% if empty images specified, ncontrasts will be greater than nimgs
ncontrasts = size(Xi,2); nimgs = size(dat2,1);
doempty = ncontrasts - nimgs;
if doempty > 0
    addtoy = zeros(doempty,1);
end

fprintf(1,'\nFROM SAVED IMAGES:\nTotal count: %3.0f\n',sum(dat2))
cnt = Xi' * dat2; for i=1:length(Xinms), fprintf(1,'%s\t%3.0f\n',Xinms{i},cnt(i));, end
fprintf(1,'\n');
wh = DB.pointind(find(dat2)); %disp(DB.Study(wh))

for i=1:length(wh)
    mydist = coord;
    fprintf(1,'%s\t%3.0f\t%s\t\n',DB.Study{wh(i)},DB.Contrast(wh(i)),DB.(testfield){wh(i)});
end
fprintf(1,'\n');

y = dat2;

[b,stats,F,p,df1,df2] = meta_analyze_data(y,'X',Xi,'chi2',DB.studyweight,'table','names',Xinms);
[b,stats,F,p,df1,df2] = meta_analyze_data(y,'X',DB.X,'logistic',DB.studyweight,'table','names',names);

%print_matrix([y condf]);
% try
%     load SETUP condf
%     [chi2,chi2p,sig,df,tab] = nonparam_chi2(y,condf);
%     fprintf(1,'%s\t%3.2f\t%3.2f\t\n','nonparametric chi2: ',chi2,chi2p);
% catch
% end

% % set up weighted average calc
% w = DB.studyweight;
% w = w ./ mean(w);
% W = diag(w); % ./ sum(w));
% w(find(dat2))   % weights
% 
% sumxiu = sum(Xi,1);
% cnt = (Xi'*y)';
% Xi = Xi' * W;           % multiply this by data to get weighted avgs
% sumxi = sum(Xi,2);
% wcnt = (Xi*y)';
% 
% avg = wcnt; avg = avg ./ sumxi';
% 
% fprintf(1,'Unweighted counts and totals\n');
% for i =1:length(Xinms), fprintf(1,'%s\t',Xinms{i});, end,fprintf(1,'\n');
% for i =1:length(Xinms), fprintf(1,'%3.2f\t',cnt(i));, end,fprintf(1,'\n');
% for i =1:length(Xinms), fprintf(1,'%3.2f\t',sumxiu(i));, end,fprintf(1,'\n');
% 
% fprintf(1,'\nWeighted counts, totals, percentage\n');
% for i =1:length(Xinms), fprintf(1,'%s\t',Xinms{i});, end,fprintf(1,'\n');
% for i =1:length(Xinms), fprintf(1,'%3.2f\t',wcnt(i));, end,fprintf(1,'\n');
% for i =1:length(Xinms), fprintf(1,'%3.2f\t',sumxi(i));, end,fprintf(1,'\n');
% for i =1:length(Xinms), fprintf(1,'%3.4f\t',100 .* avg(i));, end,fprintf(1,'\n');
% fprintf(1,'\n');

% 
% % from Meta_Logistic
% 
% [b,dev,stats]=glmfit(X,[y ones(size(y))],'binomial','logit','off',[],w); % pvals are more liberal than Fisher's Exact!
% 
% % omnibus test - R^2 change test
% % --------------------------------------------------------
% sstot = y'*y;
% r2full = (sstot - (stats.resid' * stats.resid)) ./ sstot;
% dffull = stats.dfe;
% 
% [br,devr,statsr]=glmfit(ones(size(y)),[y ones(size(y))],'binomial','logit','off',[],w,'off');
% r2red = (sstot - (statsr.resid' * statsr.resid)) ./ sstot;
% dfred = statsr.dfe;
% 
% if r2full < r2red, fprintf(1,'Warning!'); r2red = r2full;,drawnow; fprintf(1,'\b\b\b\b\b\b\b\b'); end
% [F,op,df1,df2] = compare_rsquare_noprint(r2full,r2red,dffull,dfred);
% 
% fprintf(1,'Name\tBeta\tt-value\tp\t\n');
% 
% names = [{'Intercept'} names];
% for i = 1:length(b)
% 
%     fprintf(1,'%s\t%3.2f\t%3.2f\t%3.4f\t\n', ...
%         names{i},b(i),stats.t(i),stats.p(i));
% 
% end
% 
% % create condf - from meta_logistic_design
% % return condition function (dummy codes) for chi-square test
% [i,j] = find(Xi');
% [i,s] = sort(i); 
% condf = j(s);
% 
% [chi2,df,chi2p,sig,warn,tab,e,isnonpar] = chi2test([y condf],'obs',w,1);
% 
% % print table
% wstrings = {'' ', Warning: expected counts < 5, used nonparametric p-value.'};
% warnstr = wstrings{isnonpar+1};
% fprintf(1,'Weighted: chi2(%3.0f) = %3.2f, p = %3.4f %s\n',df,chi2,chi2p,warnstr);
% fprintf(1,'\n');

return




