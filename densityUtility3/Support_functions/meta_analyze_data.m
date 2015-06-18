function varargout = Meta_analyze_data(y,varargin)
    % varargout = Meta_analyze_data(y,varargin)
    %
    % Inputs
    % =========================================================================
    % Weights:
    %'w'             % followed by weights
    % Analysis types:
    %'chi2'          % weighted chi-square
    %'logistic'      % weighted logistic regression
    %
    % Output strings:
    %'table'         % print table(s)
    % 'names'        % followed by names of regressors for table (logistic)
    %
    % Input types (for chi-square):
    %'X'             % followed by indicator matrix for chi2, design matrix for
    %                  logistic
    %'condf'         % followed by integer vector of cell assignments, for chi2
    %
    % by Tor Wager, May 2006
    %
    % Examples
    % =========================================================================
    % y = spm_get_data(DB.PP,[36 24 36]');  % get some data
    % load SETUP names                      % get regresor names
    %[b,stats,F,p,df1,df2] = meta_analyze_data(y,'X',DB.X,'logistic',DB.studyweight,'table','names',names);
    %
    % load SETUP Xi                         % get indicator matrix
    % [b,stats,F,p,df1,df2] = meta_analyze_data(y,'X',Xi,'chi2',DB.studyweight,'table');
    %
    % [chi2,df,p,sig,warn,tab,expected,isnonpar] = meta_analyze_data(dat(1:5,:)','X',Xi,'chi2','w',w,'table');
    
    % Programmers' notes
    % Tor Edited July 2011 to remove predictor rows that are zero in all conditions,
    % if any, during chi2 analysis. rare.
    %
    % Edited chi2 on Dec 2011 (tor) to check for full/double matrices and
    % enforce data type.
    
    
    dotable = 0;

    for i = 1:length(varargin)
        if isstr(varargin{i})
            switch lower(varargin{i})
                % reserved keywords
                case 'table', dotable = 1;
                case 'chi2', analysistype = 'chi2';
                case 'logistic', analysistype = 'logistic';
                case 'x', inputtype = 'X';, X = varargin{i+1};
                case 'condf', inputtype = 'condf';, condf = varargin{i+1};

                    % functional commands
                case 'w', w = varargin{i+1};
                case 'names', names = varargin{i+1};

                otherwise, warning(['Unknown input string option:' varargin{i}]);
            end
        end
    end

    if ~(exist('w')==1), w = ones(size(y, 1), 1); end
    w = w ./ mean(w);

    switch analysistype
        case 'chi2'
            % ========================================================
            % Chi-squared
            % ========================================================

            switch inputtype
                case 'X'
                    % create condf - from meta_logistic_design
                    % return condition function (dummy codes) for chi-square test
                    [i,j] = find(X);
                    [i,s] = sort(i);
                    condf = zeros(size(X, 1), 1); % tor edited july 2011 to handle zeros
                    condf(i) = j(s);
                    wh = condf ~= 0; % wh to include; usually all...
                    
                    if size(condf,1) ~= size(y,1), error('Condition function does not match data.  Design matrix X is probably not a suitable indicator matrix.');, end

                    if dotable && exist('names')==1
                        frequency_tables(y,w,X,names);
                    elseif dotable
                        disp('Frequency table omitted: no names.');
                    end

                case 'condf'
                    % do nothing
                    if dotable,
                        disp('Frequency table omitted: Enter indicator matrix to get this.');
                    end

                otherwise
                    disp('Unknown input type. Must be ''X'' or ''condf''');
            end

            nvars = size(y,2);
            doupdate = 0;
            if nvars > 100
                doupdate = 1;
                progressbar('init');
                fprintf(1,'Running %3.0f tests: 00000',nvars);

            end

            [chi2, chi2p, sig, warn, isnonpar] = deal(zeros(1, nvars));
            tab = cell(1, nvars);
            e = cell(1, nvars);
            
            for i = 1:nvars
                %[chi2(i),df,chi2p(i),sig(i),warn(i),tab{i},e{i},isnonpar(i)] = chi2test([y(:,i) condf],'obs',w,1);
                % tor edited july 2011 to omit condf == 0
                
                inputdat = double([y(wh,i) condf(wh)]);
                if issparse(inputdat), inputdat = full(inputdat); end
                
                [chi2(i),df,chi2p(i),sig(i),warn(i),tab{i},e{i},isnonpar(i)] = chi2test(inputdat,'obs',w(wh),1);
                
                if doupdate && mod(i,10)==0, 
                    fprintf('\b\b\b\b\b%05d',i); 
                    progressbar('update',100*i./nvars);
                end

            end

             if doupdate
                  fprintf('\b\b\b\b\b%05d: Done!',nvars);
             end
             
            varargout{1} = chi2; varargout{2} = df; varargout{3} = chi2p; varargout{4} = sig;
            varargout{5} = warn; varargout{6} = tab;  varargout{7} = e; varargout{8} = isnonpar;

            if dotable
                % print table
                disp('Chi-square table for first data vector')
                wstrings = {'' ', Warning: expected counts < 5.'};
                warnstr = wstrings{warn(1)+1};
                fprintf(1,'Weighted: chi2(%3.0f) = %3.2f, p = %3.4f %s\n',df(1),chi2(1),chi2p(1),warnstr);
                fprintf(1,'\n');
            end

        case 'logistic'
            % ========================================================
            % Logistic Regression
            % ========================================================

            [b,dev,stats]=glmfit(X,[y ones(size(y))],'binomial','logit','off',[],w); % pvals are more liberal than Fisher's Exact!

            % omnibus test - R^2 change test
            % --------------------------------------------------------
            sstot = y'*y;
            r2full = (sstot - (stats.resid' * stats.resid)) ./ sstot;
            dffull = stats.dfe;

            [br,devr,statsr]=glmfit(ones(size(y)),[y ones(size(y))],'binomial','logit','off',[],w,'off');
            r2red = (sstot - (statsr.resid' * statsr.resid)) ./ sstot;
            dfred = statsr.dfe;

            if r2full < r2red, fprintf(1,'Warning!'); r2red = r2full;,drawnow; fprintf(1,'\b\b\b\b\b\b\b\b'); end
            [F,op,df1,df2] = compare_rsquare_noprint(r2full,r2red,dffull,dfred);

            if dotable
                fprintf(1,'Name\tBeta\tt-value\tp\t\n');

                names = [{'Intercept'} names];
                for i = 1:length(b)

                    fprintf(1,'%s\t%3.2f\t%3.2f\t%3.4f\t\n', ...
                        names{i},b(i),stats.t(i),stats.p(i));

                end
            end

            varargout{1} = b; varargout{2} = stats; varargout{3} = F; varargout{4} = op;
            varargout{5} = df1; varargout{6} = df2;


        otherwise
            % ========================================================
            % ???
            % ========================================================
            disp('Unknown analysis type.');
    end

    return





function frequency_tables(y,w,Xi,Xinms)

    w = w ./ mean(w);
    W = diag(w); % ./ sum(w));
    %%w(find(y), :);   % weights

    sumxiu = sum(Xi,1);
    cnt = (Xi'*y)';
    Xi = Xi' * W;           % multiply this by data to get weighted avgs
    sumxi = sum(Xi,2);
    wcnt = (Xi*y)';

    avg = wcnt; avg = avg ./ sumxi';

    fprintf(1,'Unweighted counts and totals\n');
    for i =1:length(Xinms), fprintf(1,'%s\t',Xinms{i});, end,fprintf(1,'\n');
    for i =1:length(Xinms), fprintf(1,'%3.2f\t',cnt(i));, end,fprintf(1,'\n');
    for i =1:length(Xinms), fprintf(1,'%3.2f\t',sumxiu(i));, end,fprintf(1,'\n');

    fprintf(1,'\nWeighted counts, totals, percentage\n');
    for i =1:length(Xinms), fprintf(1,'%s\t',Xinms{i});, end,fprintf(1,'\n');
    for i =1:length(Xinms), fprintf(1,'%3.2f\t',wcnt(i));, end,fprintf(1,'\n');
    for i =1:length(Xinms), fprintf(1,'%3.2f\t',sumxi(i));, end,fprintf(1,'\n');
    for i =1:length(Xinms), fprintf(1,'%3.0f\t',100 .* avg(i));, end,fprintf(1,'\n');
    fprintf(1,'\n');

    return






