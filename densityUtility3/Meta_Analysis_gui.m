%
%
% Tor Wager & Brencho
%
% Thanks to Tom Nichols for the excellent GUI shell!

%-----------------------------functions-called------------------------
%
%-----------------------------functions-called------------------------

function varargout = Meta_Analysis_gui(Action,varargin)
    % global variables we need for this shell

    global DB
    global cl

    %-Format arguments
    %-----------------------------------------------------------------------
    if nargin == 0, Action='Init'; end


    switch lower(Action)

        case lower('Init')
            %=======================================================================

            %clc
            %BrainVowager_defaults;
            Meta_Analysis_gui('AsciiWelcome')
            Meta_Analysis_gui('CreateMenuWin')

            % load DB

            if(~isempty(DB))
                disp('Using DB already in memory.');
            elseif(exist('DB.mat', 'file'))
                disp('loading DB.mat: Database structure in this file (in current dir) is available.');
                load DB;
            elseif(exist('SETUP.mat', 'file'))
                disp('loading SETUP.mat: Database structure in this file (in current dir) is available.');
                load SETUP;
            else
                disp('You need to go to the directory containing meta-analysis files or results');
                disp('this directory should have information about the analysis stored in the file DB.mat or SETUP.mat')
                disp('No DB file in current directory.');
                DB = [];
                fprintf(1,'\n')
            end
            
            if(exist('MC_Info.mat', 'file'))
                load('MC_Info');
                set_FWE_results_options(MC_Setup);
            end

            varargout{1} = DB;

        case lower('AsciiWelcome')
            %=======================================================================
            disp( 'Welcome to the imaging meta-analysis gui.  Written by Tor Wager, Jan 2006')
            fprintf('\n')

        case lower('Ver')
            %=======================================================================
            varargout = {'SCNlab Meta-Analysis Menu'};



        case lower('CreateMenuWin')
            %=======================================================================
            close(findobj(get(0,'Children'),'Tag','Meta_Analysis_gui Menu'))


            %-Initialize Meta_Analysis_gui menu window
            %-----------------------------------------------------------------------
            [F, winwid, winh] = Meta_Analysis_gui('initFigure');


            % default button sizes and positions, etc.

            topbutton = winh-100;        % y location of top button
            butspace = 30;               % y spacing of buttons

            fullbutxy = [160 25];       % full-length button width and height
            halfbutxy = [80 25];        % (left-hand) half-width button w and h
            rightbutxy = halfbutxy;
            rightbutx = 110+halfbutxy(1)+5;  % right-hand button start x





            %-Frames and text
            %-----------------------------------------------------------------------
            axes('Position',[0 0 80/winwid winh/winh],'Visible','Off')
            text(0.5,0.475,'Meta-Analysis Toolbox',...
                'FontName','Times','FontSize',36,...
                'Rotation',90,...
                'VerticalAlignment','middle','HorizontalAlignment','center',...
                'Color',[1 1 1]*.6);

            text(0.2,0.96,'SCN Lab',...
                'FontName','Times','FontSize',16,'FontAngle','Italic',...
                'FontWeight','Bold',...
                'Color',[1 1 1]*.6);

            uicontrol(F,'Style','Frame','Position',[095 005 winwid-100 winh - 30],...
                'BackgroundColor',Meta_Analysis_gui('Color'));  % colored frame
            uicontrol(F,'Style','Frame','Position',[105 015 winwid-120 winh - 50]);  % inner gray frame

            %-Buttons to launch Meta_Analysis_gui functions
            %-----------------------------------------------------------------------

            % -------------------------------------------
            % Section - Random Effects Analysis
            uicontrol(F,'Style','Text',...
                'String','Setup','FontSize',14,...
                'HorizontalAlignment','Center',...
                'Position',[115 topbutton+30 fullbutxy],...
                'ForegroundColor','y','FontWeight','b');
            % -------------------------------------------

            % Read txt
            str = 'read_database;';                                    % callback function
            str = Meta_Analysis_gui('ExpandString',str);               % display then execute
            uicontrol(F,'String','Read txt',...
                'Position',[110 topbutton-(1-1)*butspace halfbutxy],...
                'CallBack',str,...
                'Interruptible','on',...
                'ForegroundColor','k','FontWeight','b');

            % Load
            spmg = 'P = spm_get(1,''*mat'',''Select DB or SETUP mat file:'');';
            str = [spmg ',load(P);'];                                  % callback function
            str = Meta_Analysis_gui('ExpandString',str);               % display then execute
            uicontrol(F,'String','Load',...
                'Position',[rightbutx topbutton-(1-1)*butspace halfbutxy],...
                'CallBack',str,...
                'Interruptible','on',...
                'ForegroundColor','k','FontWeight','b');

            % Meta_Setup
            str = 'DB = Meta_Setup(DB);';                              % callback function
            str = Meta_Analysis_gui('ExpandString', str);              % display then execute
            uicontrol(F,'String','Meta Setup',...
                'Position',[110 topbutton-(2-1)*butspace halfbutxy],...
                'CallBack',str,...
                'Interruptible','on',...
                'ForegroundColor','k','FontWeight','b');

            % Meta FWE Setup
            uicontrol(F,'String','Meta FWE Setup',...
                'Position',[rightbutx topbutton-(2-1)*butspace halfbutxy],...
                'CallBack', @FWE_setup,...
                'Interruptible','on',...
                'ForegroundColor','k','FontWeight','b');

            % Select Contrasts
            str = 'DB = Meta_Select_Contrasts(DB);';                   % callback function
            str = Meta_Analysis_gui('ExpandString', str);              % display then execute
            uicontrol(F,'String','Select Contrasts',...
                'Position',[110 topbutton-(3-1)*butspace fullbutxy],...
                'CallBack',str,...
                'Interruptible','on',...
                'ForegroundColor','k','FontWeight','b');


            % -------------------------------------------
            % Next section - Analysis
            uicontrol(F,'Style','Text',...
                'String','Analysis','FontSize',14,...
                'HorizontalAlignment','Center',...
                'Position',[115 topbutton-(4-1)*butspace fullbutxy],...
                'ForegroundColor','y','FontWeight','b');
            % -------------------------------------------


%             % Activation
%             str = 'DB = Meta_Activation(DB);';                           % callback function
%             str = Meta_Analysis_gui('ExpandString',str);                 % display then execute
%             uicontrol(F,'String','Activation',...
%                 'Position',[110 topbutton-(5-1)*butspace fullbutxy],...
%                 'CallBack',str,...
%                 'Interruptible','on',...
%                 'ForegroundColor','k','FontWeight','b');
% 
            % Monte Carlo
            str = 'Meta_Activation_FWE(''mc'',15000);';                           % callback function
            str = Meta_Analysis_gui('ExpandString',str);                 % display then execute
            uicontrol(F,'String','FWE Monte Carlo',...
                'Position',[110 topbutton-(5-1)*butspace fullbutxy],...
                'CallBack',str,...
                'Interruptible','on',...
                'ForegroundColor','k','FontWeight','b');

            % Logistic
            %     spmg = 'maskimg = spm_get(1,''*img'',''Select mask image file:''); ';
            %     str1 = 'DB = Meta_Logistic(DB,2,maskimg);';                             % callback function
            %     str = [spmg str1];
            %     str = Meta_Analysis_gui('ExpandString',str);                 % display then execute
            %     uicontrol(F,'String','Logistic',...
            %         'Position',[110 topbutton-(6-1)*butspace fullbutxy],...
            %         'CallBack',str,...
            %         'Interruptible','on',...
            %         'ForegroundColor','k','FontWeight','b');

            % Chi-square
            spmg = 'maskimg = spm_get(1,''*img'',''Select mask image file:''); ';
            str1 = 'DB = Meta_Chisq(DB,2,maskimg);';                             % callback function
            str = [spmg str1];
            str = Meta_Analysis_gui('ExpandString',str);                 % display then execute
            uicontrol(F,'String','Chi-square',...
                'Position',[110 topbutton-(6-1)*butspace fullbutxy],...
                'CallBack',str,...
                'Interruptible','on',...
                'ForegroundColor','k','FontWeight','b');

            % Specificity - IN PROGRESS
            str = 'DB = Meta_Specificity(DB)' ;                         % callback function
            str = Meta_Analysis_gui('ExpandString', str);               % display then execute
            uicontrol(F,'String','Specificity',...
                'Position',[110 topbutton-(7-1)*butspace fullbutxy],...
                'CallBack',str,...
                'Interruptible','on',...
                'ForegroundColor','k','FontWeight','b');


            % -------------------------------------------
            % Next section - Results
            uicontrol(F,'Style','Text',...
                'String','Results','FontSize',14,...
                'HorizontalAlignment','Center',...
                'Position',[115 topbutton-(8-1)*butspace fullbutxy],...
                'ForegroundColor','y','FontWeight','b');
            % -------------------------------------------

            % Get blobs
            %str1 = ['pimg = spm_get(1,''*img'',''Select p-value image.'',pwd);'];
            %str = 'cl = pmap_threshold; bar_interactive;';

            str = Meta_Analysis_gui('setupResultsView');
            %str = Meta_Analysis_gui('ExpandString', str);                  % set up display then execute string for callback
            uicontrol(F,'String','Get Blobs',...
                'Position',[110 topbutton-(9-1)*butspace fullbutxy],...
                'CallBack',str,...
                'Interruptible','on',...
                'ForegroundColor','k','FontWeight','b');

            % Blob Display
            if exist('scn_roi_gui', 'file')
                str = 'scn_roi_gui;';                                           % callback function
            else
                str1 = 'if isempty(cl), disp([''Load a clusters cl.mat file first!'']), return, end;';
                str2 = 'cluster_orthviews(cl,{[1 0 0]}); set(gcf,''WindowButtonUpFcn'',''Meta_interactive_table;'')';
                str = [str1 str2];
            end

            str = Meta_Analysis_gui('ExpandString', str);                   % display then execute
            uicontrol(F,'String','Blob Display Tool',...
                'Position',[110 topbutton-(10-1)*butspace fullbutxy],...
                'CallBack',str,...
                'Interruptible','on',...
                'ForegroundColor','k','FontWeight','b');

            % FWE results
%             str = 'Meta_Activation_FWE(''results'');';
%             str = Meta_Analysis_gui('ExpandString', str);                   % display then execute
            uicontrol(F,'String','FWE Results',...
                'Position',[110 topbutton-(11-1)*butspace fullbutxy],...
                'CallBack', @FWE_results_callback,...
                'Interruptible','on',...
                'ForegroundColor','k','FontWeight','b');

            if(exist('MC_Setup', 'var'))
                dropdown_options = FWE_results_dropdown_options(MC_Setup);
            else
                dropdown_options = {'Run Meta FWE Setup'};
            end
            uicontrol(F,'String',dropdown_options,...
                'Position',[110 topbutton-(12-1)*butspace fullbutxy],...
                'Style', 'popupmenu', ...
                'Tag', 'FWE results popup', ...
                'ForegroundColor','k','FontWeight','b');
            
            
            checkbox_indent = 20;
            uicontrol(F, 'String', 'Height', ...
                'Style', 'checkbox', ...
                'Tag', 'FWE results checkbox', ...
                'UserData', 'height', ...
                'Value', 1, ...
                'Position',[(110+checkbox_indent) topbutton-(13-1)*butspace (fullbutxy(1)-checkbox_indent) fullbutxy(2)]);
            uicontrol(F, 'String', 'Cluster 95/.05', ...
                'Style', 'checkbox', ...
                'Tag', 'FWE results checkbox', ...
                'UserData', 'lenient', ...
                'Value', 1, ...
                'Position',[(110+checkbox_indent) topbutton-(14-1)*butspace (fullbutxy(1)-checkbox_indent) fullbutxy(2)]);
            uicontrol(F, 'String', 'Cluster 99/.01', ...
                'Style', 'checkbox', ...
                'Tag', 'FWE results checkbox', ...
                'UserData', 'medium', ...
                'Value', 1, ...
                'Position',[(110+checkbox_indent) topbutton-(15-1)*butspace (fullbutxy(1)-checkbox_indent) fullbutxy(2)]);
            uicontrol(F, 'String', 'Cluster 99.9/.001', ...
                'Style', 'checkbox', ...
                'Tag', 'FWE results checkbox', ...
                'UserData', 'stringent', ...
                'Value', 1, ...
                'Position',[(110+checkbox_indent) topbutton-(16-1)*butspace (fullbutxy(1)-checkbox_indent) fullbutxy(2)]);
            

            set(F,'Pointer','Arrow','Visible','on')



        case lower('Color')
            %=======================================================================
            % Meta_Analysis_gui('Color')
            %-----------------------------------------------------------------------
            % %-Developmental livery
            % varargout = {[0.7,1.0,0.7], 'Lime Green'};
            %-Distribution livery
            varargout = {[.8 0.5 .2], 'Purple'};


        case lower('ExpandString')
            %=======================================================================
            % Meta_Analysis_gui('ExpandString')
            % Expand an action button callback string (a command string to be
            % evaluated)
            % so that it first displays the command, and then executes it
            %-----------------------------------------------------------------------
            str = varargin{1}; str2 = [];
            for i = 1:length(str)
                if str(i) == ''''
                    str2(end+1) = '''';
                    str2(end+1) = '''';
                else
                    str2(end+1) = str(i);
                end
            end

            %str = ['disp(''' char(str2) '''), ' str ];     % display then execute

            str = ['Meta_Analysis_gui(''executeCommand'',''' str2 ''');'];
            varargout = {str};


        case lower('initFigure')
            %=======================================================================
            % [F, winwid, winh] = Meta_Analysis_gui('initFigure')
            %-----------------------------------------------------------------------
            % Get the position of the main BrainVowager menu, or if
            % not available, default screen pos.

            % default sizes, etc.
            S = get(0,'ScreenSize');

            winwid = 300;               % window width
            winh = 580;                 % window height
            pos = [S(3)/2+150,S(4)/2-140,winwid,winh];  % default

            h = findobj('Tag','BrainVowager_gui Menu');
            if ~isempty(h),
                pos = get(h,'Position');
                winwid = pos(3); winh = pos(4);
                pos(1) = pos(1) + winwid;   % put next to main figure
            end

            %-Open Meta_Analysis_gui menu window
            %----------------------------------------------------------------------

            F = figure('Color',[1 1 1]*.8,...
                'Name',Meta_Analysis_gui('Ver'),...
                'NumberTitle','off',...
                'Position',pos,...
                'Resize','off',...
                'Tag','Meta_Analysis_gui Menu',...
                'Pointer','Watch',...
                'MenuBar','none',...
                'Visible','off');

            varargout{1} = F; varargout{2} = winwid; varargout{3} = winh;


        case lower('executeCommand')
            %=======================================================================
            % Meta_Analysis_gui('executeCommand',str)
            %-----------------------------------------------------------------------
            % Display a text header highlighting the command to be run

            fprintf(1,'\n______________________________________________\n')
            fprintf(1,'Meta Analysis tool   c 2006 Tor Wager\n');
            fprintf(1,'Executing command string:\n');
            fprintf(1,'%s\n',varargin{1});
            fprintf(1,'______________________________________________\n')
            eval(varargin{1});




        case lower('setupResultsView')
            %=======================================================================
            % callbackstr = Meta_Analysis_gui('setupResultsView')
            %-----------------------------------------------------------------------
            % Set up callback string for results view 'wizard'

            str1 = sprintf('class_avg_images = []; Xinms = [];\n');
            str2 = sprintf('try, load SETUP class_avg_images Xinms, \n   disp(class_avg_images);,disp(Xinms),\ncatch,\nend\n');
            str3 = sprintf('cl = pmap_threshold;\nbar_interactive(class_avg_images,Xinms);\n');
            varargout{1} = [str1 str2 str3];

        otherwise
            %=======================================================================
            error('Unknown action string')
            %======================================================================
    end
end

%=======================================================================
% Callbacks
%-----------------------------------------------------------------------
function FWE_results_callback(src, eventdata)
    useparams = {};
    
    checkboxes = findobj('Tag', 'FWE results checkbox');
    for i=1:length(checkboxes)
        if(get(checkboxes(i), 'Value') == 1)
            useparams{end+1} = get(checkboxes(i), 'UserData');
        end
    end
    
    popuph = findobj('Tag', 'FWE results popup');
    [maptype, contrast_num] = parse_results(popuph)
    if(~strcmp(maptype, 'act'))
        useparams{end+1} = maptype;
    end
    if(~isempty(contrast_num))
        useparams{end+1} = 'contrast';
        useparams{end+1} = contrast_num;
    end
    
    Meta_Activation_FWE('results', 1, useparams{:});
end

function FWE_setup(src, eventdata)
    global DB; % no getting around this without refactoring whole file... DB is not initialized when callback is setup
    MC_Setup = Meta_Activation_FWE('setup', DB);
    set_FWE_results_options(MC_Setup);
end

%=======================================================================
% Local functions
%-----------------------------------------------------------------------
function set_FWE_results_options(MC_Setup)
    popuph = findobj('Tag', 'FWE results popup');
    set(popuph, 'String', FWE_results_dropdown_options(MC_Setup));
end

function result_options = FWE_results_dropdown_options(MC_Setup)
    result_options = cell(1, 1 + 2*length(MC_Setup.connames));
    result_options{1} = 'Overall';
    for i = 1:length(MC_Setup.connames)
        result_options{2*i} = [MC_Setup.connames{i} '-pos'];
        result_options{(2*i)+1} = [MC_Setup.connames{i} '-neg'];
    end
end

function [maptype, contrast_num] = parse_results(popuph)
    selected = get(popuph, 'Value');
    if(selected == 1)
        maptype = 'act';
        contrast_num = [];
    else
        popup_options = get(popuph, 'String');
        selected_contrast_string = popup_options{selected};
        if(~isempty(strfind(selected_contrast_string, 'pos')))
            maptype = 'poscon';
        else
            maptype = 'negcon';
        end
        
        contrast_num = floor(selected / 2);
    end
end