function [count,studies,studycount] = countemstudies(varname,inputstr,varname2,inputstr2,varname3,inputstr3)
% function [count,studies,studycount] = countem(varname,inputstr,varname2,inputstr2,etc.)
% 	up to three conditions only right now
%
% reserved command: 'writein' for varname
%					then type your conditions, i.e. z > 4, without quotes.
%
% otherwise, varname is the name of a cell array of string vars, no quotes.
%
% Tor Wager, 03/06/01
%
% To load columns:
% coordfile = 'emotrev.txt'
% [study,scan,roi,coordsys,numsubj,gender,stimuli,method,emotion,valence,target,ref,other, ...
%   cogtype,cogdemand,x,y,z,zscore,region,colors,mytext] ...
%    = textread(coordfile,'%s%s%s%s%n%s%s%s%s%s%s%s%s%s%s%n%n%n%n%s%s%s');
global study
count = 0;

if isempty(study), error('Error in countemstudies: declare study as global variable in workspace/script.'),end

switch nargin
case 2
	
	if strcmp(varname,'writein')			% then inputstr will be a test string
		count = sum(inputstr);
		studies = unique(study(inputstr));
		studycount = size(studies,1);
	else
		test = strcmp(varname,inputstr);
		count = sum(test);
		studies = unique(study(test));
		studycount = size(studies,1);
	end
	
case 4
		
	if strcmp(varname,'writein')			% then inputstr will be a test string
		one = inputstr;
	else
		one = strcmp(varname,inputstr);
	end
		
	if strcmp(varname2,'writein')			% then inputstr will be a test string
		two = inputstr2;
	else
		two = strcmp(varname2,inputstr2);
	end
	
	%disp(['count for condition one is ' num2str(sum(one))])
	%disp(['count for condition two is ' num2str(sum(two))])
	count = sum(one & two);
	studies = unique(study(one & two));
	studycount = size(studies,1);
	%disp(['the AND conjunction is 	' num2str(count)])
	
	
case 6
		
	if strcmp(varname,'writein')			% then inputstr will be a test string
		one = inputstr;
	else
		one = strcmp(varname,inputstr);
	end
		
	if strcmp(varname2,'writein')			% then inputstr will be a test string
		two = inputstr2;
	else
		two = strcmp(varname2,inputstr2);
	end
	
	if strcmp(varname3,'writein')			% then inputstr will be a test string
		three = inputstr;
	else
		three = strcmp(varname3,inputstr3);
	end

	%disp(['count for condition one is 	' num2str(sum(one))])
	%disp(['count for condition two is 	' num2str(sum(two))])
	%disp(['count for condition three is ' num2str(sum(three))])
	%disp('_______________________________________')
	%disp(['one AND two are 	' num2str(sum(one & two))])
	%disp(['one AND three are 	' num2str(sum(one & three))])
	%disp(['two AND three are 	' num2str(sum(two & three))])
	%disp('_______________________________________')
	count = sum(one & two & three);
	studies = unique(study(one & two & three));
	studycount = size(studies,1);
	%disp(['the AND conjunction of all:	' num2str(count)])
	
	
otherwise error('Input 2, 4, or 6 arguments only.')
end

return
