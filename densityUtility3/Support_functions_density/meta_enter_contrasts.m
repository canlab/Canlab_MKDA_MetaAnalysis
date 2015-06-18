function [con,conname,contype,conweights] = meta_enter_contrasts(ti,tasknms)
%
% [con,conname,contype,conweights] = meta_enter_contrasts(ti,tasknms)
% [xxx,connames,xxx,contrasts] = meta_enter_contrasts(Xi,Xinms)

contype = 'contrast';

nconds = size(ti,2);

conweights = []; conname = {};
go = 1;

while go
    
    if nargin > 1
        for i = 1:length(tasknms)
            fprintf(1,'%03d  %s\n',i,tasknms{i});
        end
    end
    
    tmp = input(['Enter contrast across ' num2str(nconds) ' conditions in [], return to quit: ']);
    
    if isempty(tmp)
        go = 0;
    else
        if isempty(tmp) || length(tmp) ~= nconds, disp('Invalid contrast.');,
        else
            conweights(end+1,:) = tmp;
            conname{end+1} = input('Enter short name for this contrast, no spaces or special chars: ','s');
        end
    end
end

con = (conweights * ti')';

return
