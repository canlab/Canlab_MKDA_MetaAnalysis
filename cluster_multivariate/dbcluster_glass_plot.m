function dbcluster_glass_plot(newcl,mycol,varargin)
% dbcluster_glass_plot(newcl,colors,[select by field],[cell vector of conditions to select])
% tor wager
% newcl = clusters structure, output of database2clusters.m
%
% Creates 3 views of glass brain with points in each cluster plotted
% in different symbols
%
% Needs x, y, and z fields for each cluster, with points
% colors can be empty
%
% Example:
% dbcluster_glass_plot(EM,[],'emotion',{'happy' 'sad' 'fear' 'angry' 'disgust'})
% dbcluster_glass_plot(PAIN,{'bo' 'gs'},'Right_vs_Left',{'Left' 'Right'})
%
% Add to existing plots:
% dbcluster_glass_plot(PAIN,{'bo' 'gs'},'Right_vs_Left',{'Left' 'Right'},'add')

if isempty(mycol)
    mycol = {'bo' 'gs' 'ro' 'ms' 'yd' 'c^' 'bv' 'ro' 'gs' 'md' 'y^' 'cv' 'bs' 'r^' 'gd' 'mo' 'yv' 'cs'};
end

while length(mycol) < length(newcl), mycol = [mycol mycol];, end
    
for i = 1:length(newcl)
    newcl(i).XYZ = [newcl(i).x newcl(i).y newcl(i).z];
end
    
if length(varargin) < 3, figure('Color','w'); set(gca,'FontSize',18); hold on;, 
else
    subplot(1,3,1); hold on;
end


if length(varargin) == 0
% EACH CLUSTER IS A SEPARATE COLOR

    for i = 1:length(newcl)
        if ~isempty(newcl(i).XYZ)
            plot3(newcl(i).XYZ(:,1),newcl(i).XYZ(:,2),newcl(i).XYZ(:,3),mycol{i},'MarkerFaceColor',mycol{i}(1));
        end
    end

else
    % EACH TASK CONDITION IS A SEPARATE COLOR
    plby = varargin{1}; mlist = varargin{2};
    
    for i = 1:length(mlist)
        
        for j = 1:length(newcl)
            eval(['wh = find(strcmp(newcl(j).' plby ',''' mlist{i} '''));'])
            if ~isempty(wh)
                myh(i) = plot3(newcl(j).XYZ(wh,1),newcl(j).XYZ(wh,2),newcl(j).XYZ(wh,3),mycol{i},'MarkerFaceColor',mycol{i}(1));
            end
        end
        
    end
    
end



if length(varargin) > 2
    
    % ADD to current plot
    
    for splot = 2:3
        subplot(1,3,splot); hold on;
        
    for i = 1:length(mlist)
        
        for j = 1:length(newcl)
            eval(['wh = find(strcmp(newcl(j).' plby ',''' mlist{i} '''));'])
            if ~isempty(wh)
                myh(i) = plot3(newcl(j).XYZ(wh,1),newcl(j).XYZ(wh,2),newcl(j).XYZ(wh,3),mycol{i},'MarkerFaceColor',mycol{i}(1));
            end
        end
        
    end
    
    end
    
    
   
    
else

    % make new stuff
    
addbrain
p = ans;
set(p,'FaceAlpha',.15)
y = [-105 70 70 -105]; z = [-60 -60 78 78]; x = [0 0 0 0];
hold on; hh = fill3(x,y,z,[.9 .9 .9]); set(hh,'EdgeColor','none','FaceColor','w','FaceAlpha',1)

view(90,0)
[az,el] = view;
h = lightangle(az,el); set(gca,'FontSize',18)

[hh1,hh2,hh3,hl,a1,a2,a3] = make_figure_into_orthviews;

if exist('myh') == 1, legend(myh,mlist);,end

set(gcf,'Position',[1727         624        1060         421])

end


return





subplot(1,3,1), hold on

view(0,90)
camzoom(1.25)

subplot(1,3,2), hold on
for i = 1:length(newcl)
    plot3(newcl(i).XYZ(1,:),newcl(i).XYZ(2,:),newcl(i).XYZ(3,:),mycol{i})
end
addbrain
view(270,0)
camzoom(1.25)

subplot(1,3,3), hold on
for i = 1:length(newcl)
    plot3(newcl(i).XYZ(1,:),newcl(i).XYZ(2,:),newcl(i).XYZ(3,:),mycol{i})
end
addbrain
view(0,0)
camzoom(1.25)