load surf_single_subj_T1
figure; p = patch('Faces',faces,'Vertices',vertices)
view(135,30)
set(p,'EdgeColor','none')
lighting gouraud;
alpha(.5); lighth = lightangle(45,30); axis off;
set(p,'FaceColor',[.2 .4 .7])
% camh = camlight; 
daspect([1 1 1]);
%set(gcf,'Color',[0 0 0])
hold on

%a = ones(200,150);
%a = a + 75;
%ah = surf(a)
%set(ah,'EdgeColor','none')
%set(ah,'FaceColor','red')
%alpha(ah,.3)