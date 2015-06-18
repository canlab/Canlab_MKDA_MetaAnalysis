A = MVAR.MCA.class_avg_indicator(:, 6);
B =  MVAR.MCA.class_avg_indicator(:, 9);

%%
clear z

% Make table z of counts.
% Note: A and B should be 0 or 1
% table in same format as Patel 2006 paper
% A's activity = columns, active then inactive
% B's activity = rows, active then inactive

n = size(A,1);      % number of observations

if n ~= size(B,1), error('A and B must be same length!'); end
w = ones(n, 1);     % weights (not used now)

u1 = sort(unique(B), 'descend');
u2 = sort(unique(A), 'descend');
for i = 1:length(u1)                   % rows are 0 then 1 on the first var, "Nos" then "Yesses"
    for j = 1:length(u2)               % for each column
        z(i,j) = sum( (B == u1(i) & A == u2(j)) .* w );
    end
end

%Z is a table of frequency of activation/inactivation for the two regions  
% P is probabilities of being in each of the 4 cats, based on observed data
% P is same as theta in the paper
% each obs. can be in only one of the 4 categories, which are mutually exclusive
P = z ./ n; %probability table

% g is the elements of theta
g1 = P(1,1);  % Act , Act
g2 = P(2,1);  % Act, Inact
g3 = P(1,2);  % Inact, Act
g4 = P(2,2);  % Inact, Inact

%%
% Calculate Kappa, measure of association

% p(R1, R2) >= p(R1) * P(R2)    Is measured joint prob. of R1, R2 active
% greater than expected given independence?
% dir == 1 signifies positive dependence.  dir == 0 sigifies zero or
% negative dependence
dir = double(g1 >= (g1+g2).*(g1+g3));

%dir = double(g1 >= (g1+g2).*(g1+g3));


W = dir.*(0.5*(g1-((g1+g2).*(g1+g3)))./(g1+min(g2,g3)-(g1+g2).*(g1+g3))+0.5) + (1-dir).*(.5-0.5*(-g1+(g1+g2).*(g1+g3))./((g1+g2).*(g1+g3) - max(0,2*g1+g2+g3-1)));

% measure of association, between -1 and 1
% g is theta in Patel 2006
% max(theta1) in paper is
kappa = (g1-((g1+g2).*(g1+g3)))./(W.*(g1+ min(g2,g3)-(g1+g2).*(g1+g3)) + (1-W).*((g1+g2).*(g1+g3) - max(0,2*g1+g2+g3-1)));


%% Calculate posterior probabilities

theta = dirrnd(gamma,n);
g1 = theta(:,1);
g2 = theta(:,2);
g3 = theta(:,3);


% direction of effect; now 1 if direction of effect is "R1 ascendant to R2"
dir2 = double((g1+g2)>=(g1+g3));  

% num: either 0 or 'ascendancy' ratio; 1 - (p(R2) / p(R1)
% g1 + g3 = p(Region2 is active)
% g1 + g2 = p(Region1 is active)

% Ascendant regions are 'supersets'
% if A is ascendant, ratio P(A|B) / P(B|A) is greater than 1
tau = dir2.*(1 - (g1+g3)./(g1+g2)) + (1-dir2).*((g1+g2)./(g1+g3)-1);


% % % function out = dirrnd(a,n)
% % % out = gamrnd(ones(n,1)*a,ones(n,length(a)));
% % % total = sum(out')'*ones(1,length(a));
% % % out = out./total;

%%

% % % % This program takes in binary voxel data (active or not) from subjects
% % % % Each subject is in a file of s*vc where vx is the number of voxels in the
% % % % roi and s is the number of scans
% % % % FLAT PRIORS FOR NOW, thus prior sample size = 0;
% % % % including prior sample size
% % % % gamma*(i,j) are only non-zero when i<j and i and j
% % % % FIRST BE SURE WE ARE READY TO GO WITHIN FC.M
% % % function [gamma1,gamma2,gamma3,khigh,scans,data] = getprobs(root,tempdir,kthresh,v,subjects,vc,data,restrict)
% % % 
% % % if (isCell(data) == 0)
% % % PB1 = progressbar([],0,'Loading Data into Memory','Progress Indicator');
% % %     data = cell(size(subjects,1),1);
% % %     for i = 1:size(subjects,1)        
% % %         data{i} = load(strcat(root,tempdir,subjects(i,:)));
% % %         data{i}.retval = data{i}.retval(restrict{i},:);
% % %         progressbar(PB1,1/size(subjects,1));
% % %     end
% % % progressbar(PB1,-1);
% % % end
% % % 
% % % % Store as logicals and convert when using for multiple voxels
% % % SumV1V2 = zeros(1,vc);
% % % vgamma1 = zeros(1,vc);
% % % vgamma2 = zeros(1,vc);
% % % vgamma3 = zeros(1,vc);
% % % pass = 0;
% % % n = 0;
% % % A = ones(1,vc);
% % % % Loop here over subjects
% % % PB = progressbar([],0,'Determining significant functional connectivity...','Progress Indicator');
% % % for(i = 1:size(subjects,1))    
% % %     n = n + size(data{i}.retval,1);
% % %     SumV = sum(data{i}.retval);
% % %     A(find(SumV == 0)) = 0;
% % %     d2 = repmat(data{i}.retval(:,v),1,size(data{i}.retval,2));
% % %     SumV1V2 = sum(d2 & data{i}.retval);
% % %     vgamma1 = vgamma1 + SumV1V2;
% % %     vgamma2 = vgamma2 - SumV1V2 + SumV(v);
% % %     vgamma3 = vgamma3 - SumV1V2 + SumV;    
% % %     progressbar(PB,1/size(subjects,1));
% % % end
% % % 
% % % vgamma1 = vgamma1 + 1;
% % % vgamma2 = vgamma2 + 1;
% % % vgamma3 = vgamma3 + 1;
% % % n = n+4;
% % % 
% % % g1 = vgamma1 / n;
% % % g2 = vgamma2 / n;
% % % g3 = vgamma3 / n;
% % % g4 = (1 - g1 - g2 - g3);
% % % dir = double(g1 >= (g1+g2).*(g1+g3));
% % % W = dir.*(0.5*(g1-((g1+g2).*(g1+g3)))./(g1+min(g2,g3)-(g1+g2).*(g1+g3))+0.5) + (1-dir).*(.5-0.5*(-g1+(g1+g2).*(g1+g3))./((g1+g2).*(g1+g3) - max(0,2*g1+g2+g3-1)));
% % % 
% % % % measure of association, between -1 and 1
% % % % g is theta in Patel 2006
% % % % max(theta1) in paper is 
% % % kappa = (g1-((g1+g2).*(g1+g3)))./(W.*(g1+ min(g2,g3)-(g1+g2).*(g1+g3)) + (1-W).*((g1+g2).*(g1+g3) - max(0,2*g1+g2+g3-1)));
% % % 
% % % 
% % % % Threshold here
% % % khigh = find(abs(kappa) > kthresh & A);
% % % gamma1 = vgamma1(khigh);
% % % gamma2 = vgamma2(khigh);
% % % gamma3 = vgamma3(khigh);
% % % scans = n;
% % % progressbar(PB,-1);
% % % end




% % 
% % function [sigf, sigfn, sige, sigen] = postprobs(gamma,effect,n);
% % theta = dirrnd(gamma,n);
% % g1 = theta(:,1);
% % g2 = theta(:,2);
% % g3 = theta(:,3);
% % 
% % dir = double(g1 >= ((g1+g2).*(g1+g3)));
% % W = dir.*(0.5*(g1-((g1+g2).*(g1+g3)))./(g1+min(g2,g3)-(g1+g2).*(g1+g3))+0.5) + (1-dir).*(.5-0.5*(-g1+(g1+g2).*(g1+g3))./((g1+g2).*(g1+g3) - max(0,2*g1+g2+g3-1)));
% % kappa = (g1-((g1+g2).*(g1+g3)))./(W.*(g1+min(g2,g3)-(g1+g2).*(g1+g3)) + (1-W).*((g1+g2).*(g1+g3) - max(0,2*g1+g2+g3-1)));
% % dir2 = double((g1+g2)>=(g1+g3));
% % tau = dir2.*(1 - (g1+g3)./(g1+g2))+(1-dir2).*((g1+g2)./(g1+g3)-1);
% % sigf = sum(kappa > effect)./n;
% % sigfn = sum(kappa < -1*effect)./n;
% % sige = sum(tau > effect)./n;
% % sigen = sum(tau < -1*effect)./n;

