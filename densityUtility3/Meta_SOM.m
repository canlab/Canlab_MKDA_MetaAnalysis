function [cl,anyStudy,studyByCluster,SOMResults,theData,whSOM] = Meta_SOM(P,nGrid,theData)
% [cl,anyStudy,studyByCluster,SOMResults,theData,whSOM] = Meta_SOM(P,nGrid,[theData],[mask name])
%
% test for meta-analysis data.
% tor wager modified this from Robert Welsh's som_test_01.m
%
% [cl,anyStudy,studyByCluster,SOMResults,theData,whSOM] = Meta_SOM(PP,5);

whSOM = [];
maskInfo = [];

% go to subdirectory
[dd] = pwd;
if strcmp(dd(end-10:end),'SOM_Results')
    % do nothing
else
    fprintf('Going to directory: SOM_Results\n');
    mkdir SOM_Results
    cd('SOM_Results');
end


% Are the data broken into sessions.
nSessions  = 1;

% Create the Mask on where to actually analyze data. Writes out mask.img
%SOM_CreateMetaMask(P);     % only where ALL voxels are nonzero
mname = 'mask.img';
if exist(mname) == 2
    % do nothing
else
    maskthresh = .03;
    tor_spm_mean_ui(P,mname); % for meta: where ANY vox are nonzero
    V=spm_vol(mname); v = spm_read_vols(V);
    v(v < maskthresh) = 0;
    spm_write_vol(V,v);
end

if nargin < 4 || isempty(maskName)
    maskName = fullfile(pwd,'mask.img');
    %if nargin < 2 || isempty(maskName), maskName = spm_get(1,'*IMAGE','Select mask.');, end

    if nargin < 3 || isempty(theData)
        % Read only the masked data. Returns an array of theData(nVoxels,nImages).
        [theData maskInfo] = SOM_PrepData(P,maskName);
    else
        % data already loaded
    end
end

save SOMResults maskInfo P
% normalize data

% not implemented now.




fprintf('Starting SOM Calculation\n');


cp = cputime;
tic;


% How long to iterate for.
nIter      = 50;

% Assuming a square NET, the length of the side.
%

%nGrid = max(5,round(sqrt(size(P,1) ./10)));
disp(['nGrid: ' num2str(nGrid)]);



% Now prepare to actually calculate the SOM.

nSOM = nGrid^2;
disp(['nSOM, number of canonical study activation profiles: ' num2str(nSom)]);


% Calculate the SOM

SOMResults        = SOM_CalculateMap(theData,nSOM,nIter);

% Store the header information - needed for writing out results as images.

SOMResults.header = maskInfo.header;
SOMResults.iMask  = maskInfo.iMask;

% Organize the data into super clusters.

[SOMResults.SuperCluster SOMResults.nCluster] = SOM_SuperCluster(SOMResults.SOM);

% Final amount of time.
toc
cputime - cp

save SOMResults -append SOMResults

fprintf('Writing results images\n');


% now make results masks
[wts,indices] = SOM_WriteIMGS(SOMResults);
wts = sort(wts',1,'descend');
tor_fig; plot(wts,'k','LineWidth',3); title('Map weights');
xlabel('Map')

%SOM_ViewMap


fprintf('Getting results clusters\n');

%[S,H] = silhouette(theData, SOMResults.IDX,'Euclidean');

[cl,anyStudy,studyByCluster] = meta_SOMclusters(SOMResults,theData,[],'group');

save SOMResults -append cl anyStudy studyByCluster

return




