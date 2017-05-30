function bsc_plot_AFQ_PostNetwork_saggital_Horizontalv2_MDLF(fullFiberOutDir, fullFigureOutDir, t1path, saveHeaderData,posIndexes,wbFG)
% This function plots the arcuate, ilf, ifof and slf and mdlf and saves a saggital
% figure
%
% INPUTS
% fullFiberOutDir:  path to directory contining the relevant VPF, Arc pArc
% and VOF fibers which were generated earlier (i.e. VPFandOtherFiberSegmentWrapperHCP1.m)
% this path is a standard output from bsc_HCP_Stanford_PathGen
%
% fullFigureOutDir: path to directory where output figures will be saved
% this path is a standard output from bsc_HCP_Stanford_PathGen
%
% t1path: path to t1 file
% this path is a standard output from bsc_HCP_Stanford_PathGen
%
% OUTPUTS
%  no outputs, all figures are saved into the fullFigureOutDir
%
%
% %% create parpool compatible with qsub
%
% % create parallel cluster object
% c = parcluster;
%
% % set number of cores from arguments
% c.NumWorkers = 16;
%
% % set temporary cache directory
% t = tempname('/gpfs/home/d/n/dnbulloc/Karst/qsubParpools');
%
% % make cache dir
% OK = mkdir(t);
%
% % check and set cachedir location
% if OK
%     % set local storage for parpool
%     c.JobStorageLocation = t;
% end
%
% % start parpool
% parpool(c, 16, 'IdleTimeout', 720);

%% load up output and put into structure
%reminder
%key
%1 = 'L_Arcuate_Posterior'
%2 = 'R_Arcuate_Posterior'
%3 = 'L_VOF'
%4 = 'R_VOF'
%5 = 'Left Arcuate'
%6 = 'Right Arcuate'
%7 = 'LeftLatVPF'
%8 = 'LeftMedVPF'
%9 = 'RightLatVPF'
%10 ='RightMedVPF'
%11 = 'Left MdLF'
%12 = 'Right MdLF'
%13 = 'Left Thalamic Radiation'
%14 = 'Right Thalamic Radiation'
%15 = 'Left Corticospinal'
%16 = 'Right Corticospinal'
%17 = 'Left Cingulum Cingulate'
%18 = 'Right Cingulum Cingulate'
%19 = 'Left Cingulum Hippocampus'
%20 = 'Right Cingulum Hippocampus'
%21 = 'Callosum Forceps Major'
%22 = 'Callosum Forceps Minor'
%23 = 'Left IFOF'
%24 = 'Right IFOF'
%25 = 'Left ILF'
%26 = 'Right ILF'
%27 = 'Left SLF'
%28 = 'Right SLF'
%29 = 'Left Uncinate'
%30 = 'Right Uncinate'

%desired order
% [ 5 6 27 28 17 18 15 16 21 22 13 14 29 30 19 20 23 24 25 26 11 12 3 4 7 9
% 8 10]









load(strcat(fullFiberOutDir,saveHeaderData,'_FiberIndexes.mat'));


pos_L_Arc_Indexes=FiberIndexes.L_Arc(ismember(FiberIndexes.L_Arc,posIndexes));
pos_R_Arc_Indexes=FiberIndexes.R_Arc(ismember(FiberIndexes.R_Arc,posIndexes));
pos_L_SLF_Indexes=FiberIndexes.L_SLF(ismember(FiberIndexes.L_SLF,posIndexes));
pos_R_SLF_Indexes=FiberIndexes.R_SLF(ismember(FiberIndexes.R_SLF,posIndexes));
pos_L_ILF_Indexes=FiberIndexes.L_ILF(ismember(FiberIndexes.L_ILF,posIndexes));
pos_R_ILF_Indexes=FiberIndexes.R_ILF(ismember(FiberIndexes.R_ILF,posIndexes));
pos_L_IFOF_Indexes=FiberIndexes.L_IFOF(ismember(FiberIndexes.L_IFOF,posIndexes));
pos_R_IFOF_Indexes=FiberIndexes.R_IFOF(ismember(FiberIndexes.R_IFOF,posIndexes));
%pos_R_MdLF_Indexes=FiberIndexes.R_MdLF(ismember(FiberIndexes.R_MdLF,posIndexes));
%pos_L_MdLF_Indexes=FiberIndexes.L_MdLF(ismember(FiberIndexes.L_MdLF,posIndexes));

% Fibers.L_VOF.fibers=Fibers.L_VOF.fibers(ismember(FiberIndexes.L_VOF,posIndexes));
% Fibers.R_VOF.fibers=Fibers.R_VOF.fibers(ismember(FiberIndexes.R_VOF,posIndexes));
%
% Fibers.L_Arc.fibers=Fibers.L_Arc.fibers(ismember(FiberIndexes.L_Arc,posIndexes));
% Fibers.R_Arc.fibers=Fibers.R_Arc.fibers(ismember(FiberIndexes.R_Arc,posIndexes));
%
% Fibers.L_Lat_VPF.fibers=Fibers.L_Lat_VPF.fibers(ismember(FiberIndexes.L_Lat_VPF,posIndexes));
% Fibers.R_Lat_VPF.fibers=Fibers.R_Lat_VPF.fibers(ismember(FiberIndexes.R_Lat_VPF,posIndexes));
% Fibers.L_Med_VPF.fibers=Fibers.L_Med_VPF.fibers(ismember(FiberIndexes.L_Med_VPF,posIndexes));
% Fibers.R_Med_VPF.fibers=Fibers.R_Med_VPF.fibers(ismember(FiberIndexes.R_Med_VPF,posIndexes));
%
% Fibers.L_VOF=bsc_SegmentVOFfromAFQout(Fibers.L_VOF);
% Fibers.R_VOF=bsc_SegmentVOFfromAFQout(Fibers.R_VOF);

fiberclear=wbFG;
fiberclear.fibers=[];
%stupid matlab data structure issue
fgs=[];
for ifiberstruc=1:8
    fgs{ifiberstruc}=fiberclear;
end
%fgs{1}.fibers=wbFG.fibers(pos_L_VOF_Indexes);
fgs{1}.fibers=wbFG.fibers(pos_L_Arc_Indexes);
fgs{2}.fibers=wbFG.fibers(pos_R_Arc_Indexes);
%fgs{1}=bsc_SegmentVOFfromAFQout(fgs{1});

%left Arcuate
%fgs{3}.fibers=wbFG.fibers(pos_L_Arc_Indexes);

%right Arcuate
fgs{3}.fibers=wbFG.fibers(pos_L_SLF_Indexes);
fgs{4}.fibers=wbFG.fibers(pos_R_SLF_Indexes);

fgs{5}.fibers=wbFG.fibers(pos_L_ILF_Indexes);
fgs{6}.fibers=wbFG.fibers(pos_R_ILF_Indexes);

fgs{7}.fibers=wbFG.fibers(pos_L_IFOF_Indexes);
fgs{8}.fibers=wbFG.fibers(pos_R_IFOF_Indexes);






% fibers{7}=LfgMed;
% fibers{5}=LfgLat;
% fibers{6}=RfgLat;
% fibers{8}=RfgMed;




%% plotting stuff

% notes on standard colors
%fiber1=Left VOF [.88,0,1]
%fiber2=Right VOF [.88,0,1]
%fiber3=Arcuate 1 [1,0,.3]
%fiber4=Arcuate 2 [1,0,.3]
%fiber5=left med VPF [.6,1,0]
%fiber6=Right Med VPF [.6,1,0]
%fiber7=Left Lat VPF [1,1,0]
%fiber8=Right lat VPF [1,1,0]
%fiber9=Left pArc [.95,.6,.05]
%fiber10=Right pArc  [.95,.6,.05]

%flips fibers if the cells are arranged 1 by N, and performs outlier
%remover

%define cut deviations for later pruning algorithm

cutvar1=[3 3 3 3 3 3 3 3 3 3 3 3];
cutvar2=[3 3 3 3 3 3 3 3 3 3 3 3];

for iFiberGroups=1:length(fgs)
    if ~isempty(fgs{iFiberGroups})
        if ~isempty(fgs{iFiberGroups}.fibers)
            fiberDim=size(fgs{iFiberGroups}.fibers);
            if fiberDim(1)<fiberDim(2)
                fgs{iFiberGroups}.fibers=fgs{iFiberGroups}.fibers';
            end
            [~, keep] = mbaComputeFibersOutliers(fgs{iFiberGroups},cutvar1(iFiberGroups),cutvar2(iFiberGroups));
            fprintf('\n Found a tract with %i fibers in %s... \n',sum(keep),fgs{iFiberGroups}.name);
            fgs{iFiberGroups}.params=[];
            fgs{iFiberGroups} = fgExtract(fgs{iFiberGroups},find(keep),'keep');
        end
    end
end

fiberColors={[1.0000 0.1034 0.7241],[1.0000 0.1034 0.72417],[1.0000 0.8276 0],[1.0000 0.8276 0],[0 0.5172 0.9655],[0 0.5172 0.9655],[0 0.5517 0.2759],[0 0.5517 0.2759]};

t1          = niftiRead(t1path);
slices      = {[-1 0 0],[0 -45 0],[0 0 -1]};

fh = figure('name','HorizontalFibers','color','k','units','normalized','position',[.5 .5 .5 .5],'visible','off');
axis square
fhNum = fh.Number;
hold on

h  = mbaDisplayBrainSlice(t1, slices{1});
%h  = mbaDisplayBrainSlice(t1, slices{2});
%h  = mbaDisplayBrainSlice(t1, slices{3});

ylim([-120,60]);
zlim([-40,95]);

for itract = 1:length(fgs)
    if exist('lh','var'), delete(lh); end
    if ~isempty(fgs{itract})
        if ~isempty(fgs{itract}.fibers)
            [fh, lh] = mbaDisplayConnectome(fgs{itract}.fibers,fhNum, fiberColors{itract}, 'single');%color{itract}
            
            delete(lh)
            display (itract)
            fprintf('\n %i \n',itract)
        end
    end
end

%find Variant Name (i.e. the name of this batch of fiber/plot outputs) and
%infers subject ID
slashIndicies=strfind(fullFiberOutDir,'/');
BatchName=fullFiberOutDir(slashIndicies(end)+1:end);
SubjectID=fullFiberOutDir(slashIndicies(end-2)+1:slashIndicies(end-1)-1);




fig.views = {[90,0],[-90,0]};
light.angle = {[90,-45],[-90,45]};
fig.names = {strcat(saveHeaderData, '_PosteriorNet_saggital1_Horizontal'),strcat(saveHeaderData, '_PosteriorNet_saggital2_Horizontal')};

for iview = 1:length(fig.views)
    view(fig.views{iview}(1),fig.views{iview}(2))
    lh = camlight('left');
    
    
    feSavefig(fhNum,'verbose','yes', ...
        'figName',fig.names{iview}, ...
        'figDir',fullFigureOutDir, ...
        'figType','jpg');
    delete(lh)
end

clear fibers
clear LeftMed
clear LeftLat;
clear RightLat;
clear RightMed;
close all


poolobj = gcp('nocreate');
delete(poolobj);

end

