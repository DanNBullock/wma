function [summaryStructure] = bsc_multiTractAnalysis(tractNameList,fiberIndexList,feORwbFG,FiberDir,saveHeader,dt6,fsDIR)
%[summaryStructure] = bsc_multiTractAnalysis(tractNameList,fiberIndexList,feORwbFG,FiberDir,saveHeader,dt6,fsDIR)
%
% OVERVIEW: this function runs a bevy of analysis functions including 
% bsc_multiVirtualLesion, bsc_multiTractVolumeFE, bsc_multiTractLengthFG,
% and bsc_multiTractProfilesFE.  It combines the output from several of
% these into a single summary structure and saves it down.
%
% INPUTS:
% -tractNameList: a list of tract names.  Must be same lenght and same
% order as the fiberIndexList variable
%
% -fiberIndexList: a list of the same length as the tractNameList.  Each
% entry in this list is itself a list of indexes into the input fe
% structure, such that, for each entry, the indexed streamlines compose
% the named tract.
%
% -feORwbFG: either a path to a saved wbFG or a wbFG object.  If you pass a fe
% structure or an fe path it will (probably) extract the fg from it
% appropriately and subsequently used the positively weighted tracts when relevant.
%
% -FiberDir: directory path for the directory you would like your fiber
% indexes saved down to.
%
% -saveHeader: a string corresponding to whatever iformation you would like
% to serve as the file identifier for the fiber tract object.  Could
% include information like tractography parameter settings, subject
% identifier or group ID.  (example: 'HCP_105115_PROB_Lmax2_conn5')
%
% -dt6: either a path to a dt6 file or a dt6 object
%
% -fsDIR: path  to THIS SUBJECT'S freesurfer directory
%
% OUTPUTS:
% summaryStructure: an amalgamated data structure containing the outputs
% from bsc_multiVirtualLesion, bsc_multiTractLengthFG, and
% bsc_multiTractProfilesFE.  See those respective functions for more
% details.
%
%% preliminaries


% loads whatver was passed to it with wbFG
if ischar(feORwbFG)
    feORwbFG = load(feORwbFG);
else
end

% loads the dt6 file if a path was passed
if ischar(dt6)
    dt6 = dtiLoadDt6(dt6);
else
end

%if it is a fe structure, get the wbFG out of it
if isfield(feORwbFG, 'fg')
    positiveWeightIndexes=find(~feORwbFG.life.fit.weights==0);
    for iFieldName=1:length(tractNameList)
        posFiberIndexList{iFieldName}=intersect(positiveWeightIndexes,fiberIndexList{iFieldName});
    end
else
    fprintf('\n fe structure not detected, will not run virtual lesion \n')
end

%% virtual lesion
if exist('positiveWeightIndexes', 'var')
    [FullVL] = bsc_multiVirtualLesion(tractNameList,fiberIndexList,feORwbFG,FiberDir,saveHeader,0);
end

    %% endpoint Volume creation
    
    %can sub index so that you only creat volumes for tracts of interest
    %desired fiber indexes in the tractNameList object:
    subSelect=[23:28];
    subSelectBool(1:length(tractNameList))=false;
    subSelectBool(subSelect)=true;
    
    bsc_multiTractVolumeFE(tractNameList(subSelectBool),posFiberIndexList(subSelectBool),feORwbFG,FiberDir,saveHeader, fsDIR);
    
       
    %% Quantative Stats
    [LengthStruc] =bsc_multiTractLengthFG(tractNameList,posFiberIndexList,feORwbFG,FiberDir,saveHeader,0);
    
    %% Tract Profiles
    %doesn't work right now?
    %     [DiffusionProfile]=bsc_multiTractProfilesFE(tractNameList,posFiberIndexList,feORwbFG,FiberDir,saveHeader,dt6, 0);
    
    %% creating the amalgum
    if notDefined('DiffusionProfile'), DiffusionProfile=[];end
    summaryStructure=[DiffusionProfile LengthStruc FullVL];
    
    save (strcat(fullfile(FiberDir),saveHeader,'_summaryStructure.mat'),'summaryStructure','-v7.3');
end