function [classification] = wma_wrapper(wbFG,dt6,fsDIR)
% 
% [tracts_indexes]=wma_wrapper(wbFG,dt6path,FiberDir,saveHeader)
%
% OVERVIEW:
% Segments out fibertracts from the AFQ package, the VOF, and several other
% tracts (i.e. MdLF, pArc, and TPC)
%
% INPUTS:
% -wbFG: either a path to a saved wbFG or a wbFG object.  If you pass a fe
% structure or an fe path it will (probably) extract the fg from it.
% appropriately.
%
% -dt6: either a path to a saved dt6 file or a dt6 object
%
% -fsDIR: path  to THIS SUBJECT'S freesurfer directory
%
% OUTPUTS:
% -classification:  A structure whose field names correspond to the names
%  of the various fiber tracts that were segmented in this function.  Each
%  field stores the indexes into the wbFG corresponding to the named track.
%  i.e. tracts_indexes.R_pArc will correspond to several hundered (probably)
%  indexes into the wbFG for the right posterior arcuate.
%
%  NOTE: This function does not save down MERELY validated tracts (i.e. via
%  LiFE), nor does it run an outlier removal function (i.e
%  mbaComputeFibersOutliers).  The output indexes are "unaltered" in this
%  regard.  The only exception to this is the VOF which has had outlier
%  removal at 2 for gradient (i.e. orientation)
%
% (C) Daniel Bullock, 2017, Indiana University

%% Path generation and initialization

% loads file if a string was passed 
if ischar(wbFG)
    wbFG = load(wbFG);
    %if it is a fe structure, get the wbFG out of it
    if isfield(wbFG, 'fe')
        wbFG = feGet(wbFG.fe, 'fibers acpc');
    end
else
    if isfield(wbFG, 'fg')
        wbFG = feGet(wbFG, 'fibers acpc');
    end
end

if ischar(dt6)
    dt6 = dtiLoadDt6(dt6);
else

end

% Sets default saving behavior.  Defaults to saving.
if notDefined('nosave'), nosave=false;end

if notDefined('saveHeader'), saveHeader=[];end

%% Segmentation

disp('Segmenting Major and Associative Tracts');

% Segment the major white matter tracts in the Mori Atlas
tic
[classificationOLD] = wma_majortracts_v4(dt6, wbFG);
segmentTime=toc;
fprintf ('\n Mori Atlas segmentation run complete in %4.2f hours \n', segmentTime/(60*60))

classification=classificationOLD;

% update name field
newTractNames={'Left VOF','Right VOF','Left pArc','Right pArc','Left TPC','Right TPC','Left MdLF-SPL','Right MdLF-SPL', 'Left MdLF-Ang','Right MdLF-Ang','Left Meyer', 'Right Meyer', 'Left Baum', 'Right Baum', 'Left SLF1', 'Right SLF1', 'Left SLF2', 'Right SLF2', 'Left SLF3', 'Right SLF3', 'Left ILF', 'Right ILF' };
for iNEWtracts=length(classification.names)+1:length(classification.names)+length(newTractNames)

    
   classification.names{iNEWtracts}=newTractNames{iNEWtracts-length(classificationOLD.names)};
end

% posterior arcuate and temporo-parietal connection segmentation
[L_pArc, ~, L_pArc_Indexes,L_TPC_Indexes ,R_pArc, ~, R_pArc_Indexes,R_TPC_Indexes] = bsc_automated_roi_segment_script_neo(wbFG,fsDIR);

classification.index(L_pArc_Indexes)=find( strcmp(classification.names,'Left pArc'));
classification.index(R_pArc_Indexes)=find( strcmp(classification.names,'Right pArc'));
classification.index(L_TPC_Indexes)=find( strcmp(classification.names,'Left TPC'));
classification.index(R_TPC_Indexes)=find( strcmp(classification.names,'Right TPC'));

% Segment the Vertical Occipital Fasiculus (VOF)
[~, ~, L_VOF_Indexes, R_VOF_Indexes] =  wma_segment_vof(wbFG, fsDIR, classification, dt6); 

classification.index(L_VOF_Indexes)=find( strcmp(classification.names,'Left VOF'));
classification.index(R_VOF_Indexes)=find( strcmp(classification.names,'Right VOF'));

% Middle Longitudinal Fasiculus segmentation

[~, RightILFIndexes, ~, LeftILFIndexes, ~, RightMdLFsplIndexes, ~, LeftMdLFsplIndexes,... 
    ~, RightMdLFangIndexes, ~, LeftMdLFangIndexes] =bsc_segmentMdLF_ILF(wbFG, fsDIR);

classification.index(LeftILFIndexes)=find( strcmp(classification.names,'Left ILF'));
classification.index(RightILFIndexes)=find( strcmp(classification.names,'Right ILF'));
classification.index(LeftMdLFangIndexes)=find( strcmp(classification.names,'Left MdLF-Ang'));
classification.index(RightMdLFangIndexes)=find( strcmp(classification.names,'Right MdLF-Ang'));
classification.index(LeftMdLFsplIndexes)=find( strcmp(classification.names,'Left MdLF-SPL'));
classification.index(RightMdLFsplIndexes)=find( strcmp(classification.names,'Right MdLF-SPL'));



%[~, RightMeyerBool, ~,RightBaumBool, ~, LeftMeyerBool, ~,LeftBaumBool] =bsc_opticRadiationSeg_V3(wbfg, fsDir)
%functionally the same if outputs are indexes or bools, I believe
[~, RightMeyerInd, ~,RightBaumInd, ~, LeftMeyerInd, ~,LeftBaumInd] =bsc_opticRadiationSeg_V3(wbFG, fsDIR);

classification.index(LeftMeyerInd)=find( strcmp(classification.names,'Left Meyer'));
classification.index(RightMeyerInd)=find( strcmp(classification.names,'Right Meyer'));
classification.index(LeftBaumInd)=find( strcmp(classification.names,'Left Baum'));
classification.index(RightBaumInd)=find( strcmp(classification.names,'Right Baum'));

[~, LeftSlf1Bool, ~, RightSlf1Bool, ~, LeftSlf2Bool, ~, RightSlf2Bool,...
    ~, LeftSlf3Bool, ~, RightSlf3Bool] =wma_subsegSLF(wbFG, fsDIR);


classification.index(LeftSlf1Bool)=find( strcmp(classification.names,'Left SLF1'));
classification.index(RightSlf1Bool)=find( strcmp(classification.names,'Right SLF1'));
classification.index(LeftSlf2Bool)=find( strcmp(classification.names,'Left SLF2'));
classification.index(RightSlf2Bool)=find( strcmp(classification.names,'Right SLF2'));
classification.index(LeftSlf3Bool)=find( strcmp(classification.names,'Left SLF3'));
classification.index(RightSlf3Bool)=find( strcmp(classification.names,'Right SLF3'));


disp('\n Tracts segmentation complete');

return
