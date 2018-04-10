function [classification] = wma_wrapperDev(wbFG,fsDIR)
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
[wbFG, fe] = bsc_LoadAndParseFiberStructure(wbFG);

if  isempty(fe)
    feFlag=false;
else
    feFlag=true;
end
    

% Sets default saving behavior.  Defaults to saving.
if notDefined('nosave'), nosave=false;end

if notDefined('saveHeader'), saveHeader=[];end

%% Segmentation

disp('Segmenting Major and Associative Tracts');

% Segment the major white matter tracts in the Mori Atlas
tic
%[classificationOLD] = wma_majortracts_v2(dt6, wbFG);
classificationOLD=wma_segMajTracks_v4(wbFG, fsDIR);
classification=classificationOLD;
segmentTime=toc;
fprintf ('\n Initial major tracts complete in %4.2f hours \n', segmentTime/(60*60))
%erase cingulum



% update name field
newTractNames={'Left pArc','Right pArc','Left TPC','Right TPC','Left MdLF-SPL','Right MdLF-SPL', 'Left MdLF-Ang','Right MdLF-Ang', 'Left SLF1', 'Right SLF1', 'Left SLF2', 'Right SLF2', 'Left SLF3', 'Right SLF3', 'Left ILF', 'Right ILF','Left Aslant','Right Aslant'};
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

% Middle Longitudinal Fasiculus segmentation

[~, RightILFIndexes, ~, LeftILFIndexes, ~, RightMdLFsplIndexes, ~, LeftMdLFsplIndexes,... 
    ~, RightMdLFangIndexes, ~, LeftMdLFangIndexes] =bsc_segmentMdLF_ILF_v2(wbFG, fsDIR);

classification.index(LeftILFIndexes)=find( strcmp(classification.names,'Left ILF'));
classification.index(RightILFIndexes)=find( strcmp(classification.names,'Right ILF'));
classification.index(LeftMdLFangIndexes)=find( strcmp(classification.names,'Left MdLF-Ang'));
classification.index(RightMdLFangIndexes)=find( strcmp(classification.names,'Right MdLF-Ang'));
classification.index(LeftMdLFsplIndexes)=find( strcmp(classification.names,'Left MdLF-SPL'));
classification.index(RightMdLFsplIndexes)=find( strcmp(classification.names,'Right MdLF-SPL'));



%[~, RightMeyerBool, ~,RightBaumBool, ~, LeftMeyerBool, ~,LeftBaumBool] =bsc_opticRadiationSeg_V3(wbfg, fsDir)
%functionally the same if outputs are indexes or bools, I believe

[classificationOptic] =bsc_opticRadiationSeg_V4(wbFG, fsDIR);


[~, LeftSLF1Bool, ~, RightSLF1Bool, ~, LeftSLF2Bool, ~, RightSLF2Bool,...
    ~, LeftSLF3Bool, ~, RightSLF3Bool] =wma_subsegSLF_V3(wbFG, fsDIR);


classification.index(LeftSLF1Bool)=find( strcmp(classification.names,'Left SLF1'));
classification.index(RightSLF1Bool)=find( strcmp(classification.names,'Right SLF1'));
classification.index(LeftSLF2Bool)=find( strcmp(classification.names,'Left SLF2'));
classification.index(RightSLF2Bool)=find( strcmp(classification.names,'Right SLF2'));
classification.index(LeftSLF3Bool)=find( strcmp(classification.names,'Left SLF3'));
classification.index(RightSLF3Bool)=find( strcmp(classification.names,'Right SLF3'));


[~, RightAslantIndexes, ~, LeftAslantIndexes] =bsc_segmentAslantProto(wbFG, fsDIR);

classification.index(LeftAslantIndexes)=find( strcmp(classification.names,'Left Aslant'));
classification.index(RightAslantIndexes)=find( strcmp(classification.names,'Right Aslant'));


[classificationAdd] =bsc_cerebellarTracts(wbFG, fsDIR);

[classificationSemiFinal] = bsc_mergeClassifications(classification,classificationAdd);

[classification] = bsc_mergeClassifications(classificationSemiFinal,classificationOptic);


%[classification] = bsc_deleteClassifications(classificationFinal,[34 35 56 57 58 59 60 61 62 63 64 65 66]);

if feFlag
classification=wma_clearNonvalidClassifications(classification,fe);
end

classification=removeOutliersClassificationDev(classification,wbFG, 4, 10);

disp('\n Tracts segmentation complete');

return
