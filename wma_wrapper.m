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
[classification] = wma_majortracts(dt6, wbFG);
segmentTime=toc;
fprintf ('\n Mori Atlas segmentation run complete in %4.2f hours \n', segmentTime/(60*60))

% update name field
newTractNames={'Left VOF','Right VOF','Left pArc','Right pArc','Left TPC','Right TPC','Left MdLF','Right MdLF', 'Left Meyer', 'Right Meyer', 'Left Baum', 'Right Baum' };
for iNEWtracts=21:32
   classification.names{iNEWtracts}=newTractNames{iNEWtracts-20};
end

% posterior arcuate and temporo-parietal connection segmentation
[L_pArc, ~, L_pArc_Indexes,L_TPC_Indexes ,R_pArc, ~, R_pArc_Indexes,R_TPC_Indexes] = bsc_automated_roi_segment_script_neo(wbFG,fsDIR);

classification.index(L_pArc_Indexes)=23;
classification.index(R_pArc_Indexes)=24;
classification.index(L_TPC_Indexes)=25;
classification.index(R_TPC_Indexes)=26;

% Segment the Vertical Occipital Fasiculus (VOF)
[~, ~, L_VOF_Indexes, R_VOF_Indexes] =  wma_segment_vof(wbFG, fsDIR, classification, dt6); 

classification.index(L_VOF_Indexes)=21;
classification.index(R_VOF_Indexes)=22;

% Middle Longitudinal Fasiculus segmentation
[~, RightMdLFindexes, ...
 ~, LeftMdLFindexes] = bsc_segmentMdLF_neo(wbFG, fsDIR);
 
classification.index(LeftMdLFindexes)=27;
classification.index(RightMdLFindexes)=28;

[~, RightMeyerInd, ~,RightBaumInd, ~, LeftMeyerInd, ~,LeftBaumInd] =bsc_opticRadiationSeg_V2(wbFG, fsDIR);

classification.index(LeftMeyerInd)=29;
classification.index(RightMeyerInd)=30;
classification.index(LeftBaumInd)=31;
classification.index(RightBaumInd)=32;

% for itracts=1:length(classification.names)
%     spaceIndices=strfind(classification.names{itracts},' ');
%     classification.names{itracts}(spaceIndices)='_';
% end

disp('\n Tracts segmentation complete');

return
