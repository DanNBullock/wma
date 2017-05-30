function [FiberIndexes]=bsc_AFQseg_and_bloomtracks(wbFG,dt6path,FiberDir,saveHeader,fsDIR, nosave)
%[FiberIndexes]=bsc_AFQseg_and_bloomtracks(wbFG,dt6path,FiberDir,saveHeader)
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
% -dt6path: either a path to a saved dt6 file or a dt6 object
%
% -FiberDir: directory path for the directory you would like your fiber 
% indexes saved down to.
%
% -saveHeader: a string corresponding to whatever iformation you would like
% to serve as the file identifier for the FiberIndexes object.  Could
% include information like tractography parameter settings, subject
% identifier or group ID.  (example: 'HCP_105115_PROB_Lmax2_conn5')
%
% -fsDIR: path  to THIS SUBJECT'S freesurfer directory
%
% -nosave=flag for whether or not to save the output.  Defaults action is
% to save
%
% OUTPUTS:
% -FiberIndexes:  A structure whose field names correspond to the names
%  of the various fiber tracts that were segmented in this function.  Each
%  field stores the indexes into the wbFG corresponding to the named track.
%  i.e. FiberIndexes.R_pArc will correspond to several hundered (probably)
%  indexes into the wbFG for the right posterior arcuate.
%
%  NOTE: This function does not save down MERELY validated tracts (i.e. via
%  LiFE), nor does it run an outlier removal function (i.e
%  mbaComputeFibersOutliers).  The output indexes are "unaltered" in this
%  regard.  The only exception to this is the VOF which has had outlier
%  removal at 2 for gradient (i.e. orientation)
%
% (C) Daniel Bullock, 2017, Indiana University
%% initialization stuff

% loads file if a string was passed 
if ischar(wbFG)
    wbFG = load(wbFG);
else
end

% loads the dt6 file if a path was passed
if ischar(dt6path)
    dtiFile = dtiLoadDt6(dt6path);
else
    dtiFile=dt6path;
end

%if it is a fe structure, get the wbFG out of it
if isfield(wbFG, 'fg')
    wbFG = feGet(wbFG, 'fibers acpc');
end

%sets default saving behavior.  Defaults to saving.
if notDefined('nosave'), nosave=false;end

%% afq segmentation
tic
[fg_classified,~,classification,~]=AFQ_SegmentFiberGroups(dtiFile, wbFG);
afqTime=toc;
fprintf ('\n AFQ run complete in %i hours \n', afqTime/(60*60))
%recomend a keyboard here to inspect output of AFQ
%note that the afq outputs will be used in subsequent steps.
keyboard
%% vertical occipital fasiculus segmentation (with gradient adjustment)
%fsROIdir=fullfile(fsDIR,'/label/');
[L_VOF, R_VOF, L_VOF_Indexes, R_VOF_Indexes,  L_pArc, R_pArc, L_pArc_fibersIndices, R_pArc_fibersIndices, L_pArc_vot, R_pArc_vot, L_pArc_vot_Indexes, R_pArc_vot_Indexes]=bsc_segmentVOF(wbFG, fsDIR, FiberDir, dtiFile, fg_classified(19), fg_classified(20));
%NOTE: we are not using the posterior arcuate and vot tracts that come out
%of this function (they aren't that great?).  ALSO WE ARE NOW RELYING ON A
%NEW VERSION OF AFQ_FindVerticalFibers WHICH DOES NOT REQUIRE RUNNING
%fs_roisFromAllLabels FIRST.
fprintf('\n vof extraction complete \n')

%% posterior arcuate and temporo-parietal connection segmentation
[L_pArc, L_TPC, L_pArc_Indexes,L_TPC_Indexes ,R_pArc, R_TPC, R_pArc_Indexes, R_TPC_Indexes]=bsc_automated_roi_segment_script_neo(wbFG,fsDIR);
fprintf('\n pArc and TPC segmentation complete \n');



%% Middle Longitudinal Fasiculus segmentation
[RightMdLF, RightMdLFindexes, LeftMdLF, LeftMdLFindexes] =bsc_segmentMdLF_neo(wbFG, fsDIR);
fprintf('MdLF segmentation complete \n');

%% create fiber index structure
% afq tracts
FiberIndexes.L_ThalamicRadiation=find(classification.index==1);
FiberIndexes.R_ThalamicRadiation=find(classification.index==2);
FiberIndexes.L_CorticoSpinalTract=find(classification.index==3);
FiberIndexes.R_CorticoSpinalTract=find(classification.index==4);
FiberIndexes.L_CingulateCingulum=find(classification.index==5);
FiberIndexes.R_CingulateCingulum=find(classification.index==6);
FiberIndexes.L_CingulumHippocampus=find(classification.index==7);
FiberIndexes.R_CingulumHippocampus=find(classification.index==8);
FiberIndexes.CallosumForceMajor=find(classification.index==9);
FiberIndexes.CallosumForceMinor=find(classification.index==10);
FiberIndexes.L_IFOF=find(classification.index==11);
FiberIndexes.R_IFOF=find(classification.index==12);
FiberIndexes.L_ILF=find(classification.index==13);
FiberIndexes.R_ILF=find(classification.index==14);
FiberIndexes.L_SLF=find(classification.index==15);
FiberIndexes.R_SLF=find(classification.index==16);
FiberIndexes.L_Uncinate=find(classification.index==17);
FiberIndexes.R_Uncinate=find(classification.index==18);
FiberIndexes.L_Arcuate=find(classification.index==19);
FiberIndexes.R_Arcuate=find(classification.index==20);

%bloomington tracts + VOF
FiberIndexes.L_VOF=L_VOF_Indexes;
FiberIndexes.R_VOF=R_VOF_Indexes;
FiberIndexes.L_pArc=L_pArc_Indexes;
FiberIndexes.R_pArc=R_pArc_Indexes;
FiberIndexes.L_TPC=L_TPC_Indexes;
FiberIndexes.R_TPC=R_TPC_Indexes;
FiberIndexes.R_MdLF=RightMdLFindexes;
FiberIndexes.L_MdLF=LeftMdLFindexes;

if ~nosave
save (strcat(fullfile(FiberDir),saveHeader,'_FiberIndexes.mat'),'FiberIndexes','-v7.3');
end

fprintf(' \n Fiber Indexes saved \n');

end
