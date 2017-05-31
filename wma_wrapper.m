function [tracts_indexes] = wma_wrapper(wbFG,dti_file,FiberDir,saveHeader,fsDIR, nosave)
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
% -dt6path: either a path to a saved dt6 file or a dt6 object
%
% -FiberDir: directory path for the directory you would like your fiber 
% indexes saved down to.
%
% -saveHeader: a string corresponding to whatever iformation you would like
% to serve as the file identifier for the tracts_indexes object.  Could
% include information like tractography parameter settings, subject
% identifier or group ID.  (example: 'HCP_105115_PROB_Lmax2_conn5')
%
% -fsDIR: path  to THIS SUBJECT'S freesurfer directory
%
% -nosave=flag for whether or not to save the output.  Defaults action is
% to save
%
% OUTPUTS:
% -tracts_indexes:  A structure whose field names correspond to the names
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

% loads file if a string was passed 
if ischar(wbFG)
    wbFG = load(wbFG);
else
end

%if it is a fe structure, get the wbFG out of it
if isfield(wbFG, 'fg')
    wbFG = feGet(wbFG, 'fibers acpc');
end

% Sets default saving behavior.  Defaults to saving.
if notDefined('nosave'), nosave=false;end

disp('Segmenting Major and Associative Tracts');

% posterior arcuate and temporo-parietal connection segmentation
[L_pArc, ~, L_pArc_Indexes,L_TPC_Indexes ,R_pArc, ~, R_pArc_Indexes,R_TPC_Indexes] = bsc_automated_roi_segment_script_neo(wbFG,fsDIR);

% Segment the major white matter tracts in the Mori Atlas
[fg_classified,~,classification,~] = wma_majortracts(dti_file, wbFG);

% Segment the Vertical Occipital Fasiculus (VOF)
[~, ~, L_VOF_Indexes, R_VOF_Indexes] =  bsc_segmentVOF(wbFG, fsDIR, FiberDir, dtiFile, fg_classified(19), fg_classified(20), L_pArc, R_pArc, L_pArc_Indexes, R_pArc_Indexes);

% Middle Longitudinal Fasiculus segmentation
[~, RightMdLFindexes, ...
 ~, LeftMdLFindexes] = bsc_segmentMdLF_neo(wbFG, fsDIR);

disp('Tracts segmentation complete');

% Colelct all tracts in  signle strucutre
% Mori tracts from the JHU Atrlas (Zhang 2008)
tracts_indexes.L_ThalamicRadiation=find(  classification.index==1);
tracts_indexes.R_ThalamicRadiation=find(  classification.index==2);
tracts_indexes.L_CorticoSpinalTract=find( classification.index==3);
tracts_indexes.R_CorticoSpinalTract=find( classification.index==4);
tracts_indexes.L_CingulateCingulum=find(  classification.index==5);
tracts_indexes.R_CingulateCingulum=find(  classification.index==6);
tracts_indexes.L_CingulumHippocampus=find(classification.index==7);
tracts_indexes.R_CingulumHippocampus=find(classification.index==8);
tracts_indexes.CallosumForceMajor=find(   classification.index==9);
tracts_indexes.CallosumForceMinor=find(   classification.index==10);
tracts_indexes.L_IFOF=find(classification.index==11);
tracts_indexes.R_IFOF=find(classification.index==12);
tracts_indexes.L_ILF=find( classification.index==13);
tracts_indexes.R_ILF=find( classification.index==14);
tracts_indexes.L_SLF=find( classification.index==15);
tracts_indexes.R_SLF=find( classification.index==16);
tracts_indexes.L_Uncinate=find(classification.index==17);
tracts_indexes.R_Uncinate=find(classification.index==18);
tracts_indexes.L_Arcuate=find( classification.index==19);
tracts_indexes.R_Arcuate=find( classification.index==20);

% Vertical associative tracts
tracts_indexes.L_VOF=L_VOF_Indexes;
tracts_indexes.R_VOF=R_VOF_Indexes;
tracts_indexes.L_pArc=L_pArc_Indexes;
tracts_indexes.R_pArc=R_pArc_Indexes;
tracts_indexes.L_TPC=L_TPC_Indexes;
tracts_indexes.R_TPC=R_TPC_Indexes;
tracts_indexes.R_MdLF=RightMdLFindexes;
tracts_indexes.L_MdLF=LeftMdLFindexes;

if ~nosave
   save (strcat(fullfile(FiberDir),saveHeader,'_FiberIndexes.mat'), ...
        'tracts_indexes','-v7.3');
end

fprintf(' \n Fiber Indexes saved \n');

return
