function [tracts_indexes] = wma_wrapper(wbFG,dt6,FiberDir,saveHeader,fsDIR, nosave)
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

%if it is a fe structure, get the wbFG out of it
if isfield(wbFG, 'fg')
    wbFG = feGet(wbFG, 'fibers acpc');
end

if ischar(dt6)
    dt6 = dtiLoadDt6(dt6);
else
    dti6=dt6;
end

% Sets default saving behavior.  Defaults to saving.
if notDefined('nosave'), nosave=false;end

%% Segmentation

disp('Segmenting Major and Associative Tracts');

% Segment the major white matter tracts in the Mori Atlas
tic
[fg_classified,~,classification,~] = wma_majortracts(dti_file, wbFG);
toc=segmentTime;
fprintf ('\n Mori Atlas segmentation run complete in %4.2f hours \n', segmentTime/(60*60))

% update name field
newTractNames={'L_VOF','R_VOF','L_pArc','R_pArc','L_TPC','R_TPC','L_MdLF','R_MdLF'};
for iNEWtracts=21:28
   classification.names{iNEWtracts}=newTractNames{iNEWtracts-20};
end

% posterior arcuate and temporo-parietal connection segmentation
[L_pArc, ~, L_pArc_Indexes,L_TPC_Indexes ,R_pArc, ~, R_pArc_Indexes,R_TPC_Indexes] = bsc_automated_roi_segment_script_neo(wbFG,fsDIR);

classification(L_pArc_Indexes)=23;
classification(R_pArc_Indexes)=24;
classification(L_TPC_Indexes)=25;
classification(R_TPC_Indexes)=26;

% Segment the Vertical Occipital Fasiculus (VOF)
[~, ~, L_VOF_Indexes, R_VOF_Indexes] =  bsc_segmentVOF(wbFG, fsDIR, FiberDir, dtiFile, fg_classified(19), fg_classified(20), L_pArc, R_pArc, L_pArc_Indexes, R_pArc_Indexes);

classification(L_VOF_Indexes)=21;
classification(R_VOF_Indexes)=22;

% Middle Longitudinal Fasiculus segmentation
[~, RightMdLFindexes, ...
 ~, LeftMdLFindexes] = bsc_segmentMdLF_neo(wbFG, fsDIR);
 
classification(LeftMdLFindexes)=27;
classification(RightMdLFindexes)=28;

disp('Tracts segmentation complete');

if ~nosave
    
    save (strcat(fullfile(FiberDir),saveHeader,'classification.mat'),'classification','-v7.3');
    fprintf(' \n classification saved \n')
end

return
