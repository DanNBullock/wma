function [FullVL] =bsc_multiVirtualLesion(tractNameList,fiberIndexList,fe,FiberDir,saveHeader, nosave)
% function [FullVL] =bsc_multiVirtualLesion(tractNameList,fiberIndexList,fe)
%
% OVERVIEW:  This function generates a data structure wherein each field
% corresponds to a virtual lesion output for a fiber tract that has been
% entered into the function
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
% fe:  the FE structure that the virtual lesioning will be conducted upon.
% it is from the FE's tractome that the fiber tracts will be removed and
% the effect computed.
%
% -FiberDir: directory path for the directory you would like your fiber 
% indexes saved down to.
%
% -saveHeader: a string corresponding to whatever iformation you would like
% to serve as the file identifier for the FiberIndexes object.  Could
% include information like tractography parameter settings, subject
% identifier or group ID.  (example: 'HCP_105115_PROB_Lmax2_conn5')
%
% -nosave: flag for whether or not to save output.  Default = false, and
% therefore will save.
%
%
% OUTPUTS:
% -FullVL:  a structure wherein each field is named after a fiber tract.
% Each of these fields contains information like evidence, headerData (from
% saveHeader), tract name (from the tract list), and other data
% corresponding to the virtual lesion.  If there are no indexed fibers in
% the corresponding fiberIndexList this will be blank (i.e. if segmentation
% failed).
%
% % (C) Daniel Bullock 2017 Bloomington

%% Preliminaries
% loads file if a string was passed 
if ischar(fe)
    fe = load(fe);
else
end

%sets default saving behavior.  Defaults to saving.
if notDefined('nosave'), nosave=false;end


%% Virtual lesion loop
for iFGs = 1:length(fiberIndexList)
    if ~isempty(fiberIndexList{iFGs})
        [ vl.rmse_wVL, vl.rmse_woVL, vl.nFib_tract, vl.nFib_PN, vl.nVoxels ] = feComputeVirtualLesion(fe, fiberIndexList{iFGs});
        % feVirtualLesion for plots
        fevl = feComputeEvidence(vl.rmse_woVL, vl.rmse_wVL);
        vl.evidence=fevl;
        vl.headerData=saveHeader;
        vl.fiberName=tractNameList{iFGs};
        FullVL.tractNameList{iFGs}=vl;
        clear vl
    else
        FullVL=[];
    end
end

%the FullVL output structure is saved down
if ~nosave
save (strcat(fullfile(FiberDir),saveHeader,'_VLOutput.mat'),'FullVL','-v7.3');
end
fprintf ('\n multi lesion copmputation performed')
end