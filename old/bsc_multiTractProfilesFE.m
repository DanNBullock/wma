function [DiffusionProfile]=bsc_multiTractProfilesFE(tractNameList,fiberIndexList,wbFG,FiberDir,saveHeader,dt6, nosave)
% function [DiffusionProfile]=bsc_multiTractProfilesFE(tractNameList,fiberIndexList,wbFG,FiberDir,saveHeader,dt6, nosave)
%
% OVERVIEW:  This function runs dtiComputeDiffusionPropertiesAlongFG for all fibers that are
% passed to it.  Implicitly it assumes that the fiber indexes are for
% positive only fibers (seeing as how the diffusivity metrics are a quantative
% measure of sorts)
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
% -wbFG: either a path to a saved wbFG or a wbFG object.  If you pass a fe
% structure or an fe path it will (probably) extract the fg from it.
% appropriately.
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
% -nosave: flag for whether or not to save output.  Default = false, and
% therefore will save.
%
% -OUTPUTS
% -DiffusionProfile: a structure wherein each field is named after a fiber tract.
% Each of these fields contains information like fractional anisotropy,
% mean diffusivity, radial diffusivity, and axial diffusivity.  If there 
% are no indexed fibers in % the corresponding fiberIndexList this will be
% blank (i.e. if segmentation failed).
%
% % (C) Daniel Bullock 2017 Bloomington

%% preliminaries
% set nosave flag
if notDefined('nosave'), nosave=false;end

% loads whatver was passed to it with wbFG
if ischar(wbFG)
    wbFG = load(wbFG);
else
end

% loads the dt6 file if a path was passed
if ischar(dt6)
    dtiFile = dtiLoadDt6(dt6);
else
    dtiFile=dt6;
end

%if it is a fe structure, get the wbFG out of it
if isfield(wbFG, 'fg')
    wbFG = feGet(wbFG, 'fibers acpc');
end
%% ComputeDiffusion loop
for iFGs = 1:length(fiberIndexList)
    if ~isempty(fiberIndexList{iFGs})
    fg=wbFG;
    fg.fibers=wbFG.fibers(fiberIndexList{iFGs});
    fg.name=strcat(saveHeader,'_',tractNameList{iFGs});
    
    [fa, md, rd, ad, cl, SuperFiber, fgClipped, cp, cs, fgResampled]=dtiComputeDiffusionPropertiesAlongFG(fg, dt6, [] ,[] , 100);
    DiffusionProfile.(tractNameList{iFGs}).FractionalAnisotropy=fa;
    DiffusionProfile.(tractNameList{iFGs}).MeanDiffusivity=md;
    DiffusionProfile.(tractNameList{iFGs}).RadialDiffusivity=rd;
    DiffusionProfile.(tractNameList{iFGs}).AxialDiffusivity=ad;
    else
        DiffusionProfile.(tractNameList{iFGs})=[];
    end
end

%saving, if desired
if ~nosave
    save (strcat(fullfile(FiberDir),saveHeader,'_DiffusionProfile.mat'),'DiffusionProfile','-v7.3');
end
fprintf ('\n DiffusionProfile copmputations performed \n')