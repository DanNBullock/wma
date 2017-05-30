function bsc_multiTractVolumeFE(tractNameList,fiberIndexList,wbFG,FiberDir,saveHeader, fsDIR)
% function bsc_multiTractVolumeFE(tractNameList,fiberIndexList,fe,FiberDir,saveHeader, fsDIR)
%
% OVERVIEW:  This function creates endpoint volumes for all fibers that are
% passed to it.  Implicitly it assumes that the fiber indexes are for
% positive only fibers (seeing as how the endpoint volume is a quantative
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
% -wbFG: the wbFG structure from which the streamlines will be extracted. 
% If you pass a fe structure or an fe path it will (probably) extract the
% wbFG from it appropriately. 
%
% -FiberDir: directory path for the directory you would like your fiber 
% indexes saved down to.
%
% -saveHeader: a string corresponding to whatever iformation you would like
% to serve as the file identifier for the fiber tract object.  Could
% include information like tractography parameter settings, subject
% identifier or group ID.  (example: 'HCP_105115_PROB_Lmax2_conn5')
% -nosave: flag for whether or not to save output.  Default = false, and
% therefore will save.
%
% OUTPUTS:
% none, it saves the volumes down.
%
% % (C) Daniel Bullock 2017 Bloomington
%
%% initalization stuff
% loads file if a string was passed 
if ischar(wbFG)
    wbFG = load(wbFG);
else
end

%if it is a fe structure, get the wbFG out of it
if isfield(wbFG, 'fg')
    wbFG = feGet(wbFG, 'fibers acpc');
end

%% volume creation loop
for iFGs = 1:length(fiberIndexList)
    fg=wbFG;
    fg.fibers=wbFG.fibers(fiberIndexList{iFGs});
    fg.name=strcat(saveHeader,'_',tractNameList{iFGs});
    %not really catching any of this, as it is being saved down
    [nii, nii2, nii_normalized, nii2_normalized]=csc_endpointMaps_Decay_v3(fg,fsDIR,FiberDir, 3, 0, 'uniform');
end
    
end


