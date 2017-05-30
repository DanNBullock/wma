function [LengthStruc] =bsc_multiTractLengthFG(tractNameList,fiberIndexList,wbFG,FiberDir,saveHeader, nosave)
% function [LengthStruc] =bsc_multiTractLengthFG(tractNameList,fiberIndexList,wbFG,FiberDir,saveHeader, nosave)
%
% OVERVIEW:  This function generates a data structure wherein each field
% corresponds to a fiber name, under which a number of quantative details
% are stored, including average length, volume and streamline count.
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
% wbFG:  a whole brain connectome, which the fiber index list corresponds to.
% If you pass a fe structure or an fe path it will (probably) extract 
% the wbFG from it.
%
% -FiberDir: directory path for the directory you would like your fiber 
% indexes saved down to.
%
% -saveHeader: a string corresponding to whatever iformation you would like
% to serve as the file identifier for the FiberIndexes object.  Could
% include information like tractography parameter settings, subject
% identifier or group ID.  (example: 'HCP_105115_PROB_Lmax2_conn5')
%
% OUTPUTS:
% -LengthStruc:  a structure wherein each field is named after a fiber tract.
% Each of these fields contains information like mean, sd , count and 
% a list of fiber lengths corresponding to the virtual lesion.  If there 
% are no indexed fibers in % the corresponding fiberIndexList this will be
% blank (i.e. if segmentation failed).
%
% % (C) Daniel Bullock 2017 Bloomington, Indiana
%% initialization stuff

% loads file if a string was passed 
if ischar(wbFG)
    wbFG = load(wbFG);
else
end

%if it is a fe structure, get the wbFG out of it
if isfield(wbFG, 'fg')
    wbFG = feGet(wbFG, 'fibers acpc');
end

%sets default saving behavior.  Defaults to saving.
if notDefined('nosave'), nosave=false;end

%% Virtual lesion loop
for iFGs = 1:length(fiberIndexList)
    if ~isempty(fiberIndexList{iFGs})
        volVec=[];
        tractFibers=wbFG.fibers(fiberIndexList{iFGs});
        for istreamlines=1:length(tractFibers)
            streamLengths(istreamlines)=sum(sqrt(sum(diff(tractFibers{istreamlines},1,2).^2)));
            volVec=horzcat(volVec,tractFibers{istreamlines});
        end
        LengthStruc.(tractNameList{iFGs}).streamLengths=streamLengths;
        LengthStruc.(tractNameList{iFGs}).mean=mean(streamLengths);
        LengthStruc.(tractNameList{iFGs}).std=std(streamLengths);
        LengthStruc.(tractNameList{iFGs}).count=length(streamLengths);
        % volume, as indicated by number of 1mm voxels occupied by at least 1 tract
        % node.  In theory if you wanted to resample this computation (i.e. how
        % many .5 mm voxels are occpupied by at least 1 tract node) you could
        % just multiply the volVec by a factor.  I.e. .5mm -> 2 or .25 mm -> 4
        LengthStruc.(tractNameList{iFGs}).volume=length(unique(floor(volVec'),'rows'));
        
        clear streamLengths
        
    else
        LengthStruc.(tractNameList{iFGs})=[];
    end

%the LengthStruc output structure is saved down
if ~nosave
save (strcat(fullfile(FiberDir),saveHeader,'_LengthStruc.mat'),'LengthStruc','-v7.3');
end
fprintf ('\n length copmputation performed \n')
end
