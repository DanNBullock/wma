function FiberIndexes= removeOutliersBSC(FiberIndexes,wbFG, centroidSD, lengthSD)
%
% FiberIndexes= removeOutliersBSC(FiberIndexes,wbFG)
%
% Performs outlier removal given input parameters defauts to 4 and 4 if
% nothing put in
%
% INPUTS: 
%
% -FiberIndexes: Either the path to the bsc FiberIndexes structure, or the
%  structure itself
%
% -wbFG: Either the path to the fe or WBFG, or the
%  structure itself
%
% -centroid SD: paramater for the sd deviance cut off from distance from
%  centroid of fibers.  Fibers outside of this distance will be pruned.
%
% -length SD: paramater for the sd deviance cut off for fiber length
%  fibers exceeding the threshold will be pruned.
%% preliminaries

if ischar(bscFiberIndexesFile)
    load(bscFiberIndexesFile);
end

% loads file if a string was passed 
if ischar(wbFG)
    wbFG = load(wbFG);
else
end

%if it is a fe structure, get the wbFG out of it
if isfield(wbFG, 'fe')
    wbFG = feGet(wbFG.fe, 'fibers acpc');
end

if notDefined('centroidSD')
    centroidSD=4;
end

if notDefined('lengthSD')
    lengthSD=4;
end

%% begin outlier removal
tractNames=fieldnames(FiberIndexes);
tractFG=wbFG;
tractFG.fibers=[];
for itracts=1:length(tractNames)
    tractFG.name=tractNames{itracts};
    tractFG.fibers={wbFG.fibers{FiberIndexes.(tractFG.name)}};
    tractFG.fibers=tractFG.fibers';
    [fg, keep]=mbaComputeFibersOutliers(tractFG,3,3,40);
    tractIndexes=FiberIndexes.(tractFG.name);
    FiberIndexes.(tractFG.name)=tractIndexes(keep);
end

end
    
    
