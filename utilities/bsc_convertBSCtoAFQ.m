function classification = bsc_convertBSCtoAFQ (bscFiberIndexesFile, wbFG)
%
% classification = bsc_convertBSCtoAFQ (bscFiberIndexesFile, wbFG)
%
% Converts bsc tract indexing format to AFQ type structure.
%
% bscFiberIndexesFile= either the FiberIndexes structure or the path to it
%
% wbFG= Either the FE or the WBFG, path or object

%%Preliminaries
if ischar(bscFiberIndexesFile)
    load(bscFiberIndexesFile);
    bscFiberIndexesFile=FiberIndexes;
end

% loads file if a string was passed 
if ischar(wbFG)
    wbFG = load(wbFG);
    testName=fieldnames(wbFG);
    if length(fie
else
end


%if it is a fe structure, get the wbFG out of it
if isfield(wbFG, 'fe')
    wbFG = feGet(wbFG.fe, 'fibers acpc');
end

%% Start Conversion

classification.index=zeros(length(wbFG.fibers),1);

fiberNames=fieldnames(FiberIndexes);
classification.names=fiberNames;

for itracts=1:length(classification.names)
    indexes=FiberIndexes.(classification.names{itracts});
    classification.index(indexes)=itracts;
end

end