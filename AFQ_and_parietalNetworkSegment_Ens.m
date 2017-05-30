function AFQ_and_parietalNetworkSegment_Ens(groupOptionNum,subjectOptionNum)
%  Quick and dirty loop for multi segmentation

%% path generation
%GROUP OPTIONS IS SPECIFIC TO PESTILLILAB DATA SETUP
groupOptions={'HCP','Stanford'};
lmaxOptions=(2:2:12);
%CON OPTIONS IS SPECIFIC TO PESTILLILAB DATA SETUP
connOptions=1:10;
trackingOptions={'STREAM','PROB'}; %save 'tensor' for another day

group=groupOptions{groupOptionNum};
%outDir simply names the directory which the
outDir='VPF_Master_Script';
% lmax=lmaxOptions(lmaxOptionNum);
tracking=trackingOptions{2};
% conn=connOptions(1);
isubjects=subjectOptionNum;
%start loop here

for iconn=1:10 %:length(connOptions)
    %we use try here because if we run this in a loop we don't want the
    %whole attempt to crash if one of the reconstructions fails.
    
    %generate the file paths to all of
    [subjects, dt6path, t1path, wbFGpath, freeSurfROIsPath, fullFiberOutDir, fullFigureOutDir,  dwidata, fePath] = bsc_HCP_Stanford_PathGen_Permute4(group,outDir, 10, tracking,iconn);
    
    %generate the filename tag which will append the front of all
    %output files.  Identifies all of the various input parameters
    %assocaited with the input data.
    saveHeader=strcat(group,'_',subjects{isubjects},'_Ens','_conn',num2str(iconn));
    %basically fprintf to tell us what iteration we are on
    saveHeader
    %more quick and dirty to get fsDIR
    fsROIDirSlashIndexes=strfind(freeSurfROIsPath{isubjects}, '/');
    fsDIR=freeSurfROIsPath{isubjects}(fsROIDirSlashIndexes(1):fsROIDirSlashIndexes(end-1));
    
    
    
    %% segmentation stuff
    %load FE and get fibers
    load(strcat(fullFiberOutDir{isubjects},'subj_',subjects{isubjects},'_probabalistic_ensemble_conn',num2str(iconn),'.mat'));
    %not necessary to get the wbFG out here.
%     wbFG = feGet(fe, 'fg acpc');
    fprintf('\n wbFG extracted')
    
    [FiberIndexes]=bsc_AFQseg_and_bloomtracks(fe,dt6path{isubjects},fullFiberOutDir{isubjects},saveHeader,fsDIR);
    %% build lists for later functions
    tractNameList=fieldnames(FiberIndexes);
    
    % building the fiberIndexList object
    for iFieldName=1:length(tractNameList)
        fiberIndexList{iFieldName}=FiberIndexes.(tractNameList{iFieldName});
    end
    
    [summaryStructure] = bsc_multiTractAnalysis(tractNameList,fiberIndexList,fe,fullFiberOutDir{isubjects},saveHeader,dt6path{isubjects},fsDIR);
    
end
end
