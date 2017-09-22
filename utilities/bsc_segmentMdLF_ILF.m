function [RightILF, RightILFIndexes, LeftILF, LeftILFIndexes, RightMdLFspl, RightMdLFsplIndexes, LeftMdLFspl, LeftMdLFsplIndexes,... 
    RightMdLFang, RightMdLFangIndexes, LeftMdLFang, LeftMdLFangIndexes] =bsc_segmentMdLF_ILF(wbfg, fsDir)
%
%[RightILF, RightILFIndexes, LeftILF, LeftILFIndexes, RightMdLFspl, RightMdLFsplIndexes, LeftMdLFspl, LeftMdLFsplIndexes,... 
%    RightMdLFang, RightMdLFangIndexes, LeftMdLFang, LeftMdLFangIndexes] =bsc_segmentMdLF_ILF(wbfg, fsDir)
%
% This function automatedly segments the middle longitudinal fasiculus
% from a given whole brain fiber group using the subject's 2009 DK
% freesurfer parcellation.

% Inputs:
% -wbfg: a whole brain fiber group structure
% -fsDir: path to THIS SUBJECT'S freesurfer directory

% Outputs:
% -RightMdLF: fiber structure for right middle longitudinal fasiculus
% -LeftMdLF: fiber structure for left middle longitudinal fasiculus

% -RightMdLFIndexes: fiber indexes into the given wbfg for the right middle longitudinal fasiculus
% -LeftMdLFIndexes: fiber indexes into the given wbfg for the left middle longitudinal fasiculus
%
%  Same for the other tracts
% (C) Daniel Bullock, 2017, Indiana University

%% parameter note & initialization

%maybe play with this if something isn't to your liking, it corresponds to
%the smoothing kernel used for the parietal and temporal ROIs
%smoothParameter=5;
% Unused
caudateLUT=[11,50];
for iFibers=1:length(wbfg.fibers)
    fiberNodeNum=round(length(wbfg.fibers{iFibers})/2);
    curStreamline=wbfg.fibers{iFibers};
    midpoints(iFibers,:)=curStreamline(:,fiberNodeNum);
    %endpoints1(iFibers,:)=curStreamline(:,1);
    %endpoints2(iFibers,:)=curStreamline(:,end);
    %streamLengths(iFibers)=sum(sqrt(sum(diff(wbfg.fibers{iFibers},1,2).^2)));
end

%iterates through left and right sides
for leftright= [1,2]
    
    %sidenum is basically a way of switching  between the left and right
    %hemispheres of the brain in accordance with freesurfer's ROI
    %numbering scheme. left = 1, right = 2
    sidenum=10000+leftright*1000;
    

    
    %% occipito roi
    %generates the roi for the occipito-parietal regions corresponding to
    %MdLF
    
    [mergedOCTROI] =bsc_roiFromFSnums(fsDir,[120,111,166, 143,119,158,145,122,102,119,159]+sidenum,1,5);
    mergedOCTROI.name='occipito';  
    
    %curtail somewhat
    [ipsROI] =bsc_roiFromFSnums(fsDir,157+sidenum,0);
    ipsPost=max(ipsROI.coords(:,2));
    posteriorPoints=find(mergedOCTROI.coords(:,2)<ipsPost);
    mergedOCTROI.coords=mergedOCTROI.coords(posteriorPoints ,: );
    
    %% parietal roi
    [mergedParieto] =bsc_roiFromFSnums(fsDir,[157,127,130,172]+sidenum,1,5);
    mergedParieto.name='parieto';
    

    
    %%  Cutoff
        insSup =bsc_roiFromFSnums(fsDir,[150]+sidenum,0);
        insSupcutoff=max(insSup.coords(:,3))-2;
        octpri=bsc_roiFromFSnums(fsDir,[166]+sidenum,0);
        octpriCut=max(octpri.coords(:,2));
        supmarg=bsc_roiFromFSnums(fsDir,[126]+sidenum,0);
        supmargCut=max(supmarg.coords(:,2));
        
    
    %% lateral temporal roi
    %generates the roi for the lateral-temporal regions corresponding to
    %MdLF
    
    [mergedLatTempROI] =bsc_roiFromFSnums(fsDir,[134, 144, 174,135]+sidenum,1,9);
    mergedLatTempROI.name='lateral-temporalROI'; 
    
    %impliments a cutoff to ensure that the temporal roi coordinates are
    %anterior of the y=-15 plane
    %[amygdlaROI] =bsc_roiFromFSnums(fsDir,[134, 122, 144, 174, 114, 135]+sidenum,1,9);
    amygdlaIDs=[18,54];
    [amygdlaROI] =bsc_roiFromFSnums(fsDir,amygdlaIDs(leftright));
    amygdalaPost=min(amygdlaROI.coords(:,2));
    
    
    anteriorPoints=find(mergedLatTempROI.coords(:,2)>amygdalaPost);
    mergedLatTempROI.coords=mergedLatTempROI.coords(anteriorPoints ,: );
    %try and use posterior of amygdla as cut off
    
    %% IPL
   [mergedIPLROI] =bsc_roiFromFSnums(fsDir,[126, 168,125]+sidenum,1,3);
    mergedIPLROI.name='IPLROI'; 
    
    %%
    %mdlfAngCorrection

    %mdlfAngCorrect=bsc_roiFromFSnums(fsDir,[118,124,164,112,169,104,151]+sidenum,1,3);
  

     %mdlfAngCorrection

    %mdlfSplCorrect=bsc_roiFromFSnums(fsDir,[118,124,164,112,169,104,151]+sidenum,1,3);
    
    %ILF correct
    %ILFCorrect=bsc_roiFromFSnums(fsDir,[[118,124,164,112,169,104,151]+sidenum ,49,10,53,17,60,28],1,3); %,125,165,148,153,146,128,168,129,157,147,108,116,130,107,123,130,172
    
    SubCortCorrect=bsc_roiFromFSnums(fsDir,[49,10,53,17,60,28],1,3); 

    %% segmenting
    
    
    %switch for correct name
    if leftright==2
        sideflag='R';
    else
        sideflag='L';
    end
    %currentFascicleName=strcat(sideflag,'_MdLF');
    
    %create object containing all rois
    currentROIs1= [{mergedLatTempROI} {mergedOCTROI} ];
    currentROIs2= [{mergedLatTempROI} {mergedParieto} ];
    currentROIs3= [{mergedLatTempROI} {mergedIPLROI}];
    
    
    
    %actually segment
    [ILF, ILFFiberBoolVec]=bsc_tractByEndpointROIs(wbfg, currentROIs1);
    [ILF, removeILFIND]=wma_SegmentFascicleFromConnectome(ILF, [{SubCortCorrect}], {'not'}, 'dud');
    ILFPreIndex=find(ILFFiberBoolVec);
    removeThese=ILFPreIndex(~removeILFIND);
    ILFFiberBoolVec(removeThese)=false;
    ILFFiberBoolVec=ILFFiberBoolVec' & midpoints(:,3)<insSupcutoff &  midpoints(:,2)>octpriCut & midpoints(:,2)<supmargCut;
    ILF.fibers=wbfg.fibers(ILFFiberBoolVec);
    
    [MDLFANG, MdLFangFiberBoolVec]=bsc_tractByEndpointROIs(wbfg, currentROIs3);
    [MDLFANG, removeMDLFANGIND]=wma_SegmentFascicleFromConnectome(MDLFANG, [{SubCortCorrect}], {'not'}, 'dud');
    MDLFANGPreIndex=find(MdLFangFiberBoolVec);
    removeThese=MDLFANGPreIndex(~removeMDLFANGIND);
    MdLFangFiberBoolVec(removeThese)=false;
    MdLFangFiberBoolVec=MdLFangFiberBoolVec' & midpoints(:,3)<insSupcutoff &  midpoints(:,2)>octpriCut & midpoints(:,2)<supmargCut;
    MDLFANG.fibers=wbfg.fibers(MdLFangFiberBoolVec);
    
    
    [MDLFSPL, MdLFsplFiberBoolVec]=bsc_tractByEndpointROIs(wbfg, currentROIs2);
    [MDLFSPL, removeMDLFSPLIND]=wma_SegmentFascicleFromConnectome(MDLFSPL, [{SubCortCorrect}], {'not'}, 'dud');
    MDLFSPlPreIndex=find(MdLFsplFiberBoolVec);
    removeThese=MDLFSPlPreIndex(~removeMDLFSPLIND);
    MdLFsplFiberBoolVec(removeThese)=false;
    MdLFsplFiberBoolVec=MdLFsplFiberBoolVec' & midpoints(:,3)<insSupcutoff &  midpoints(:,2)>octpriCut & midpoints(:,2)<supmargCut;
    MDLFSPL.fibers=wbfg.fibers(MdLFsplFiberBoolVec);

    %directs segmentation output to correct function output holder
    if leftright == 2
        RightILF=ILF;
        RightILF.name='Right ILF';
        RightILFIndexes=ILFFiberBoolVec';
        
        RightMdLFspl=MDLFSPL;
        RightMdLFspl.name='Right MdLF-SPL';
        RightMdLFsplIndexes=MdLFsplFiberBoolVec';
        
        RightMdLFang=MDLFANG;
        RightMdLFang.name='Right MdLF-ANG';
        RightMdLFangIndexes=MdLFangFiberBoolVec';
        
        
    else
        LeftILF=ILF;
        LeftILF.name='Left ILF';
        LeftILFIndexes=ILFFiberBoolVec;
        
        LeftMdLFspl=MDLFSPL;
        LeftMdLFspl.name='Left MdLF-SPL';
        LeftMdLFsplIndexes=MdLFsplFiberBoolVec;
        
        LeftMdLFang=MDLFANG;
        LeftMdLFang.name='Left MdLF-ANG';
        LeftMdLFangIndexes=MdLFangFiberBoolVec;
    end
    
end



end