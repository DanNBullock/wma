function [RightAslant, RightAslantIndexes, LeftAslant, LeftAslantIndexes] =bsc_segmentAslantProto(wbfg, fsDir)
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


%iterates through left and right sides
for leftright= [1,2]
    
    %sidenum is basically a way of switching  between the left and right
    %hemispheres of the brain in accordance with freesurfer's ROI
    %numbering scheme. left = 1, right = 2
    sidenum=10000+leftright*1000;
    

    
    %% occipito roi
    %generates the roi for the occipito-parietal regions corresponding to
    %MdLF
    
    [superiorROI] =bsc_roiFromFSnums(fsDir,[116]+sidenum,1,7);
    superiorROI.name='superior';  
    

    
    %% parietal roi
    [lateralROI] =bsc_roiFromFSnums(fsDir,[18]+sidenum-10000,1,5);
    lateralROI.name='lateral';
    

    
    %%  Cutoff
    [triangular] =bsc_roiFromFSnums(fsDir,[114]+sidenum,0,[]);
    posteriorTriangular=min(triangular.coords(:,2));
    atlasNifti = wma_getAsegFile(fsDir , '2009');
        
    [anteriorBorder] =bsc_makePlanarROI(atlasNifti,posteriorTriangular,'y');
    
    %posterior border
    [paracentral] =bsc_roiFromFSnums(fsDir,[103]+sidenum,0,[]);
    posteriorBorderPoint=max(paracentral.coords(:,2));
    atlasNifti = wma_getAsegFile(fsDir , '2009');
        
    [posteriorBorder] =bsc_makePlanarROI(atlasNifti,posteriorBorderPoint,'y');
    
    [notCross] =bsc_makePlanarROI(atlasNifti,0,'x');
    
     [triangularNot] =bsc_roiFromFSnums(fsDir,[114]+sidenum,1,9);
    
    
     [notROI] =bsc_roiFromFSnums(fsDir,[114 150 155 118 115]+sidenum,1,3);
    %% segmenting
    
    
    %switch for correct name
    if leftright==2
        sideflag='R';
    else
        sideflag='L';
    end
    %currentFascicleName=strcat(sideflag,'_MdLF');
    
    %create object containing all rois
    currentROIs1= [{superiorROI} {lateralROI} ];

    
    
    %actually segment
    [aslant, aslantFiberBoolVec]=bsc_tractByEndpointROIs(wbfg, currentROIs1);
    
    [fascicle, FiberBoolVec] = wma_SegmentFascicleFromConnectome(wbfg, [{anteriorBorder} {notROI} {posteriorBorder} {triangularNot} {notCross}], {'not', 'not', 'not', 'not', 'not'}, 'blank');
    
    
    for ifibers=1:length(aslant.fibers)
        midpointbool(ifibers)=aslant.fibers{ifibers}(2,round(length(aslant.fibers{ifibers})/2))<posteriorTriangular;
    end
    eliminateIndexes=find(~midpointbool);
    aslantIndexes=find(aslantFiberBoolVec);
    aslantFiberBoolVec(aslantIndexes(eliminateIndexes))=0;
    
    
    aslant.fibers=wbfg.fibers(FiberBoolVec &aslantFiberBoolVec');
    


    %directs segmentation output to correct function output holder
    if leftright == 2
        RightAslant=aslant;
        RightAslant.name='Right Aslant';
        RightAslantIndexes=FiberBoolVec &aslantFiberBoolVec';

        
        
    else
        LeftAslant=aslant;
        LeftAslant.name='Left Aslant';
        LeftAslantIndexes=FiberBoolVec &aslantFiberBoolVec';

    end
    
    clear midpointbool
    clear eliminateIndexes
    clear aslantFiberBoolvec
    
end



end