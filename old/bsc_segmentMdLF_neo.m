function [RightMdLF, RightMdLFindexes, LeftMdLF, LeftMdLFindexes] =bsc_segmentMdLF_neo(wbfg, fsDir)
%[RightMdLF, RightMdLFindexes, LeftMdLF, LeftMdLFindexes] =bsc_segmentMdLF_neo(wbfg, fsDir)
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

% (C) Daniel Bullock, 2017, Indiana University

%% parameter note & initialization

%maybe play with this if something isn't to your liking, it corresponds to
%the smoothing kernel used for the parietal and temporal ROIs
smoothParameter=5;

%this version relies on the aparc.a2009s+aseg file.  Here we make a
%nii.gz version of one doesn't already exist.  keyboards out if there's
%an error doing this
if ~exist(strcat(fsDir,'/mri/aparc.a2009s+aseg.nii.gz'))
    %apaprently necessary for matlab?
    spaceChar={' '};
    [status result] = system(strcat('mri_convert',spaceChar,fsDir,'/mri/aparc.a2009s+aseg.mgz',spaceChar, fsDir, '/mri/aparc.a2009s+aseg.nii.gz'));
    if status~=0
        error('/n Error generating aseg nifti file.  There may be a problem finding the aparc.a2009s+aseg file.')
        
    end
end

%reads in label data
labelNifti=niftiRead(strcat(fsDir, '/mri/aparc.a2009s+aseg.nii.gz'));
%generates a blank boolean array that corresponds to the dimensions of
%the fs aparc + aseg nifti data.
sizeLabelNifti=size(labelNifti.data);
blankLabelNifti(1:sizeLabelNifti(1),1:sizeLabelNifti(2),1:sizeLabelNifti(3))=false;

%iterates through left and right sides
for leftright= [1,2]
    
    %sidenum is basically a way of switching  between the left and right
    %hemispheres of the brain in accordance with freesurfer's ROI
    %numbering scheme. left = 1, right = 2
    sidenum=10000+leftright*1000;
    
    %% occipito-parietal roi
    %generates the roi for the occipito-parietal regions corresponding to
    %MdLF
    
    [mergedOCTPROI] =bsc_roiFromFSnums(fsDir,[120,111,166, 127]+sidenum,1,smoothParameter);
    mergedOCTPROI.name='occipito-parietalROI';  
    
    %impliments a cutoff to prevent the occiptito-parietal roi from
    %going too high up on the cortex
    inferiorPoints=find(mergedOCTPROI.coords(:,2)<45);
    mergedOCTPROI.coords=mergedOCTPROI.coords(inferiorPoints ,: );
    
    %% lateral temporal roi
    %generates the roi for the lateral-temporal regions corresponding to
    %MdLF
    
    [mergedLatTempROI] =bsc_roiFromFSnums(fsDir,[134, 122, 144, 174]+sidenum,1,smoothParameter);
    mergedLatTempROI.name='lateral-temporalROI'; 
    
    %impliments a cutoff to ensure that the temporal roi coordinates are
    %anterior of the y=-15 plane
    anteriorPoints=find(mergedLatTempROI.coords(:,2)>-15);
    mergedLatTempROI.coords=mergedLatTempROI.coords(anteriorPoints ,: );
    
    
    %% not ROI
    %excludes the lingual gyrus which is too lateral for the MdLF
    % 122=oc-temp_med-Lingual
    
    [mergedNotROI] =bsc_roiFromFSnums(fsDir,[122]+sidenum,1,smoothParameter);
    mergedNotROI.name='NotROI'; 

    %% segmenting
    
    %set operands for ROIS
    operands={'and','and', 'not'};
    
    %switch for correct name
    if leftright==2
        sideflag='R';
    else
        sideflag='L';
    end
    currentFascicleName=strcat(sideflag,'_MdLF');
    
    %create object containing all rois
    currentROIs= [{mergedLatTempROI} {mergedOCTPROI} {mergedNotROI}];
    
    %actually segment
    [fascicle, FiberBoolVec] = feSegmentFascicleFromConnectome(wbfg, currentROIs, operands, currentFascicleName);
    
    %obtain fiber indexes corresponding to MdLF
    FiberIndexes=find(FiberBoolVec);
    
    %directs segmentation output to correct function output holder
    if leftright == 2
        RightMdLF=fascicle;
        RightMdLFindexes=FiberIndexes;
    else
        LeftMdLF=fascicle;
        LeftMdLFindexes=FiberIndexes;
    end
end
end
