function [classification] =bsc_opticRadiationSeg_V4(wbfg, fsDir)
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
for iFibers=1:length(wbfg.fibers)
    fiberNodeNum=round(length(wbfg.fibers{iFibers})/2);
    curStreamline=wbfg.fibers{iFibers};
    midpoints(iFibers,:)=curStreamline(:,fiberNodeNum);
    endpoints1(iFibers,:)=curStreamline(:,1);
    endpoints2(iFibers,:)=curStreamline(:,end);
    %streamLengths(iFibers)=sum(sqrt(sum(diff(wbfg.fibers{iFibers},1,2).^2)));
end
%maybe play with this if something isn't to your liking, it corresponds to
%the smoothing kernel used for the parietal and temporal ROIs
smoothParameter=5;
thalamusLut=[10 49];
choroidLut=[31 63];
latVentLut=[4 43];
wmLut=[2 41];
hippocamLut=[17 53];
classification=[];
 classification.index=zeros(length(wbfg.fibers),1);
     
   
    classification.names={'Left Meyers Loop' , 'Right Meyers Loop', 'Left Baums Loop', 'Right Baums Loop', ...
        'Left Optic Unclassified', 'Right Optic Unclassified', 'Left Optic Accesory', 'Right Optic Accesory'};
    

%iterates through left and right sides
for leftright= [1,2]
    
    %sidenum is basically a way of switching  between the left and right
    %hemispheres of the brain in accordance with freesurfer's ROI
    %numbering scheme. left = 1, right = 2
    sidenum=10000+leftright*1000;
    
    
    
    %% thalamic ROI
    % we can be generous here, as we will be limiting our selection in
    % subsequent steps
    
    
    
    [thalamicROI] =wma_roiFromFSnums(fsDir,thalamusLut(leftright),0);
    thalamicROI.name='thalamicROI';
    
    anteriorThalamicPoint=max(thalamicROI.coords(:,2));
    latThalBorder=bsc_planeFromROI(thalamicROI,'lateral',fsDir);
    lateralThalamicPoint=unique(latThalBorder.coords(:,1));
    postThalamicPoint=min(thalamicROI.coords(:,2));
    anteriorPoints=find(thalamicROI.coords(:,2)<anteriorThalamicPoint-5);
    topThalamicPoint=max(thalamicROI.coords(:,3));
    thalamicMidAxial=mean(unique(thalamicROI.coords(:,1)));
    
    inferiorThalNot=bsc_planeFromROI(thalamicROI,'inferior',fsDir);
    [thalamicROI] =wma_roiFromFSnums(fsDir,thalamusLut(leftright),1,7);
    strangePrevent=bsc_planeFromROI(thalamicROI,'posterior',fsDir);
    strangePrevent.coords=strangePrevent.coords( strangePrevent.coords(:,3)<topThalamicPoint &  abs(strangePrevent.coords(:,1))<abs(lateralThalamicPoint)-5,:);
    hippocampalNOT=bsc_planeFromROI(hippocamLut(leftright),'posterior',fsDir);
    hippocampalNOT.coords=hippocampalNOT.coords(   abs(hippocampalNOT.coords(:,1))<abs(lateralThalamicPoint),:);
    
    %thalamicROI.coords=thalamicROI.coords(anteriorPoints ,: );
    
    
    
    
    %% posterior Occpital
    %generates the roi posterior occipital regions
    
    
    %119
    fullOccpROI=bsc_roiFromFSnums(fsDir,[120 119 111 158 166 143 145 159 152 122 162 161 121]+sidenum,1,9);
    occSupPlane=bsc_planeFromROI( [120]+sidenum,'anterior',fsDir);
    calcerineROI=bsc_roiFromFSnums(fsDir,[145]+sidenum,1,9);
    calcerineROI.coords=calcerineROI.coords(calcerineROI.coords(:,2)<unique(occSupPlane.coords(:,2)),:);
    
    %fullOccpROI=bsc_mergeROIs(fullOccpROI,calcerineROI);
    
    
    
    
    postOccpROI =bsc_roiFromFSnums(fsDir,[111 122 120]+sidenum,1,9);
    postOccpROI=bsc_mergeROIs(postOccpROI,calcerineROI);
    postOccpROI.coords=postOccpROI.coords(postOccpROI.coords(:,2)<unique(occSupPlane.coords(:,2)),:);
    
    
    postOccpROI.name='posteriorOccipitalROI';
    
    
    
    
    %% not ROIs
    % Make planar ROI 10 mm above the top voxel of the thalamic ROI
    
    atlasNifti = wma_getAsegFile(fsDir , '2009');
    
    [topNotROI] =bsc_makePlanarROI(atlasNifti,topThalamicPoint+10,'z');
    topNotROI.name='topNotROI';
    
    % posterior bound of thalamus
    postTopExclude=bsc_planeFromROI(thalamicROI,'posterior',fsDir);
    postTopExclude.coords= postTopExclude.coords( postTopExclude.coords(:,3)>topThalamicPoint,:);
    
    %interMeyerNot
    interMeyerExclude=bsc_makePlanarROI(atlasNifti,(thalamicMidAxial+lateralThalamicPoint)/2,'x');
    interMeyerExclude.coords= interMeyerExclude.coords( interMeyerExclude.coords(:,2)<postThalamicPoint &  interMeyerExclude.coords(:,2)>unique(occSupPlane.coords(:,2))+10,:);
    
    interMeyerExcludeHalf=interMeyerExclude;
    interMeyerExcludeHalf.coords=interMeyerExcludeHalf.coords( interMeyerExcludeHalf.coords(:,2) >mean(unique(interMeyerExcludeHalf.coords(:,2))),:);
    
    
    % prevent cross hemispheric fibers
    [midsaggitalNot] =bsc_makePlanarROI(atlasNifti,0,'x');
    midsaggitalNot.name='midsaggitalNot';
    
    %prevent anterior stretching fibers
    
    preventOccip=bsc_planeFromROI(120+sidenum,'anterior',fsDir);
    
    
    [antNotROI] =bsc_makePlanarROI(atlasNifti,anteriorThalamicPoint-15,'y');
    antNotROI.name='antNotROI';
    
    periPost=bsc_planeFromROI(167 + sidenum,'posterior',fsDir);
    postCingInf=bsc_planeFromROI(109 + sidenum,'inferior',fsDir);
    periPostUpper=periPost;
    periPostLower=periPost;
    
    periPostUpper.coords=periPostUpper.coords( periPostUpper.coords(:,3) >unique(postCingInf.coords(:,3)),:);
    periPostLower.coords=periPostLower.coords( periPostLower.coords(:,3) <unique(postCingInf.coords(:,3))-2,:);
    
    postCingInf.coords=postCingInf.coords( postCingInf.coords(:,2) <unique(periPostLower.coords(:,2))-6,:);
    
    
    
    
    %% other ROIS
    
    
    
    %% segmenting
    
    %set operands for ROIS
    
    
    %switch for correct name
    if leftright==2
        sideflag='R';
    else
        sideflag='L';
    end
    currentFascicleName=strcat(sideflag,'_OR');
    
    %create object containing all rois
    
    
    %actually segment
    [fascicle, FiberBoolVec] = wma_SegmentFascicleFromConnectome(wbfg, [{postOccpROI} {thalamicROI} {midsaggitalNot} {antNotROI} {postTopExclude} {interMeyerExclude} {strangePrevent} {hippocampalNOT}], {'endpoints','endpoints' 'not', 'not' , 'not','not','not','not'}, currentFascicleName);
    
    [test, testIND] = wma_SegmentFascicleFromConnectome(wbfg, [{fullOccpROI} {thalamicROI} {midsaggitalNot} {antNotROI} {preventOccip} {strangePrevent} {postTopExclude} {periPostUpper} {periPostLower} {postCingInf}], {'endpoints','endpoints' 'not', 'not', 'not','not' ,'not','and','not','not'}, currentFascicleName);
    
    critearia1=or(endpoints1(:,2)<postThalamicPoint & abs(endpoints1(:,1))<abs(lateralThalamicPoint),endpoints2(:,2)<postThalamicPoint & abs(endpoints2(:,1))<abs(lateralThalamicPoint));
    %critearia2=or(endpoints1(:,3)<topThalamicPoint & endpoints1(:,2)>unique(occSupPlane.coords(:,2)) & endpoints2(:,2)<thalPostCoord, endpoints2(:,3)<topThalamicPoint & endpoints2(:,2)>unique(occSupPlane.coords(:,2)) & endpoints1(:,2)<thalPostCoord);
    %critearia2=or(endpoints1(:,3)<midpoints(:,3) & endpoints2(:,3)>midpoints(:,3) , endpoints2(:,3)<midpoints(:,3) & endpoints1(:,3)>midpoints(:,3));
    critearia2=~(endpoints1(:,2)<postThalamicPoint & endpoints2(:,2)<postThalamicPoint);
    
    
    [test2, test2IND]=wma_SegmentFascicleFromConnectome(wbfg, [{fullOccpROI} {thalamicROI} {midsaggitalNot} {antNotROI} {interMeyerExclude} { occSupPlane} {strangePrevent} {postTopExclude} {hippocampalNOT}], {'endpoints','endpoints' 'not', 'not','not', 'and', 'not', 'not','not'}, currentFascicleName);
    % critearia3=or(endpoints1(:,2)<thalPostCoord-10 & abs(endpoints1(:,1))<abs(lateralThalamicPoint) & endpoints1(:,3)<topThalamicPoint,endpoints2(:,2)<thalPostCoord-10 & abs(endpoints2(:,1))<abs(lateralThalamicPoint) & endpoints2(:,3)<topThalamicPoint);
    critearia7=~( endpoints1(:,2)<midpoints(:,2) & endpoints2(:,2)<midpoints(:,2)) ;
    
    
    [test3, test3IND]=wma_SegmentFascicleFromConnectome(wbfg, [{fullOccpROI} {thalamicROI} {midsaggitalNot} {antNotROI} {postCingInf}  {periPostUpper}  {periPostLower} {inferiorThalNot} {interMeyerExclude}], {'endpoints','endpoints' 'not', 'not', 'not','not','and','not','and'}, currentFascicleName);
    critearia4= midpoints(:,2)<postThalamicPoint;
    critearia5= abs(endpoints1(:,1))<abs(midpoints(:,1)) & abs(endpoints2(:,1))<abs(midpoints(:,1)) ;
    critearia6= endpoints1(:,3)<unique(postCingInf.coords(:,3)) & endpoints2(:,3)<unique(postCingInf.coords(:,3)) ;
    
    test3IND=test3IND & ~test2IND & ~FiberBoolVec & ~testIND &critearia4 &  critearia5 &   critearia6;
    
    
    %test2IND=test2IND & critearia4 & critearia3 & critearia5;
    
    
    testIND=testIND & ~FiberBoolVec & critearia1 &critearia2 ;
    
    test2IND=test2IND & ~FiberBoolVec & ~testIND & critearia7 ;
    
    %     test.fibers=wbfg.fibers(testIND  );
    %      test2.fibers=wbfg.fibers(  test2IND );
    %         test3.fibers=wbfg.fibers(  test3IND );
    
    
    
    


    
    if leftright==1
        classification.index(FiberBoolVec)=find( strcmp(classification.names,'Left Meyers Loop'));
        classification.index(testIND)=find( strcmp(classification.names,'Left Baums Loop'));
        
        classification.index(test2IND)=find( strcmp(classification.names,'Left Optic Unclassified'));
        classification.index(test3IND)=find( strcmp(classification.names,'Left Optic Accesory'));
        
    else
        
        classification.index(FiberBoolVec)=find( strcmp(classification.names,'Right Meyers Loop'));
        classification.index(testIND)=find( strcmp(classification.names,'Right Baums Loop'));
        
        classification.index(test2IND)=find( strcmp(classification.names,'Right Optic Unclassified'));
        classification.index(test3IND)=find( strcmp(classification.names,'Right Optic Accesory'));
        
        
        
    end
    
end


end