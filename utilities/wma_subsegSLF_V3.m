function [LeftSLF1, LeftSLF1Bool, RightSLF1, RightSLF1Bool, LeftSLF2, LeftSLF2Bool, RightSLF2, RightSLF2Bool,...
    LeftSLF3, LeftSLF3Bool, RightSLF3, RightSLF3Bool] =wma_subsegSLF_V3(wbfg, fsDir)
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
%   fg structures and boolean vectors for the left and right variants of
%   each of the slf subcomponents
%
%
%  Same for the other tracts
% (C) Daniel Bullock, 2017, Indiana University



%%  compute relevant stats for all tracts
for iFibers=1:length(wbfg.fibers)
    fiberNodeNum=round(length(wbfg.fibers{iFibers})/2);
    curStreamline=wbfg.fibers{iFibers};
    midpoints(iFibers,:)=curStreamline(:,fiberNodeNum);
    endpoints1(iFibers,:)=curStreamline(:,1);
    endpoints2(iFibers,:)=curStreamline(:,end);
    streamLengths(iFibers)=sum(sqrt(sum(diff(wbfg.fibers{iFibers},1,2).^2)));
end

leftBool=midpoints(:,1)<0;
lengths=hist(streamLengths,1:ceil(max(streamLengths)));
lengthBool=streamLengths>80 & streamLengths<150;

labelNifti= wma_getAsegFile(fsDir , '2009');

[Plane1]=bsc_makePlanarROI(labelNifti,0, 'y');
[Plane2]=bsc_makePlanarROI(labelNifti,-10, 'y');
[Plane3]=bsc_makePlanarROI(labelNifti,10, 'y');

    putCaudNot=bsc_roiFromFSnums(fsDir,[11, 12 , 50 ,51, 31,63, 17,53, ],1,11);%12145,11145, 11162, 12162,11161, 12161, 11102,12102
    cingNot=bsc_roiFromFSnums(fsDir,[11167,12167],1,11);

%% SLF1
for isides=[11000,12000]
%     SLF1FrontalEndpointROI=bsc_roiFromFSnums(fsDir,isides+[106,116],1,3);
%     SLF1PosteriorEndpointROI=bsc_roiFromFSnums(fsDir,isides+[130,127,157],1,3);
%     
%     [fg, keep]=bsc_tractByEndpointROIs(wbfg,{SLF1FrontalEndpointROI, SLF1PosteriorEndpointROI});
%     
%     

    noCross=bsc_makePlanarROI(labelNifti,0, 'x');
    operands={'and','and','and','not','not','not'};
    currentROIs= [{Plane1} {Plane2} {Plane3} {noCross} {putCaudNot} {cingNot}];
    
    [fascicle1, FiberBoolVec1] = wma_SegmentFascicleFromConnectome(wbfg, currentROIs, operands, 'slf');
    
    
    
    slf1MedRoi=wma_roiFromFSnums(fsDir,[isides+154], 0);
    
    
    if isides==11000
        meanFrontX=max(slf1MedRoi.coords(:,1));
        LeftSLF1=fascicle1;
        LeftSLF1.fibers=wbfg.fibers(FiberBoolVec1 & midpoints(:,1)<0 & midpoints(:,1)>meanFrontX);
        LeftSLF1Bool=FiberBoolVec1 & midpoints(:,1)<0 & midpoints(:,1)>meanFrontX;
    else
        meanFrontX=min(slf1MedRoi.coords(:,1));
        RightSLF1=fascicle1;
        RightSLF1.fibers=wbfg.fibers(FiberBoolVec1 & midpoints(:,1)>0 & midpoints(:,1)<meanFrontX);
        RightSLF1Bool=FiberBoolVec1 & midpoints(:,1)>0 & midpoints(:,1)<meanFrontX;
    end
end
    
%% slf2
for isides=[11000,12000]
    SLF2FrontalEndpointROI=bsc_roiFromFSnums(fsDir,isides+[115],1,7);
    SLF2PosteriorEndpointROI=bsc_roiFromFSnums(fsDir,isides+[125,156],1,7);
    
    if isides==11000
        slf2AntLimitROI=bsc_roiFromFSnums(fsDir,[12],1,3);
    else
        slf2AntLimitROI=bsc_roiFromFSnums(fsDir,[51],1,3);
    end
    slf2AntLimit=max(slf2AntLimitROI.coords(:,2));
    SLF2FrontalEndpointROI.coords=SLF2FrontalEndpointROI.coords(SLF2FrontalEndpointROI.coords(:,2)<slf2AntLimit,:);
    
    [fg, keep]=bsc_tractByEndpointROIs(wbfg,{SLF2FrontalEndpointROI, SLF2PosteriorEndpointROI});
    
    if isides==11000
        LeftSLF2=fg;
        LeftSLF2Bool=keep;
    else
        RightSLF2=fg;
        RightSLF2Bool=keep;
    end
    
    
end



%% slf3

for isides=[11000,12000]
    SLF3FrontalEndpointROI=bsc_roiFromFSnums(fsDir,isides+[112],1,3);
    SLF3PosteriorEndpointROI=bsc_roiFromFSnums(fsDir,isides+[126],1,1);
    
    
    [fg, keep]=bsc_tractByEndpointROIs(wbfg,{SLF3FrontalEndpointROI, SLF3PosteriorEndpointROI});
    
 
    
    if isides==11000
        LeftSLF3=fg;
        LeftSLF3Bool=keep;
    else
        RightSLF3=fg;
        RightSLF3Bool=keep;
    end
    
    
end


   
    
    

%% Cingulum Cingulate
    putCaudNot=bsc_roiFromFSnums(fsDir,[10,49,51,12 ],1,3);

% for isides=[11000,12000]
%     
%     
%     if isides==11000
%         thallROI=bsc_roiFromFSnums(fsDir,[10],0);
%         ThalLatBound=min (thallROI.coords(:,1));
%         latPlane=bsc_makePlanarROI(labelNifti,ThalLatBound, 'x');
%     else
%         thallROI=bsc_roiFromFSnums(fsDir,[49],0);
%         ThalLatBound=max (thallROI.coords(:,1));
%         latPlane=bsc_makePlanarROI(labelNifti,ThalLatBound, 'x');
%     end
%     
%     
%     
%     cing=bsc_roiFromFSnums(fsDir,[isides+167],1,5);
%     operands={'and','and','and','not','not','and', 'not'};
%     currentROIs= [{Plane1} {Plane2} {Plane3} {noCross} {putCaudNot } {cing} {latPlane}];
%     
%     if isides==11000
%         [LeftCing, LeftCingBoolVec] = wma_SegmentFascicleFromConnectome(wbfg, currentROIs, operands, 'slf');
%     else
%         [RightCing, RightCingBoolVec] = wma_SegmentFascicleFromConnectome(wbfg, currentROIs, operands, 'slf');
%     end
% 
% 
% 
% 
% end
% 
% 
% 







% 
%  
% 
% figure
% bsc_quickPlot(LeftCing)
% bsc_quickPlot(RightCing)
% bsc_quickPlot(RightSLF1)
% %bsc_quickPlot(fascicle1)
% bsc_quickPlot(LeftSLF1)
% bsc_quickPlot(RightSLF1)
% bsc_quickPlot(LeftSLF2)
% bsc_quickPlot(RightSLF2)
% bsc_quickPlot(LeftSLF3)
% bsc_quickPlot(RightSLF3)
% bsc_quickPlot(LeftSLF6)
% bsc_quickPlot(RightSLF6)
% bsc_quickPlot(LeftSLF7)
% bsc_quickPlot(RightSLF7)
% view([90,0])
% %figure
% 
% testRoi=wma_roiFromFSnums(fsDir,[isides]+153, 0);
% 
% 
% figure
% bsc_plotROIEndpointsOnFG(RightSLF3, testRoi)
% 
% 
% figure
% bsc_plotEndpointsOnFG(LeftSLF3)




end

