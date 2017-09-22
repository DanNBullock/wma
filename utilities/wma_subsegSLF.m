function [LeftSlf1, LeftSlf1Bool, RightSlf1, RightSlf1Bool, LeftSlf2, LeftSlf2Bool, RightSlf2, RightSlf2Bool,...
    LeftSlf3, LeftSlf3Bool, RightSlf3, RightSlf3Bool] =wma_subsegSLF(wbfg, fsDir)
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



%% SLF1

putCaudNot=bsc_roiFromFSnums(fsDir,[11, 12 , 50 ,51, 31,63, 17,53, ],1,13);%12145,11145, 11162, 12162,11161, 12161, 11102,12102
noCross=bsc_makePlanarROI(labelNifti,0, 'x');
operands={'and','and','and','not','not'};
currentROIs= [{Plane1} {Plane2} {Plane3} {noCross} {putCaudNot}];

[fascicle1, FiberBoolVec1] = wma_SegmentFascicleFromConnectome(wbfg, currentROIs, operands, 'slf');

%medial border of 149
slf1XCriterion=abs(midpoints(:,1))<30;

slf1YCriterion=(midpoints(:,2))>-30 & (midpoints(:,2))<0;
slf1ZCriterion=(midpoints(:,3))>35; %45
%slf1EndpointsCriterion=or(endpoints1(:,2)<-30,  endpoints2(:,2)<-30) &  and(or(endpoints1(:,3)>25,endpoints1(:,2)>20),or(endpoints2(:,3)>25,endpoints2(:,2)>20));

fascicle1.fibers=wbfg.fibers(FiberBoolVec1 & slf1XCriterion & slf1YCriterion & slf1ZCriterion & lengthBool');% & slf1EndpointsCriterion)   ;% slf1ZCriterion
fascicle1Bool=FiberBoolVec1 & slf1XCriterion & slf1YCriterion & slf1ZCriterion & lengthBool';

%% slf2

slf2XCriterion=abs(midpoints(:,1))<35 & abs(midpoints(:,1))>25;
slf2YCriterion=(midpoints(:,2))>-30 & (midpoints(:,2))<0;
slf2ZCriterion=(midpoints(:,3))<40 & (midpoints(:,3))>20;                   
slf2EndpointsCriterion=and(endpoints1(:,3)>10,endpoints2(:,3)>10);
slf2EndpointsCriterion=slf2EndpointsCriterion & or(and(endpoints1(:,2)<0,endpoints1(:,2)<-50),and(endpoints2(:,2)<0,endpoints2(:,2)<-50));

fascicle2.fibers=wbfg.fibers(FiberBoolVec1 & slf2XCriterion & slf2YCriterion & slf2ZCriterion & slf2EndpointsCriterion);

%% slf3 

[Planea3]=bsc_makePlanarROI(labelNifti,35, 'z');
Planea3.coords=Planea3.coords(abs(Planea3.coords(:,1))<40,:);

[~, excludeBoolVec3] = wma_SegmentFascicleFromConnectome(wbfg, {Planea3}, {'and'}, 'exclusion');

slf3XCriterion=abs(midpoints(:,1))>30;

slf3YCriterion=(midpoints(:,2))>-15 & (midpoints(:,2))<15;

slf3ZCriterion=(midpoints(:,3))<25 & (midpoints(:,3))>5;

%Y criterion: anterior posterior cutoff
slf3EndpointsCriterion=and(and(endpoints1(:,2)>-40,endpoints1(:,2)<45),and(endpoints2(:,2)>-40,endpoints2(:,2)<45));

%z criterion: no hi:hi or lo:lo relative to midpoint
slf3EndpointsCriterion=slf3EndpointsCriterion & ~or(and((midpoints(:,3)-endpoints1(:,3))>0,(midpoints(:,3)-endpoints2(:,3))>0),and((midpoints(:,3)-endpoints1(:,3))<0,(midpoints(:,3)-endpoints2(:,3))<0));

%anterior endpoints should be no more than 5 mm higher than posterior ones
slf3EndpointsCriterion=slf3EndpointsCriterion & or(and(endpoints1(:,2)<0,endpoints1(:,3)-endpoints2(:,3)>-5),and(endpoints2(:,2)<0,endpoints2(:,3)-endpoints1(:,3)>-5));

%remove fibers that go too medial after the 0 y plane
slf3EndpointsCriterion=slf3EndpointsCriterion & or(and(endpoints1(:,2)<0,abs(endpoints1(:,1))>abs(midpoints(:,1))-5),and(endpoints2(:,2)<0,abs(endpoints2(:,1))>abs(midpoints(:,1))-5));

%Left_putamen exclusion

putNot=bsc_roiFromFSnums(fsDir,[12, 51],1,3);

noCross=bsc_makePlanarROI(labelNifti,0, 'x');

operands={'and','and','and','not','not'};
currentROIs= [{Plane1} {Plane2} {Plane3},{putNot},{noCross}];


[fascicle3, FiberBoolVec3] = wma_SegmentFascicleFromConnectome(wbfg, currentROIs, operands, 'slf3');

fascicle3.fibers=wbfg.fibers(FiberBoolVec3 & slf3XCriterion & slf3YCriterion & slf3ZCriterion & slf3EndpointsCriterion & ~excludeBoolVec3);

fascicle3Bool=FiberBoolVec3 & slf3XCriterion & slf3YCriterion & slf3ZCriterion & slf3EndpointsCriterion & ~excludeBoolVec3;

%% cleanup
LeftSlf1=fascicle1;
RightSlf1=fascicle1;
LeftSlf1.fibers=wbfg.fibers(FiberBoolVec1 & slf1XCriterion & slf1YCriterion & slf1ZCriterion & lengthBool'& leftBool);
RightSlf1.fibers=wbfg.fibers(FiberBoolVec1 & slf1XCriterion & slf1YCriterion & slf1ZCriterion & lengthBool'& ~leftBool);
LeftSlf1Bool=FiberBoolVec1 & slf1XCriterion & slf1YCriterion & slf1ZCriterion & lengthBool'& leftBool;
RightSlf1Bool=FiberBoolVec1 & slf1XCriterion & slf1YCriterion & slf1ZCriterion & lengthBool'& ~leftBool;

LeftSlf2=fascicle2;
RightSlf2=fascicle2;
LeftSlf2.fibers=wbfg.fibers(FiberBoolVec1 & slf2XCriterion & slf2YCriterion & slf2ZCriterion & slf2EndpointsCriterion & leftBool);
RightSlf2.fibers=wbfg.fibers(FiberBoolVec1 & slf2XCriterion & slf2YCriterion & slf2ZCriterion & slf2EndpointsCriterion & ~leftBool);
LeftSlf2Bool=FiberBoolVec1 & slf1XCriterion & slf1YCriterion & slf1ZCriterion & lengthBool'& leftBool;
RightSlf2Bool=FiberBoolVec1 & slf1XCriterion & slf1YCriterion & slf1ZCriterion & lengthBool'& ~leftBool;

LeftSlf3=fascicle3;
RightSlf3=fascicle3;
LeftSlf3.fibers=wbfg.fibers(FiberBoolVec3 & slf3XCriterion & slf3YCriterion & slf3ZCriterion & slf3EndpointsCriterion & ~excludeBoolVec3 & leftBool);
RightSlf3.fibers=wbfg.fibers(FiberBoolVec3 & slf3XCriterion & slf3YCriterion & slf3ZCriterion & slf3EndpointsCriterion & ~excludeBoolVec3 & ~leftBool);
LeftSlf3Bool=FiberBoolVec3 & slf3XCriterion & slf3YCriterion & slf3ZCriterion & slf3EndpointsCriterion & ~excludeBoolVec3;
RightSlf3Bool=FiberBoolVec3 & slf3XCriterion & slf3YCriterion & slf3ZCriterion & slf3EndpointsCriterion & ~excludeBoolVec3;
end

