function [classification] =wma_segMajTracks_v4(wbfg, fsDir)
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
for iFibers=1:length(wbfg.fibers)
    fiberNodeNum=round(length(wbfg.fibers{iFibers})/2);
    curStreamline=wbfg.fibers{iFibers};
    midpoints(iFibers,:)=curStreamline(:,fiberNodeNum);
    endpoints1(iFibers,:)=curStreamline(:,1);
    endpoints2(iFibers,:)=curStreamline(:,end);
    %streamLengths(iFibers)=sum(sqrt(sum(diff(wbfg.fibers{iFibers},1,2).^2)));
end

fprintf('\n Midpoints and endpoints computed')

wmLut=[2 41];
putLut=[12 51];
hippocampusLut=[17 53];
caudateLut=[11 50];
amigdalaLut=[18 54];
amigLut=[12 18; 51 54];
thalamusLut=[10 49];
lentiLut=[12 13; 51 52];
vesselLut=[30 62];
ventricleLut=[4 43];

labelNifti= wma_getAsegFile(fsDir , '2009');
LeftStreams=midpoints(:,1)<0;

classification=[];

classification.index=zeros(length(wbfg.fibers),1);
classification.names={'Corpus Callosum' , 'Forceps Major', 'Forceps Minor', 'Left IFOF', 'Right IFOF', 'Left Cingulum1', ...
    'Right Cingulum1', 'Left Cingulum2', 'Right Cingulum2', 'Left Cingulum3', 'Right Cingulum3', ...
    'Left SFOF', 'Right SFOF', 'Left Uncinate', 'Right Uncinate',...
    'Left Arcuate','Right Arcuate', 'Left cortico-spinal','Right cortico-spinal', ...
    'Left anterior thalamic','Right anterior thalamic','Left inferior thalamic','Right inferior thalamic',...
    'Left superior thalamic','Right superior thalamic','Left posterior thalamic','Right posterior thalamic', ...
    'Left VOF', 'Right VOF', 'Left AntVOF', 'Right AntVOF', 'Left VisionU', 'Right VisionU','Left OccipitalU', 'Right OccipitalU',...
    'Left Thalamico-spinal','Right Thalamico-spinal'};







%% non-lateral ROIs, planes and coordinates

%ROIs
stemnot=bsc_roiFromFSnums(fsDir,[28 60 16],0);
cerebellumRoi=bsc_roiFromFSnums(fsDir,[47,8],0);
stemroi=bsc_roiFromFSnums(fsDir,[16],0);
thalTotal=bsc_roiFromFSnums(fsDir,[10 49],0);
AnteriorCCROI=bsc_roiFromFSnums(fsDir,[255],0);
periC2ROI=bsc_roiFromFSnums(fsDir,[12167 11167],0);
preCuS2ROI=bsc_roiFromFSnums(fsDir,[12130 11130],0);
PosteriorCCROI1=bsc_roiFromFSnums(fsDir,[ 251],0);
PosteriorCCROI2=bsc_roiFromFSnums(fsDir,[ 251 252],0);
parietalROI=bsc_roiFromFSnums(fsDir,[121230 11120],0);
csf=bsc_roiFromFSnums(fsDir,[24],0);
arcInt=bsc_planeFromROI(thalTotal,'superior', fsDir);

%coordinates
stemant=max(stemroi.coords(:,2));
cerebellumTop=max(cerebellumRoi.coords(:,3));
csfMid=median(csf.coords(:,2));
csfTop=max(csf.coords(:,3));
bstemAnt=max(stemnot.coords(:,2));
AnteriorCCcutoffP=min(AnteriorCCROI.coords(:,2));
AnteriorCCcutoffA=max(periC2ROI.coords(:,2));

%planes
%[superiorCutPlane] = bsc_planeFromROI(bstemAnt, 'anterior',fsDir)
superiorCutPlane=bsc_makePlanarROI(labelNifti,bstemAnt, 'y');
superiorCutPlane.coords=superiorCutPlane.coords(superiorCutPlane.coords(:,3)>cerebellumTop,:);
superiorCutPlane2=bsc_makePlanarROI(labelNifti,stemant, 'y');
superiorCutPlane2.coords=superiorCutPlane2.coords(superiorCutPlane2.coords(:,3)>cerebellumTop,:);
superiorCutPlane3=superiorCutPlane;
superiorCutPlane3.coords=vertcat(superiorCutPlane.coords(superiorCutPlane.coords(:,3)>cerebellumTop,:),superiorCutPlane2.coords(superiorCutPlane2.coords(:,3)>cerebellumTop,:));
transvPlane=bsc_planeFromROI(thalTotal,'superior', fsDir);
thaltop=unique(transvPlane.coords(:,3));
interHemisphere=bsc_makePlanarROI(labelNifti,0, 'x');

AnteriorCCPlaneP=bsc_makePlanarROI(labelNifti,AnteriorCCcutoffP -5, 'y');
AnteriorCCPlaneA=bsc_makePlanarROI(labelNifti,AnteriorCCcutoffA , 'y');
AnteriorCCPlaneAL=AnteriorCCPlaneA;
AnteriorCCPlaneAL.coords=AnteriorCCPlaneAL.coords(AnteriorCCPlaneAL.coords(:,1)<-3,:);
AnteriorCCPlaneAR=AnteriorCCPlaneA;
AnteriorCCPlaneAR.coords=AnteriorCCPlaneAR.coords(AnteriorCCPlaneAR.coords(:,1)>3,:);












[AnteriorCC, AnteriorCCIND]=wma_SegmentFascicleFromConnectome(wbfg, [  {AnteriorCCROI} {AnteriorCCPlaneP} {AnteriorCCPlaneAL} {AnteriorCCPlaneAR}], { 'and' , 'not' , 'and', 'and'}, 'AnteriorCC');
AnteriorEndpointCriteria=endpoints1(:,2) > AnteriorCCcutoffP+(AnteriorCCcutoffA-AnteriorCCcutoffP)/1.5  & endpoints2(:,2) > AnteriorCCcutoffP+(AnteriorCCcutoffA-AnteriorCCcutoffP)/1.5 ;
AnteriorCCIND=AnteriorCCIND & AnteriorEndpointCriteria;
%AnteriorCCm=AnteriorCC;

AnteriorCC.fibers=wbfg.fibers(AnteriorCCIND & AnteriorEndpointCriteria,:);
%AnteriorCCm.fibers=wbfg.fibers(AnteriorCCIND & ~AnteriorEndpointCriteria,:);




PosteriorCCcutoffP=min(preCuS2ROI.coords(:,2));
PosteriorCCcutoffA=max(PosteriorCCROI1.coords(:,2));

parietalAnt=max(parietalROI.coords(:,2));
PosteriorCCPlaneP=bsc_makePlanarROI(labelNifti,PosteriorCCcutoffP , 'y');
PosteriorCCPlanePL=PosteriorCCPlaneP;
PosteriorCCPlanePL.coords=PosteriorCCPlanePL.coords(PosteriorCCPlanePL.coords(:,1)<-3,:);
PosteriorCCPlanePR=PosteriorCCPlaneP;
PosteriorCCPlanePR.coords=PosteriorCCPlanePR.coords(PosteriorCCPlanePR.coords(:,1)>3,:);
PosteriorCCPlaneA=bsc_makePlanarROI(labelNifti,PosteriorCCcutoffA +5, 'y');
%[PosteriorCC, PosteriorCCIND]=wma_SegmentFascicleFromConnectome(wbfg, [  {PosteriorCCROI} {PosteriorCCPlaneA} {PosteriorCCPlanePL} {PosteriorCCPlanePR}], { 'and' , 'not','and' , 'and', }, 'AnteriorCC');
[PosteriorCC, PosteriorCCIND]=wma_SegmentFascicleFromConnectome(wbfg, [  {PosteriorCCROI2} {PosteriorCCPlaneA} ], { 'and' , 'not'}, 'AnteriorCC');

%posteriorEndpointCriteria=endpoints1(:,2) < PosteriorCCcutoffA+ (PosteriorCCcutoffP-PosteriorCCcutoffA)/1.5  & endpoints2(:,2) < PosteriorCCcutoffA +(PosteriorCCcutoffP-PosteriorCCcutoffA)/1.5 ;
%posteriorEndpointCriteria=or(endpoints1(:,2) < parietalAnt  & endpoints2(:,2) < PosteriorCCcutoffA +(PosteriorCCcutoffP-PosteriorCCcutoffA)/1.25,...
%endpoints1(:,2) < PosteriorCCcutoffA+ (PosteriorCCcutoffP-PosteriorCCcutoffA)/1.25  & endpoints2(:,2) < parietalAnt);


posteriorEndpointCriteria=endpoints1(:,2) < parietalAnt  & endpoints2(:,2) < parietalAnt ;
%posteriorEndpointCriteria=or(endpoints1(:,2) < parietalAnt  , endpoints2(:,2) < parietalAnt) ;


PosteriorCCIND=PosteriorCCIND & posteriorEndpointCriteria;
%PosteriorCCm=PosteriorCC;
PosteriorCC.fibers=wbfg.fibers(PosteriorCCIND & posteriorEndpointCriteria,:);
%PosteriorCCm.fibers=wbfg.fibers(PosteriorCCIND & ~posteriorEndpointCriteria,:);


BodyCCROI=bsc_roiFromFSnums(fsDir,[254 253 252],0);
CCant=max(BodyCCROI.coords(:,2));
CCpost=min(BodyCCROI.coords(:,2));
infCCcutoff=min(BodyCCROI.coords(:,3));
infCCPlane=bsc_makePlanarROI(labelNifti,infCCcutoff , 'z');
BodyCCtop=max(BodyCCROI.coords(:,3));

latROIs=bsc_roiFromFSnums(fsDir,[12141 11141 11150 12150 11104 12104],1,15);

subparROI=bsc_roiFromFSnums(fsDir,[11172 12172],0);
subparCutPost=min(subparROI.coords(:,2));
posteriorPlane=bsc_makePlanarROI(labelNifti,subparCutPost, 'y');


%[BodyCC, BodyCCIND]=wma_SegmentFascicleFromConnectome(wbfg, [  {BodyCCROI} {stemnot}  {infCCPlane} {AnteriorCCPlaneA} {posteriorPlane} {latROIs}], { 'and' ,'not' , 'not', 'not', 'not', 'not'}, 'BodyCC');
[BodyCC, BodyCCIND]=wma_SegmentFascicleFromConnectome(wbfg, [  {BodyCCROI} {stemnot}   {latROIs}], { 'and' ,'not' ,  'not'}, 'BodyCC');
CCEndpointCriteria= endpoints1(:,3) > midpoints(:,3) & endpoints2(:,3) > midpoints(:,3) & midpoints(:,2) < CCant & midpoints(:,2) > CCpost & endpoints1(:,3) > infCCcutoff & endpoints2(:,3) > infCCcutoff & ~PosteriorCCIND & ~ AnteriorCCIND;
BodyCCIND=BodyCCIND & CCEndpointCriteria;
%BodyCC.fibers=wbfg.fibers(BodyCCIND);

classification.index(BodyCCIND)=find( strcmp(classification.names,'Corpus Callosum'));
classification.index(AnteriorCCIND)=find( strcmp(classification.names,'Forceps Minor'));
classification.index(PosteriorCCIND)=find( strcmp(classification.names,'Forceps Major'));

fprintf('\n Bihemispheric tracts segmented');

for leftright= [1,2]
    
    
    sidenum=10000+leftright*1000;
    
    
    
    
    
    putROI=bsc_roiFromFSnums(fsDir,putLut(leftright),0);
    putTOP=max(putROI.coords(:,3));
    transverseArcCutPlane=bsc_planeFromROI(251,'inferior', fsDir);
    
    periCROI=bsc_roiFromFSnums(fsDir,[167]+sidenum,1,13);
    %hippocampusROI=bsc_roiFromFSnums(fsDir,[17 18],1,5);
    caudate=bsc_roiFromFSnums(fsDir,caudateLut(leftright),0);
    subparROI=bsc_roiFromFSnums(fsDir,[172]+sidenum,0);
    subparCutPost=min(subparROI.coords(:,2));
    posteriorPlane=bsc_makePlanarROI(labelNifti,subparCutPost, 'y');
    caudateCutAnt=max(caudate.coords(:,2));
    
    
    transverseArcCutPlane.coords=transverseArcCutPlane.coords(transverseArcCutPlane.coords(:,2)<caudateCutAnt,:);
    
    periCantCut=max(periCROI.coords(:,2));
    periCantPlane=bsc_makePlanarROI(labelNifti,periCantCut, 'y');
    %caudatePlaneCut=bsc_makePlanarROI(labelNifti,caudateCutLat, 'x');
    caudatePlaneCut=bsc_planeFromROI(caudateLut(leftright),'lateral',fsDir);
    caudateCutLat=unique(caudatePlaneCut.coords(:,1));
    caudateCutLatUpper=caudatePlaneCut;
    anteriorPlane=bsc_makePlanarROI(labelNifti,caudateCutAnt, 'y');
    
    occLimit=bsc_roiFromFSnums(fsDir,[111]+sidenum,0);
    occLimitVal=max(occLimit.coords(:,2));
    %calcerineROI=bsc_roiFromFSnums(fsDir,[145]+sidenum,1,5);
    
    occipitalROI=bsc_roiFromFSnums(fsDir,[120 119 111 158 166 143 145 159 152 122 162 161 121]+sidenum,1,5);
    %parietalExclude=bsc_roiFromFSnums(fsDir,[120 119 111 158 166 143 145 159 152 122 162 161 121]+sidenum,1,5)
    %occipitalROI.coords=occipitalROI.coords(occipitalROI.coords(:,2)<occLimitVal,:);
    
    
    cingPostROI=bsc_roiFromFSnums(fsDir,[110]+sidenum,21);
    
    upperCaudateCut=max(cingPostROI.coords(:,3));
    antArcCut=max(cingPostROI.coords(:,2));
    antArcPlane=bsc_makePlanarROI(labelNifti,antArcCut, 'y');
    antArcPlane.coords=antArcPlane.coords(antArcPlane.coords(:,3)<thaltop,:);
    
    periCTop=max(periCROI.coords(:,3));
    caudateCutLatUpper.coords=caudateCutLatUpper.coords(caudateCutLatUpper.coords(:,3)>upperCaudateCut,:);
    
    fprintf('\n Initial ROIs created');
    
    
    
    [FOF, FOFIND]=wma_SegmentFascicleFromConnectome(wbfg, [{anteriorPlane} {posteriorPlane} {stemnot} {thalTotal} {superiorCutPlane3} {interHemisphere}], {'and', 'and', 'not', 'not', 'not','not'}, 'FOF');
    IFOFEndpointCriteria=or(endpoints1(:,2) < subparCutPost & endpoints2(:,2) > caudateCutAnt, endpoints2(:,2) < subparCutPost & endpoints1(:,2) > caudateCutAnt) ;
    FOFIND=FOFIND & IFOFEndpointCriteria;
    %FOF.fibers=wbfg.fibers(FOFIND);
    
    ensureArcPlane=bsc_planeFromROI(253,'anterior',fsDir);
    ensureArcLimit=bsc_planeFromROI(253,'superior',fsDir);
     ensureArcPlane.coords=ensureArcPlane.coords(ensureArcPlane.coords(:,3)< unique(ensureArcLimit.coords(:,3)),:);
     
       CingCurveLimit=bsc_planeFromROI(254,'inferior',fsDir);
    
    LeftFOFIND=LeftStreams & FOFIND;
    
    [cingulum1, cingulum1IND]=wma_SegmentFascicleFromConnectome(wbfg, [{periCROI} {periCantPlane} {caudatePlaneCut} {superiorCutPlane} {superiorCutPlane2} {interHemisphere} , {ensureArcPlane}], {'and', 'and','not','and','and', 'not','not'}, 'cingulum1');
    cing1EndpointCriteria=or(endpoints1(:,2) < antArcCut & endpoints2(:,2) > periCantCut, endpoints2(:,2) < antArcCut & endpoints1(:,2) > periCantCut) ;
    cingulum1IND=cingulum1IND & cing1EndpointCriteria;
    
    periNOTPlane=bsc_makePlanarROI(labelNifti,periCTop+5, 'z');
    periNOTPlane.coords=periNOTPlane.coords(periNOTPlane.coords(:,2)>PosteriorCCcutoffA,:);
    
    [cingulum2, cingulum2IND]=wma_SegmentFascicleFromConnectome(wbfg, [{periCROI} {periCantPlane} {caudatePlaneCut} {superiorCutPlane} {superiorCutPlane2}  {interHemisphere} {periNOTPlane} {ensureArcPlane} {CingCurveLimit}], {'and', 'not','not','and','and', 'not', 'not','not','and'}, 'cingulum2');
    cing2Criteria1=midpoints(:,3)>thaltop & ~or( endpoints1(:,3)> periCTop & endpoints1(:,2)> bstemAnt, endpoints2(:,3)> periCTop & endpoints2(:,2)> bstemAnt );
    cing2Criteria2=midpoints(:,3)>thaltop & midpoints(:,2)> PosteriorCCcutoffA;
    cingulum2IND=cingulum2IND & cing2Criteria1 & cing2Criteria2 & ~cingulum1IND;

    
    [cingulum3, cingulum3IND]=wma_SegmentFascicleFromConnectome(wbfg, [  {cingPostROI} {superiorCutPlane} {superiorCutPlane2} {interHemisphere} {posteriorPlane} {caudateCutLatUpper}], { 'and', 'not','not', 'not', 'not', 'not'}, 'cingulum3');
    
    fprintf('\n cingulum segmented');
    
    %LeftCingulum1IND=LeftStreams & cingulum1IND;
    %LeftCingulum2IND=LeftStreams & cingulum2IND;
    %LeftCingulum3IND=LeftStreams & cingulum3IND;
    hugePutamen=bsc_roiFromFSnums(fsDir,lentiLut(leftright,:),1,21);
    %[sFOF, sFOFIND]=wma_SegmentFascicleFromConnectome(wbfg, [{anteriorPlane} {stemnot} {thalTotal} {superiorCutPlane} {superiorCutPlane2} {interHemisphere}  {posteriorPlane} {occipitalROI} {antArcPlane} {periCROI} {transverseArcCutPlane}], {'and', 'not', 'not', 'and','and','not', 'and','and', 'not', 'not', 'not'}, 'sFOF');
    [sFOF, sFOFIND]=wma_SegmentFascicleFromConnectome(wbfg, [{anteriorPlane} {stemnot} {thalTotal} {superiorCutPlane} {superiorCutPlane2} {interHemisphere}  {posteriorPlane} {occipitalROI}  {antArcPlane} {ensureArcPlane} {hugePutamen} {transverseArcCutPlane}], {'and', 'not', 'not', 'and','and','not', 'and','and','not','not','not','not'}, 'sFOF');
  
    
    SFOFEndpointCriteria=or(endpoints1(:,2) < antArcCut & endpoints2(:,2) > caudateCutAnt, endpoints2(:,2) < antArcCut & endpoints1(:,2) > caudateCutAnt) & midpoints(:,3) >periCTop;
    sFOFIND=sFOFIND & SFOFEndpointCriteria & ~cingulum1IND & ~cingulum2IND;
    

    
  
  
    
   
    
    tempPoleROI=bsc_roiFromFSnums(fsDir,[144]+sidenum,0);
    planPolareROI=bsc_roiFromFSnums(fsDir,[135]+sidenum,1,5);
    insROI=bsc_roiFromFSnums(fsDir,[118]+sidenum,1,5);
    intROI=bsc_intersectROIs(planPolareROI,insROI);
    intROImin=min(intROI.coords(:,3));
    IntPlane=bsc_planeFromROI(intROI,'medial',fsDir);
    intROIlat=unique(IntPlane.coords(:,1));
    intPlanPost=min(intROI.coords(intROI.coords(:,3)==intROImin,2));
    uncPlane=     bsc_makePlanarROI(labelNifti,intROImin , 'z');
    uncPlane.coords=uncPlane.coords(uncPlane.coords(:,2)<intPlanPost,:);
    uncPlane.coords=uncPlane.coords(abs(uncPlane.coords(:,1))>abs(intROIlat),:);
    %arcPlaneCoronal=     bsc_makePlanarROI(labelNifti,intPlanPost , 'y');
    %arcPlaneCoronal.coords=arcPlaneCoronal.coords(arcPlaneCoronal.coords(:,3)<intROImin,:);
    
    
    
    tempPoleMax=max(tempPoleROI.coords(:,2));
    tempPolePlane=     bsc_makePlanarROI(labelNifti,tempPoleMax , 'y');
    amigRoi=bsc_roiFromFSnums(fsDir,amigLut(leftright,:),0);
    amigMin=min(amigRoi.coords(:,2));
    amigMinPlane=     bsc_makePlanarROI(labelNifti,amigMin , 'y');
    [Unc, UncIND]=wma_SegmentFascicleFromConnectome(wbfg, [ {anteriorPlane} {tempPolePlane} {transvPlane} {amigMinPlane} {uncPlane} {interHemisphere}], { 'and' ,'and' , 'not', 'not', 'and', 'not' }, 'Unc');
    UncIND=UncIND & midpoints(:,3)<thaltop;
    %Unc.fibers=wbfg.fibers(UncIND & midpoints(:,3)<thaltop);
    %LeftUncIND=UncIND & LeftStreams;
    
    insNOT=bsc_roiFromFSnums(fsDir,[  118]+sidenum,0);
    
    [arc, arcIND]=wma_SegmentFascicleFromConnectome(wbfg, [{uncPlane} {superiorCutPlane} {superiorCutPlane2} {amigRoi} {insNOT},{interHemisphere} {thalTotal}], {'and', 'and', 'and', 'not', 'not','not', 'not'}, 'arc');
    %arcEndpointCriteria= abs(endpoints1(:,1)) > abs(caudateCutLat) & abs(endpoints2(:,1)) < abs(caudateCutLat);
    %arc.fibers=wbfg.fibers(arcIND & arcEndpointCriteria & midpoints(:,2)<PosteriorCCcutoffA);
    %arcIND=arcIND & arcEndpointCriteria & midpoints(:,2)<PosteriorCCcutoffA;
   
    
    fprintf('\n Arcuate segmented');
    
    %LeftArcIND=arcIND & LeftStreams;
    
    motorAreaROI= bsc_roiFromFSnums(fsDir,[  128 146 129 168 104 126]+sidenum,1,11);
    
    [corticoSpine, corticoSpineInd]=  bsc_tractByEndpointROIs(wbfg, [{motorAreaROI} {stemroi}]);
    
    
    
    %corticoSpineInd=corticoSpineInd & midpoints(:,1) <0 ;
    
    %LeftCorticoSpineInd=corticoSpineInd' & LeftStreams;
    
    
    thalROI=bsc_roiFromFSnums(fsDir,thalamusLut(leftright),1,3);
    lentiROI=bsc_roiFromFSnums(fsDir,lentiLut(leftright,:),1,3);
    
    thalPostCoord=min(thalROI.coords(:,2));
    thalAntCoord=max(thalROI.coords(:,2));
    thalInfCoord=min(thalROI.coords(:,3));
    thalPostPlane=     bsc_makePlanarROI(labelNifti,thalPostCoord, 'y');
    thalAntPlane=     bsc_makePlanarROI(labelNifti,thalAntCoord, 'y');
    thalInfPlane= bsc_planeFromROI(thalROI,'inferior',fsDir);
    ccTop=max(BodyCCROI.coords(:,3));
    ccTopPlane=     bsc_makePlanarROI(labelNifti,ccTop, 'z');
    
    
    caudInflate=bsc_roiFromFSnums(fsDir,caudateLut(leftright),1,11);
    lentiInflate=bsc_roiFromFSnums(fsDir,lentiLut(leftright,:),1,11);
    thalAntWM=bsc_intersectROIs(caudInflate,lentiInflate);
    thalAntWM.coords=thalAntWM.coords(thalAntWM.coords(:,2)>thalAntCoord,:);
    
    
    [antThalRad, antThalRadIND]=wma_SegmentFascicleFromConnectome(wbfg, [{thalROI} {interHemisphere} {thalPostPlane} {lentiROI} {thalAntWM}], {'endpoints', 'not', 'not', 'not', 'and'}, 'antThal');
    antThalCriteria=endpoints2(:,2) < midpoints(:,2) & endpoints1(:,2) < midpoints(:,2) ;
    antThalRadIND=antThalRadIND & ~antThalCriteria;
    
    %LeftAntThalRadIND=antThalRadIND & LeftStreams;
    
    
    
    
    
    [supThalRad, supThalRadIND]=wma_SegmentFascicleFromConnectome(wbfg, [{thalROI} {interHemisphere} {anteriorPlane}  {lentiROI} {ccTopPlane} {occipitalROI} {thalInfPlane}], {'endpoints', 'not', 'not', 'not', 'and', 'not','not'}, 'supThal');
    supThalCriteria=~(endpoints1(:,3) < midpoints(:,3) & endpoints2(:,3) < midpoints(:,3));
    
    periCInflateROI=bsc_roiFromFSnums(fsDir,[167]+sidenum,1,15);
    [PeriCStreams, PeriCStreamsIND]=wma_SegmentFascicleFromConnectome(wbfg, [{periCInflateROI}], {'endpoints'}, 'periCExclude');
    
    supThalRadIND=supThalRadIND & supThalCriteria & ~PeriCStreamsIND;
    

    
    
    thalPostHalfROI=thalROI;
    thalPostHalfROIUniqueCoords=unique(thalPostHalfROI.coords(:,2));
    thalPostHalfROI.coords=thalPostHalfROI.coords(thalPostHalfROI.coords(:,2)<mean(thalPostHalfROIUniqueCoords),:);
    
    
    
    temporalROI1=bsc_roiFromFSnums(fsDir,[  151]+sidenum,1,15);
    temporalROI2=bsc_roiFromFSnums(fsDir,[  149]+sidenum,1,15);
    temporalWMROI=bsc_intersectROIs(temporalROI1,temporalROI2);
    

    
    wmNotROI=bsc_roiFromFSnums(fsDir,wmLut(leftright),0);
    vesselNotROI=bsc_roiFromFSnums(fsDir,vesselLut(leftright),1,7);
    vesselWMROI=bsc_intersectROIs(wmNotROI,vesselNotROI);
    
    test7= bsc_roiFromFSnums(fsDir,[  174]+sidenum,1,3);
    
    [infThalRad, infThalRadIND]=wma_SegmentFascicleFromConnectome(wbfg, [{thalPostHalfROI} {interHemisphere} {stemroi} {temporalWMROI} {thalPostPlane} {thalAntWM} {vesselWMROI} {test7}], {'endpoints', 'not',  'not', 'and', 'not', 'not','not','not'}, 'infThal');
    
    %LeftInfThalRadIND= LeftStreams & infThalRadIND;
    
   
    postCriteria1= or(endpoints1(:,2) > thalPostCoord &  endpoints2(:,2) < thalPostCoord,endpoints1(:,2) < thalPostCoord &  endpoints2(:,2) > thalPostCoord);
    postCriteria2= ~or(endpoints1(:,2) > thalPostCoord &  endpoints1(:,3) < thalInfCoord,endpoints2(:,2) > thalPostCoord &  endpoints2(:,3) < thalInfCoord);
    
    thalPostPlane.coords=thalPostPlane.coords(thalPostPlane.coords(:,3)<thaltop,:);
    
    [postThalRad, postThalRadIND]=wma_SegmentFascicleFromConnectome(wbfg, [{thalROI} {interHemisphere}  {anteriorPlane}  {thalPostPlane} {stemroi} {thalAntPlane} {occipitalROI}], {'endpoints', 'not',  'not', 'and', 'not', 'not', 'endpoints'}, 'postThal');
    
    [excludeThalRad, excludeThalRadIND]=wma_SegmentFascicleFromConnectome(wbfg, [{thalROI} ], {'both_endpoints'}, 'notThal');
    %postThalRad.fibers=wbfg.fibers(  ~excludeThalRadIND & postThalRadIND & postCriteria1 & postCriteria2);
    
    postThalRadIND= ~excludeThalRadIND & postThalRadIND & postCriteria1 & postCriteria2;
    
    %LeftPostThalRadIND=postThalRadIND & LeftStreams;
    
    
    
    thalAntNot=bsc_planeFromROI(thalTotal,'anterior',fsDir);
    stemPostNot=bsc_planeFromROI(stemroi,'posterior',fsDir);
    
    [~, thalamicoSpinalIndPre]=wma_SegmentFascicleFromConnectome(wbfg, [{transvPlane}, {thalAntNot} {stemPostNot} ], {'not','not','not'}, 'postThal');
    
    [thalamicoSpinal, thalamicoSpinalInd]=bsc_tractByEndpointROIs(wbfg,[{stemroi},{thalROI} ]);
    thalamicoSpinalInd= thalamicoSpinalIndPre' & thalamicoSpinalInd;
    
    %LeftThalamicoSpinalIND=thalamicoSpinalInd & LeftStreams;
    fprintf('\n thalamic tracts segmented');
    
    %%
    %occipitalROISplit=occipitalROI;
    %occipitalROISplit.coords=occipitalROISplit.coords(or(occipitalROISplit.coords(:,3)>15,occipitalROISplit.coords(:,3)<-5),:);
    
    OccAntLimit=bsc_planeFromROI(occipitalROI,'anterior',fsDir);
    
    
    
    testOcc1inf=bsc_roiFromFSnums(fsDir,[ 120]+sidenum,1,5);
    wmROI=bsc_roiFromFSnums(fsDir,[ 119]+sidenum,1,5);
    parOccp=bsc_roiFromFSnums(fsDir,[ 111 166]+sidenum,0);
    
    
    
    occLimit=bsc_makePlanarROI(labelNifti,occLimitVal, 'y');
    
    [VOF1, VOFIND1]=wma_SegmentFascicleFromConnectome(wbfg, [{testOcc1inf} {wmROI} {occLimit}  {interHemisphere} {parOccp}], {'endpoints', 'endpoints', 'not','not', 'not'}, 'vof');
    
    
    
    ventricleInflate=bsc_roiFromFSnums(fsDir,ventricleLut(leftright),1,11);
    ventricleLimit=bsc_planeFromROI([160]+sidenum, 'anterior', fsDir);
    ventricleInflate.coords=ventricleInflate.coords(ventricleInflate.coords(:,2)<unique(ventricleLimit.coords(:,2)),:);
    
    
    lateralOccip1=bsc_roiFromFSnums(fsDir,[145 ]+sidenum,1,7);
    medialOccip1=bsc_roiFromFSnums(fsDir,[162 ]+sidenum,1,7);
    occWM1=bsc_intersectROIs(lateralOccip1,medialOccip1);
    
    
    occAntPlane=bsc_planeFromROI(160+sidenum,'anterior',fsDir);
    
    % top of calcerine  bottom of IPS t
    calcTopPlane=bsc_planeFromROI(145+sidenum,'superior',fsDir);
    
    occAntPlane.coords=occAntPlane.coords(occAntPlane.coords(:,3)>unique(calcTopPlane.coords(:,3)),:);
    
    IPSbottomPlane=bsc_planeFromROI(157+sidenum,'inferior',fsDir);
    lingualPlane=bsc_planeFromROI(122+sidenum,'superior',fsDir);
    
    
    testCriteria= endpoints1(:,3)>unique(calcTopPlane.coords(:,3)) & endpoints2(:,3)>unique(calcTopPlane.coords(:,3)) & midpoints(:,3)>unique(calcTopPlane.coords(:,3)) & midpoints(:,3)<unique(IPSbottomPlane.coords(:,3));
    
    
    
    
    VOFIND1=VOFIND1 & testCriteria ;
    transVSulc=bsc_planeFromROI(152+sidenum,'superior',fsDir);
    fullOccipitalPlane=bsc_planeFromROI(160+sidenum,'anterior',fsDir);
    bottomOccipitalPlane=bsc_planeFromROI(160+sidenum,'anterior',fsDir);
    bottomOccipitalPlane.coords=bottomOccipitalPlane.coords(bottomOccipitalPlane.coords(:,3)>unique(transVSulc.coords(:,3)),:);
    
    subtractROI=bsc_roiFromFSnums(fsDir,143+sidenum,1,7);
    
    
    
    %OccAntLimit
    %occWM
    %[VOF2, VOFIND2]=wma_SegmentFascicleFromConnectome(wbfg, [{occipitalROI} {interHemisphere} {OccAntLimit} {stemroi} {testplaneb}], {'both_endpoints', 'not', 'not', 'not','and'}, 'vof');
    %[VOF2, VOFIND2]=wma_SegmentFascicleFromConnectome(wbfg, [{VOFoccipSup} {VOFoccipInf} {interHemisphere} {OccAntLimit} {stemroi}], {'endpoints','endpoints', 'not', 'not', 'not'}, 'vof');
    [VOF2, VOFIND2]=wma_SegmentFascicleFromConnectome(wbfg, [{occipitalROI} {interHemisphere} {stemroi} {occAntPlane} {occWM1} {fullOccipitalPlane} {subtractROI}], {'both_endpoints', 'not', 'not','not', 'not', 'not','not'}, 'vof');
    
    [VOFAnt, VOFAntIND]=wma_SegmentFascicleFromConnectome(wbfg, [{occipitalROI} {interHemisphere} {stemroi} {occAntPlane} {occWM1} {bottomOccipitalPlane} {subtractROI}], {'both_endpoints', 'not', 'not','not', 'not', 'not','not'}, 'vof');
    
    VOFAntIND=VOFAntIND & ~VOFIND2;
    
    OccipitalMedial=bsc_planeFromROI(152+sidenum,'medial',fsDir);
    OccipitalSupLimit=bsc_planeFromROI(109+sidenum,'inferior',fsDir);
    OccIncPlane=OccipitalSupLimit;
    OccIncPlane.coords=OccIncPlane.coords(abs(OccIncPlane.coords(:,1))<abs(unique(OccipitalMedial.coords(:,1))),:);
    
     [VOF3, VOFIND3]=wma_SegmentFascicleFromConnectome(wbfg, [{occipitalROI} {interHemisphere} {stemroi} {occAntPlane} {occWM1} {fullOccipitalPlane} {subtractROI} {OccIncPlane}], {'both_endpoints', 'not', 'not','not', 'not', 'not','not','and'}, 'vof');
    
 
    
    testCriteria2= or(endpoints1(:,3)<unique(lingualPlane.coords(:,3)) & endpoints2(:,3)>unique(IPSbottomPlane.coords(:,3)) , endpoints2(:,3)<unique(lingualPlane.coords(:,3)) & endpoints1(:,3)>unique(IPSbottomPlane.coords(:,3)));
    testCriteria3= endpoints1(:,2)<unique(OccAntLimit.coords(:,2)) & endpoints2(:,2)<unique(OccAntLimit.coords(:,2));
    testCriteria4= midpoints(:,2)<unique(bottomOccipitalPlane.coords(:,2));
    %VOF2.fibers=wbfg.fibers(VOFIND2 & testCriteria2 & testCriteria3 & testCriteria4 & ~VOFIND1);
    
        VOFIND3=VOFIND3 & testCriteria2 & testCriteria3 & testCriteria4 & ~VOFIND1 ;
    VOFIND2=VOFIND2 & testCriteria2 & testCriteria3 & testCriteria4 & ~VOFIND1 &~VOFIND3 ;
    VOFAntIND= VOFAntIND & testCriteria2 & testCriteria3 & testCriteria4 & ~VOFIND2 &~VOFIND3;
 
    
    fprintf('\n Occipital tract segmentation complete')
    

    if leftright==1
        fprintf('\n Left segmentation complete')
        classification.index(FOFIND & LeftStreams)=find( strcmp(classification.names,'Left IFOF'));
        classification.index(cingulum1IND & LeftStreams)=find( strcmp(classification.names,'Left Cingulum1'));
        classification.index(cingulum2IND & LeftStreams)=find( strcmp(classification.names,'Left Cingulum2'));
        classification.index(cingulum3IND & LeftStreams)=find( strcmp(classification.names,'Left Cingulum3'));
        classification.index(sFOFIND & LeftStreams)=find( strcmp(classification.names,'Left SFOF'));
        classification.index(UncIND & LeftStreams)=find( strcmp(classification.names,'Left Uncinate'));
        classification.index(arcIND & LeftStreams)=find( strcmp(classification.names,'Left Arcuate'));
        classification.index(corticoSpineInd' & LeftStreams)=find( strcmp(classification.names,'Left cortico-spinal'));
        classification.index(antThalRadIND & LeftStreams)=find( strcmp(classification.names,'Left anterior thalamic'));
        classification.index(infThalRadIND & LeftStreams)=find( strcmp(classification.names,'Left inferior thalamic'));
        classification.index(supThalRadIND & LeftStreams)=find( strcmp(classification.names,'Left superior thalamic'));
        classification.index(postThalRadIND & LeftStreams)=find( strcmp(classification.names,'Left posterior thalamic'));
        classification.index(thalamicoSpinalInd' & LeftStreams)=find( strcmp(classification.names,'Left Thalamico-spinal'));
        classification.index(VOFIND1 & LeftStreams)=find( strcmp(classification.names,'Left OccipitalU'));
        classification.index(VOFIND2 & LeftStreams)=find( strcmp(classification.names,'Left VOF'));
        classification.index(VOFIND3 & LeftStreams)=find( strcmp(classification.names,'Left VisionU'));
        classification.index(VOFAntIND & LeftStreams)=find( strcmp(classification.names,'Left AntVOF'));
        
        
        
    else
              fprintf('\n Right segmentation complete')
        classification.index(FOFIND & ~LeftStreams)=find( strcmp(classification.names,'Right IFOF'));
        classification.index(cingulum1IND & ~LeftStreams)=find( strcmp(classification.names,'Right Cingulum1'));
        classification.index(cingulum2IND & ~LeftStreams)=find( strcmp(classification.names,'Right Cingulum2'));
        classification.index(cingulum3IND & ~LeftStreams)=find( strcmp(classification.names,'Right Cingulum3'));
        classification.index(sFOFIND & ~LeftStreams)=find( strcmp(classification.names,'Right SFOF'));
        classification.index(UncIND & ~LeftStreams)=find( strcmp(classification.names,'Right Uncinate'));
        classification.index(arcIND & ~LeftStreams)=find( strcmp(classification.names,'Right Arcuate'));
        classification.index(corticoSpineInd' & ~LeftStreams)=find( strcmp(classification.names,'Right cortico-spinal'));
        classification.index(antThalRadIND & ~LeftStreams)=find( strcmp(classification.names,'Right anterior thalamic'));
        classification.index(infThalRadIND & ~LeftStreams)=find( strcmp(classification.names,'Right inferior thalamic'));
        classification.index(supThalRadIND & ~LeftStreams)=find( strcmp(classification.names,'Right superior thalamic'));
        classification.index(postThalRadIND & ~LeftStreams)=find( strcmp(classification.names,'Right posterior thalamic'));
        classification.index(thalamicoSpinalInd' & ~LeftStreams)=find( strcmp(classification.names,'Right Thalamico-spinal'));
        classification.index(VOFIND1 & ~LeftStreams)=find( strcmp(classification.names,'Right OccipitalU'));
        classification.index(VOFIND2 & ~LeftStreams)=find( strcmp(classification.names,'Right VOF'));
        classification.index(VOFIND3 & ~LeftStreams)=find( strcmp(classification.names,'Right VisionU'));
        classification.index(VOFAntIND & ~LeftStreams)=find( strcmp(classification.names,'Right AntVOF'));
        
    end
    
end
end
