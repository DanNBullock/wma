function [classification] =bsc_cerebellarTracts(wbfg, fsDir)
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
    
    
    sidenum=10000+leftright*1000;
    
    
    
    if leftright==2
        [cerebellar] =bsc_roiFromFSnums(fsDir,[46],0);
        [cerebellarCortex] =bsc_roiFromFSnums(fsDir,[47],0);
        cerebWhole=bsc_roiFromFSnums(fsDir,[47,46],0);
        Lenti=bsc_roiFromFSnums(fsDir,[12,13],0);
        contramotor=bsc_roiFromFSnums(fsDir,[128,103,168]+sidenum-1000,7);
        contrafrontal=bsc_roiFromFSnums(fsDir,[115,101,105,116,124,155,154]+sidenum-1000,7);
        ipsimotor=bsc_roiFromFSnums(fsDir,[128,103,168]+sidenum,7);
        ipsifrontal=bsc_roiFromFSnums(fsDir,[115,101,105,116,124,155,154]+sidenum,7);
        
        tempBottomCut=bsc_roiFromFSnums(fsDir,[123]+sidenum-1000,0);
        medCut=bsc_roiFromFSnums(fsDir,[15],1,5);
        
        mergedOCTROI.name='Rcerebellar';
        redN=bsc_roiFromFSnums(fsDir,[28],0);
        thal=bsc_roiFromFSnums(fsDir,[49],0);
    else
        [cerebellar] =bsc_roiFromFSnums(fsDir,[7],0);
        [cerebellarCortex] =bsc_roiFromFSnums(fsDir,[8],0);
        cerebWhole=bsc_roiFromFSnums(fsDir,[7,8],0);
        
        medCut=bsc_roiFromFSnums(fsDir,[15],1,5);
        
        mergedOCTROI.name='Lcerebellar';
        thal=bsc_roiFromFSnums(fsDir,[10],0);
        redN=bsc_roiFromFSnums(fsDir,[60],0);
        Lenti=bsc_roiFromFSnums(fsDir,[51,52],0);
        contrafrontal=bsc_roiFromFSnums(fsDir,[115,101,105,116,124,155,154]+sidenum+1000,1,7);
        contramotor=bsc_roiFromFSnums(fsDir,[128,103,168]+sidenum+1000,1,7);
        ipsifrontal=bsc_roiFromFSnums(fsDir,[115,101,105,116,124,155,154]+sidenum,1,7);
        ipsimotor=bsc_roiFromFSnums(fsDir,[128,103,168]+sidenum,1,7);
         tempBottomCut=bsc_roiFromFSnums(fsDir,[123]+sidenum+1000,0);
        
    end
    
    fprintf('\n Hemispheric ROIs created')
    
       bstem=bsc_roiFromFSnums(fsDir,[16],0);
       corpus=bsc_roiFromFSnums(fsDir,[255,254,253,252,251],1,5);
       

         tempBottomCutCoord=min(tempBottomCut.coords(:,3));
         bstem.coords=bstem.coords(bstem.coords(:,3)<tempBottomCutCoord,:);
    
    atlasNifti = wma_getAsegFile(fsDir , '2009');
    LentiTopCoord=max(Lenti.coords(:,3));
    [topNotROI] =bsc_makePlanarROI(atlasNifti,LentiTopCoord,'z');
    
  
    
    
    [red2ceb, red2cebIND]=bsc_tractByEndpointROIs(wbfg, [{cerebWhole} {redN}]);
    [thal2cebtex, thal2cebtexIND]=bsc_tractByEndpointROIs(wbfg, [{cerebWhole} {thal}]);
    [Lenti2cebtex, Lenti2cebtexIND]=bsc_tractByEndpointROIs(wbfg, [{cerebWhole} {Lenti}]);
    [contracorticopontine, contracorticopontineIND]=bsc_tractByEndpointROIs(wbfg, [{cerebWhole} {contramotor}]);
    [contrafrontopontine, contrafrontopontineIND]=bsc_tractByEndpointROIs(wbfg, [{cerebWhole} {contrafrontal}]);
    [ipsicorticopontine, ipsicorticopontineIND]=bsc_tractByEndpointROIs(wbfg, [{cerebWhole} {ipsimotor}]);
    [ipsifrontopontine, ipsifrontopontineIND]=bsc_tractByEndpointROIs(wbfg, [{cerebWhole} {ipsifrontal}]);
    

    [~, CereballarFibersIND]=wma_SegmentFascicleFromConnectome(wbfg, [{medCut}], {'and'}, 'ceb');
    [~, lentiOutIND]=wma_SegmentFascicleFromConnectome(wbfg, [{topNotROI}], {'not'}, 'dud');
        [~, bstemOutIND]=wma_SegmentFascicleFromConnectome(wbfg, [{bstem}], {'not'}, 'dud');
                [~, corpusOutIND]=wma_SegmentFascicleFromConnectome(wbfg, [{corpus}], {'not'}, 'dud');
        
    
    Lenti2cebtex.fibers= wbfg.fibers(lentiOutIND & Lenti2cebtexIND' & bstemOutIND );
    
    red2ceb1=red2ceb;
    
    red2ceb2=red2ceb;
    
    
    red2ceb1.fibers=wbfg.fibers(CereballarFibersIND & red2cebIND'  & bstemOutIND);
    red2ceb2.fibers=wbfg.fibers(~CereballarFibersIND & red2cebIND' &  bstemOutIND);
    
    thal2cebtex.fibers=wbfg.fibers(thal2cebtexIND'  & bstemOutIND);
    
    contracorticopontine.fibers=wbfg.fibers(contracorticopontineIND' & corpusOutIND & bstemOutIND);
    contrafrontopontine.fibers=wbfg.fibers(contrafrontopontineIND' & corpusOutIND & bstemOutIND);
    ipsicorticopontine.fibers=wbfg.fibers(ipsicorticopontineIND' & corpusOutIND & bstemOutIND);
    ipsifrontopontine.fibers=wbfg.fibers(ipsifrontopontineIND' & corpusOutIND & bstemOutIND);
    
    fprintf('\n Hemispheric segmentation complete.');
    

    
    if leftright==1
        %named from cerebellar source
        LeftLenti2ceb=Lenti2cebtex;
        LeftLenti2cebIndexes=lentiOutIND & Lenti2cebtexIND' & bstemOutIND;
        
        LeftRed2Ceb1=red2ceb1;
        LeftRed2Ceb1Indexes=CereballarFibersIND & red2cebIND'  & bstemOutIND;
        
        LeftRed2Ceb2=red2ceb2;
        LeftRed2Ceb2Indexes=~CereballarFibersIND & red2cebIND'  & bstemOutIND;
        
        LeftThal2ceb=thal2cebtex;
        LeftThal2cebIndexes=thal2cebtexIND'  & bstemOutIND;
        
        LeftContraCorticoPontine=contracorticopontine;
        LeftContraCorticoPontineIndexes=contracorticopontineIND' & corpusOutIND & bstemOutIND;
        
        LeftContraFrontoPontine=contrafrontopontine;
        LeftContraFrontoPontineIndexes=contrafrontopontineIND' & corpusOutIND & bstemOutIND;
        
        LeftIpsiCorticoPontine=ipsicorticopontine;
        LeftIpsiCorticoPontineIndexes=ipsicorticopontineIND' & corpusOutIND & bstemOutIND;
        
        LeftIpsiFrontoPontine=ipsifrontopontine;
        LeftIpsiFrontoPontineIndexes=ipsifrontopontineIND' & corpusOutIND & bstemOutIND;
        
    else
        
        RightLenti2ceb=Lenti2cebtex;
        RightLenti2cebIndexes=lentiOutIND & Lenti2cebtexIND' & bstemOutIND;
        
        RightRed2Ceb1=red2ceb1;
        RightRed2Ceb1Indexes=CereballarFibersIND & red2cebIND'  & bstemOutIND;
        
        RightRed2Ceb2=red2ceb2;
        RightRed2Ceb2Indexes=~CereballarFibersIND & red2cebIND'  & bstemOutIND;
        
        RightThal2ceb=thal2cebtex;
        RightThal2cebIndexes=thal2cebtexIND'  & bstemOutIND;
        
        RightContraCorticoPontine=contracorticopontine;
        RightContraCorticoPontineIndexes=contracorticopontineIND' & corpusOutIND & bstemOutIND;
        
        RightContraFrontoPontine=contrafrontopontine;
        RightContraFrontoPontineIndexes=contrafrontopontineIND' & corpusOutIND & bstemOutIND;
        
        RightIpsiCorticoPontine=ipsicorticopontine;
        RightIpsiCorticoPontineIndexes=ipsicorticopontineIND' & corpusOutIND & bstemOutIND;
        
        RightIpsiFrontoPontine=ipsifrontopontine;
        RightIpsiFrontoPontineIndexes=ipsifrontopontineIND' & corpusOutIND & bstemOutIND;
    end

end
    
    classification=[];

    classification.index=zeros(length(wbfg.fibers),1);
    classification.names={'Right Lenti2ceb' , 'Left Lenti2ceb', 'Right Red2Ceb1', 'Left Red2Ceb1', 'Right Red2Ceb2', 'Left Red2Ceb2', ...
        'Right Thal2ceb', 'Left Thal2ceb', 'Right ContraCorticoPontine', 'Left ContraCorticoPontine', ...
        'Right ContraFrontoPontine', 'Left ContraFrontoPontine', 'Right IpsiCorticoPontine', 'Left IpsiCorticoPontine',...
        'Right IpsiFrontoPontine','Left IpsiFrontoPontine'};
    classification.index(RightLenti2cebIndexes)=find( strcmp(classification.names,'Right Lenti2ceb'));
    classification.index(LeftLenti2cebIndexes)=find( strcmp(classification.names,'Left Lenti2ceb'));
    
    classification.index(RightRed2Ceb1Indexes)=find( strcmp(classification.names,'Right Red2Ceb1'));
    classification.index(LeftRed2Ceb1Indexes)=find( strcmp(classification.names,'Left Red2Ceb1'));
    
    classification.index(RightRed2Ceb2Indexes)=find( strcmp(classification.names,'Right Red2Ceb2'));
    classification.index(LeftRed2Ceb2Indexes)=find( strcmp(classification.names,'Left Red2Ceb2'));
    
    classification.index(RightThal2cebIndexes)=find( strcmp(classification.names,'Right Thal2ceb'));
    classification.index(LeftThal2cebIndexes)=find( strcmp(classification.names,'Left Thal2ceb'));
    
    classification.index(RightContraCorticoPontineIndexes)=find( strcmp(classification.names,'Right ContraCorticoPontine'));
    classification.index(LeftContraCorticoPontineIndexes)=find( strcmp(classification.names,'Left ContraCorticoPontine'));
    
    classification.index(RightContraFrontoPontineIndexes)=find( strcmp(classification.names,'Right ContraFrontoPontine'));
    classification.index(LeftContraFrontoPontineIndexes)=find( strcmp(classification.names,'Left ContraFrontoPontine'));
    
    classification.index(RightIpsiCorticoPontineIndexes)=find( strcmp(classification.names,'Right IpsiCorticoPontine'));
    classification.index(LeftIpsiCorticoPontineIndexes)=find( strcmp(classification.names,'Left IpsiCorticoPontine'));
    
    classification.index(RightIpsiFrontoPontineIndexes)=find( strcmp(classification.names,'Right IpsiFrontoPontine'));
    classification.index(LeftIpsiFrontoPontineIndexes)=find( strcmp(classification.names,'Left IpsiFrontoPontine'));
    

end
