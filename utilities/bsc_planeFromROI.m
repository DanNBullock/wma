function [roiOUT] = bsc_planeFromROI(roiIN, location,fsDir)

labelNifti= wma_getAsegFile(fsDir , '2009');


if isnumeric(roiIN)
    
    roiIN=bsc_roiFromFSnums(fsDir,[roiIN ],0);
    
elseif isstruct(roiIN)
    
    %no action
else
    keyboard
end

LRflag=mean(roiIN.coords(:,1))<0;

switch lower(location)
    case 'top'
        roiCoord=max(roiIN.coords(:,3));
        roiOUT=bsc_makePlanarROI(labelNifti,roiCoord, 'z');
    case 'superior'
        roiCoord=max(roiIN.coords(:,3));
        roiOUT=bsc_makePlanarROI(labelNifti,roiCoord, 'z');
    case 'bottom'
        roiCoord=min(roiIN.coords(:,3));
        roiOUT=bsc_makePlanarROI(labelNifti,roiCoord, 'z');
    case 'inferior'
        roiCoord=min(roiIN.coords(:,3));
        roiOUT=bsc_makePlanarROI(labelNifti,roiCoord, 'z');
    case 'anterior'
        roiCoord=max(roiIN.coords(:,2));
        roiOUT=bsc_makePlanarROI(labelNifti,roiCoord, 'y');
    case 'posterior'
        roiCoord=min(roiIN.coords(:,2));
        roiOUT=bsc_makePlanarROI(labelNifti,roiCoord, 'y');
    case 'medial'
        if LRflag
            roiCoord=max(roiIN.coords(:,1));
            roiOUT=bsc_makePlanarROI(labelNifti,roiCoord, 'x');
        else
            roiCoord=min(roiIN.coords(:,1));
            roiOUT=bsc_makePlanarROI(labelNifti,roiCoord, 'x');
        end
    case 'lateral'
        if LRflag
            roiCoord=min(roiIN.coords(:,1));
            roiOUT=bsc_makePlanarROI(labelNifti,roiCoord, 'x');
        else
            roiCoord=max(roiIN.coords(:,1));
            roiOUT=bsc_makePlanarROI(labelNifti,roiCoord, 'x');
        end
        
end

end
            