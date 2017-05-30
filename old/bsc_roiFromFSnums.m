function mergedROI = wma_roiFromFSnums(fsDir,fsROInums, smoothFlag, smoothKernel)
%
% Build a VISTASOFT ROI from a FreeSurfer segmentations and a series of
% indeces to regions in the FreeSurfer segmentaion
%
%  INPUTS:
%  -fsDir:  path to THIS SUBJECT'S freesurfer rectory.
%
%  -fsROInums:  a list of region ID numbers from freesurfer that you would
%  like merged into a single ROI.  i.e. [ 3, 4 , ... 41,42]
%
%  OUTPUTS:
%  -mergedROI: a merged ROI which has coordinates corresponding to each of
%  the fs regions associated with a number entered in the fsROInums.  For
%  example [3,42] would create a grey matter mask roi, while [2, 41] would
%  create a white matter mask roi
%
%  NOTE:  Makes a call to mri_convert, so requires that FreeSurfer be
%  installed and set up properly.
%
% (C) Daniel Bullock 2017 Indiana University


% set smooth flag, default = off
if notDefined('smoothFlag'),   smoothFlag = false;end
if notDefined('smoothKernel'), smoothFlag = false;end
if length( num2str(fsROInums(1))) < 5
    niiPath=fullfile(fsDir,'/mri/','aparc+aseg.nii.gz');
    mgzfile='aparc+aseg.mgz';
else
    niiPath=fullfile(fsDir,'/mri/','aparc.a2009s+aseg.nii.gz');
    mgzfile='aparc.a2009s+aseg.mgz';
end

% Make the aparc aseg file if it does not exist
if ~exist(niiPath,'file')
    spaceChar={' '};
    [status, ~] = system(strcat('mri_convert',spaceChar, mgzfile ,spaceChar, niiPath));
    if status==0
        error('/n Error generating aseg nifti file.  There may be a problem finding the file.')
        keyboard
    end
end

% read in the appropriate aseg niftifile
atlasNifti = niftiRead(niiPath);

% get size of atlasNifti.data and make a blank matrix mask for it
atlasDataSize = size(atlasNifti.data);
blankLabelNifti(1:atlasDataSize(1),1:atlasDataSize(2),1:atlasDataSize(3))=false;

ROImask=blankLabelNifti;
roiNameString=[];
for iRois = fsROInums
    ROImask=or(ROImask,atlasNifti.data==iRois);
    if iRois==fsROInums(end)
        roiNameString=strcat(roiNameString,num2str(iRois));
    else
        roiNameString=strcat(roiNameString,num2str(iRois),'_');
    end
end

% converts nifti formatting to roi formatting
%smooth if you want to
if smoothFlag
   ROImask =~ ((smooth3(ROImask,'box', smoothKernel)) == 0 );
end

%get index coordinates from nifti mask volume
Roi = find(ROImask);
[x1,y1,z1] = ind2sub(size(atlasNifti.data), Roi);

%create new ROI structure
mergedROI = dtiNewRoi(['fs_', roiNameString], 'r');
%apply the qto_xyz transform to transform the mask coordinates to acpc
%cordinates
mergedROI.coords = mrAnatXformCoords(atlasNifti.qto_xyz, [x1,y1,z1]);
end
