function [niisave] = csc_normalizedensityNifti(nii,saveFileName)

% Read nifti file of fiber density data and then normalize.
% The fiber density value of each voxel is divided by maximum value.

% (C) Hiromasa Takemura, Stanford VISTA team
%Edit by Dan Bullock in 2017

%% Read nifti
if ischar(nii)
    nii = niftiRead(nii);
else
end


niisave = nii;

%% Find out max, and then normalize
%% is this really normalizing?  or is this just flattening?

maxnum = max(nii.data);
maxmax = max(maxnum);
maxmaxmax = max(maxmax);
niisave.data = nii.data/maxmaxmax;

%% Save file
niftiWrite(niisave, saveFileName);