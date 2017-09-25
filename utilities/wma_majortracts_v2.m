function [classification] = wma_majortracts_v2(dt, wbfg)
%
% Segments 20 tracts defined in the JHU White Matter Atlas (Zhang et al., 2008).
%
% THis is modified version of AFQ_SegmentFiberGroups.m adapted to be
% compatible to the ENCODE framework and LiFE.
%
%  [fg_classified, fg_unclassified, classification, fg] = ...
%      wma_majortracts(dt6_file, fg);
%
%  Fibers are segmented in two steps. Fibers become candidates for a fiber
%  group if the pass through the 2 waypoint ROIs that define the
%  tracjectory of the tract. Then each fiber is compared to a fiber
%  proability map and high probability fibers are retained in the group.
%  The segmentation alogrithm is based on:
%
%  Hua K, Zhang J, Wakana S, Jiang H, Li X, Reich DS, Calabresi PA, Pekar
%  JJ, van Zijl PC, Mori S. 2008. Tract probability maps in stereotaxic
%  spaces: analyses of white matter anatomy and tract-specific
%  quantification. Neuroimage 39(1):336-47.
%
%  Zhang, W., Olivi, A., Hertig, S. J., van Zijl, P., & Mori, S. (2008).
%  Automated fiber tracking of human brain white matter using diffusion
%  tensor imaging. Neuroimage, 42(2), 771-777.
%
% Input parameters:
% dt6_file                 - A path to a dt6 file or a dt6 object
%
% wbfg                     - A file with previously tracked elsewhere fibers
%                          to be categorized.

% Output parameters:
% fg_ classified  - fibers structure containing all fibers assigned to
%                   one of Mori Groups. Respective group labeles are stored
%                   in fg.subgroups field.
% fg_unclassified - fiber structure containing the rest of the (not Mori) fibers.
% classification  - This variable gives the fiber group names and the group
%                   fiber group number for each fiber in the input group
%                   fg.  classification is a structure with two fields.
%                   classification.names is a cell array where each cell is
%                   the name of that fiber group. For example
%                   classification.names{3} = 'Corticospinal tract L'.
%                   classification.index is a vector that defines which
%                   group number each fiber in the origional fiber group
%                   was assigned to. For example
%                   classification.index(150)=3 means that fg.fibers(150)
%                   is part of the corticospinal tract fiber group.  The
%                   values in classification may not match the origional
%                   fiber group because of pre-processing.  However they
%                   will match the output fg which is the origional group
%                   with preprocessing.
% fg              - This is the origional pre-segmented fiber group.  It
%                   may differ slightly from the input due to preprocessing
%                   (eg splitting fibers that cross at the pons, removing
%                   fibers that are too short)
%
% Example:
%    data = '/home/jyeatman/matlab/svn/vistadata/AFQ';
%    dt6_file = fullfile(data, 'dt6.mat');
%    fg      = wma_majortracts(dt6_file);
%
% See also: dtiSplitInterhemisphericFibers
%
% (c) Vistalab
%
% Modified by Franco Pestilli Indiana University 2017

%% Initialize SPM defualt parameters for normalization
spm_get_defaults; global defaults;
defaults.normalise.estimate.smosrc  = 8;
defaults.normalise.estimate.smoref  = 0;
defaults.normalise.estimate.regtype = 'mni';
defaults.normalise.estimate.weight  = '';
defaults.normalise.estimate.cutoff  = 25;
defaults.normalise.estimate.nits    = 16;
defaults.normalise.estimate.reg     = 1;
params = defaults.normalise.estimate;

%% Check arguments

minDist = 2; % defualt is .89;
%is it though?  This was set to 2 by default.

    Atlas = 'MNI_JHU_tracts_prob.nii.gz';
    
% Set the directory where templates can be found
tdir = fullfile(fileparts(which('mrDiffusion.m')), 'templates');

%% Read the data
% Load the dt6 file from disk
if ischar(dt)
    dt = dtiLoadDt6(dt);
end

%% Spatially normalize diffusion data with the MNI (ICBM) template
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% If you wanted to inverse-normalize the maps to this subject's brain:
% invDef.outMat = moriTracts.qto_ijk;
% bb = mrAnatXformCoords(dt.xformToAcpc,[1 1 1; size(dt.b0)]);
% tprob = mrAnatResliceSpm(tprob, invDef, bb, dt.mmPerVoxel, [1 1 1 0 0 0]);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
template = fullfile(tdir,'MNI_JHU_T2.nii.gz');

% Rescale image valueds to get better gray/white/CSF contrast
alignIm = mrAnatHistogramClip(double(dt.b0),0.3,0.99);

% Compute normalization
[sn, Vtemplate, invDef] = mrAnatComputeSpmSpatialNorm(alignIm, dt.xformToAcpc, template, params);

% check the normalization
% mm = diag(chol(Vtemplate.mat(1:3,1:3)'*Vtemplate.mat(1:3,1:3)))';
bb = mrAnatXformCoords(Vtemplate.mat,[1 1 1; Vtemplate.dim]);
alignIm_sn = mrAnatResliceSpm(alignIm, sn, bb, [2 2 2], [1 1 1 0 0 0], 0);
tIm = mrAnatResliceSpm(double(Vtemplate.dat), inv(Vtemplate.mat), bb, [2 2 2], [1 1 1 0 0 0], 0);

im(:,:,:,1) = uint8(tIm);
im(:,:,:,2) = uint8(round(clip(alignIm_sn)*255));
im(:,:,:,3) = im(:,:,:,2);


%% Compute fiber group probabilities using the atlas proceedure of Hua.2008
%
% Load the Mori atlas maps these are saved in nifti images
moriTracts = readFileNifti(fullfile(tdir, Atlas));


% 15 is a subregion of 19 and 16 a subregion of 20. To better separate them,
% we subtract 19 from 15 and 20 from 16.
moriTracts.data(:,:,:,15) = moriTracts.data(:,:,:,15)-moriTracts.data(:,:,:,19);
moriTracts.data(:,:,:,16) = moriTracts.data(:,:,:,16)-moriTracts.data(:,:,:,20);

% Load the fiber group labels
% labels = readTab(fullfile(tdir,'MNI_JHU_tracts_prob.txt'),',',false);
% labels = labels(1:20,2);
labels = {'Left_Thalamic_Radiation','Right_Thalamic_Radiation', ...
    'Left_Corticospinal','Right_Corticospinal', ...
    'Left_Cingulum Cingulate', 'Right Cingulum_Cingulate', ...
    'Left_Cingulum Hippocampus','Right Cingulum_Hippocampus', ...
    'Callosum Forceps_Major', 'Callosum_Forceps_Minor', ...
    'Left_IFOF','Right_IFOF', ...
    'Left_ILF','Right_ILF', ...
    'Left_SLF','Right_SLF', ...
    'Left_Uncinate','Right_Uncinate', ...
    'Left_Arcuate','Right_Arcuate'};

% Warp the fibers to the MNI standard space so they can be compared to the
% template
fg_sn = dtiXformFiberCoords(wbfg, invDef);

% moriTracts.data is a an XxYxZx20 array contianing the 20 Mori probability
% atlases (range is 0-100 where 100 represents p(1)).
sz = size(moriTracts.data);

% fg_sn fiber coords are in MNI space- now convert them to atlas space by
% applying the affine xform from the atlas NIFTI header. Since the atlas is
% already in MNI space, this transform will just account for any
% translation and scale differences between the atlas maps and the MNI
% template used to compute our sn.
fgCoords = mrAnatXformCoords(moriTracts.qto_ijk, horzcat(fg_sn.fibers{:}));
fgLen    = cellfun('size',wbfg.fibers,2);

clear fg_sn;   % all we need from fg_sn is now stored in fgCoords

% Now loop over the 20 atlases and get the atlas probability score for each
% fiber point. We collapse the scores across all points in a fiber by taking
% the mean. Below, we will use these 20 mean scores to categorize the fibers.
fp = zeros(sz(4),numel(wbfg.fibers));
fprintf('\n Obtaining probabality atlases')
for(ii=1:sz(4))
    % Get the Mori atlas score for each point in the fibers using
    % trilinear interpolation.
    p = myCinterp3(double(moriTracts.data(:,:,:,ii))/100, sz([1,2]), sz(3), fgCoords(:,[2,1,3]));
    
    % The previous line interpolated one giant array with all fiber points
    % concatenated. The next loop will separate the coordinates back into
    % fibers and take the mean score for the points within each fiber.
    fiberCoord = 1;
    for(jj=1:numel(wbfg.fibers))
        fp(ii,jj)  = nanmean(p([fiberCoord:fiberCoord + fgLen(jj)-1]));
        fiberCoord = fiberCoord+fgLen(jj);
    end
end
clear p fgCoords;

%% Find fibers that pass through both waypoint ROIs eg. Wakana 2007
%
% Warp Mori ROIs to individual space; collect candidates for each fiber
% group based on protocol of 2 or > ROIs a fiber should travel through. The
% following ROIs are saved within
% trunk/mrDiffusion/templates/MNI_JHU_tracts_ROIs folder and are created
% using MNI template as described in Wakana et al.(2007) Neuroimage 36 with
% a single modification: For SLFt Roi2, they recommend drawing the ROI at
% the AC level, whereas we use a slice just inferior of CC splenium. The
% reason for this modification is that Wakana et al. ACPC aligned images
% appear different from MNI images (the latter we use for defininng ROIs).
% If defining SLFt-Roi2 on a slice actually at the AC level (althought
% highly consistently across human raters), many SLFt fibers were not
% correctly labeled as they extend laterally into temporal lobe just above
% the aforementioned ROI plane.
% A 20x2 cell array containing the names of both waypoint ROIs for each
% of 20 fiber groups
moriRois={'ATR_roi1_L.nii.gz',  'ATR_roi2_L.nii.gz'; 'ATR_roi1_R.nii.gz', 'ATR_roi2_R.nii.gz'; ...
    'CST_roi1_L.nii.gz', 'CST_roi2_L.nii.gz'; 'CST_roi1_R.nii.gz',  'CST_roi2_R.nii.gz'; ...
    'CGC_roi1_L.nii.gz', 'CGC_roi2_L.nii.gz'; 'CGC_roi1_R.nii.gz', 'CGC_roi2_R.nii.gz'; ...
    'HCC_roi1_L.nii.gz', 'HCC_roi2_L.nii.gz'; 'HCC_roi1_R.nii.gz', 'HCC_roi2_R.nii.gz';...
    'FP_R.nii.gz', 'FP_L.nii.gz'; ...
    'FA_L.nii.gz', 'FA_R.nii.gz'; ...
    'IFO_roi1_L.nii.gz', 'IFO_roi2_L.nii.gz'; 'IFO_roi1_R.nii.gz', 'IFO_roi2_R.nii.gz'; ...
    'ILF_roi1_L.nii.gz', 'ILF_roi2_L.nii.gz'; 'ILF_roi1_R.nii.gz', 'ILF_roi2_R.nii.gz'; ...
    'SLF_roi1_L.nii.gz', 'SLF_roi2_L.nii.gz'; 'SLF_roi1_R.nii.gz', 'SLF_roi2_R.nii.gz'; ...
    'UNC_roi1_L.nii.gz', 'UNC_roi2_L.nii.gz'; 'UNC_roi1_R.nii.gz', 'UNC_roi2_R.nii.gz'; ...
    'SLF_roi1_L.nii.gz', 'SLFt_roi2_L.nii.gz'; 'SLF_roi1_R.nii.gz', 'SLFt_roi2_R.nii.gz'};

% Make an ROI for the mid saggital plane
midSaggitalRoi = dtiRoiMakePlane([0, dt.bb(1, 2), dt.bb(1, 3); ...
    0 , dt.bb(2, 2) , dt.bb(2, 3)], ...
    'midsaggital', 'g');
keep1 = zeros(length(wbfg.fibers), size(moriRois, 1));
keep2 = zeros(length(wbfg.fibers), size(moriRois, 1));

% Find fibers that cross mid saggital plane
[~, ~, InterHemisphericFibers] = dtiIntersectFibersWithRoi([], 'not', [], midSaggitalRoi, wbfg);
%NOTICE: ~keep3 (not "keep3") will mark fibers that DO NOT cross
%midSaggitalRoi.

keep3  = repmat(InterHemisphericFibers, [1 size(moriRois, 1)]);
fgCopy = wbfg;
fgCopy.subgroup=[];
fprintf('\n beginning segmentation of Mori tracts \n')
for roiID = 1:size(moriRois, 1)
    % Load the nifti image containing ROI-1 in MNI space
    fprintf('\n %s to %s \n', moriRois{roiID, 1}, moriRois{roiID, 2})
    ROI_img_file=fullfile(tdir, 'MNI_JHU_tracts_ROIs',  [moriRois{roiID, 1}]);
    % Transform ROI-1 to an individuals native space
    
    % Default is to use the spm normalization unless a superior
    % ANTS normalization was passed in
    
    [invDef, roi1(roiID)]=wma_CreateRoiFromMniNifti(dt, ROI_img_file, invDef);
    
    
    % Find fibers that intersect the ROI
    [~,~, keep1(:, roiID)] = dtiIntersectFibersWithRoi([], 'and', minDist, roi1(roiID), wbfg);
    keepID1 = find(keep1(:, roiID));
    
    % Load the nifit image containing ROI-2 in MNI space
    ROI_img_file = fullfile(tdir, 'MNI_JHU_tracts_ROIs',  [moriRois{roiID, 2}]);
    
    % Transform ROI-2 to an individuals native space
    
    [invDef, ~] = wma_CreateRoiFromMniNifti(dt, ROI_img_file, invDef);
    % Default is to use the spm normalization unless a superior
    % ANTS normalization was passed in
    
    [invDef, roi2(roiID)] = wma_CreateRoiFromMniNifti(dt, ROI_img_file, invDef);
    
    
    %To speed up the function, we intersect with the second ROI not all the
    %fibers, but only those that passed first ROI.
    fgCopy.fibers = wbfg.fibers(keepID1(keepID1>0));
    [~,~, keep2given1] = dtiIntersectFibersWithRoi([], 'and', minDist, roi2(roiID), fgCopy);
    keep2(keepID1(keep2given1), roiID) = true;
    % end roi loop
end

clear fgOut contentiousFibers keepID
%Note: forceps major and minor should NOT have interhemipsheric fibers
% excluded
keep3(:, 9:10)=keep3(:, 9:10).*0;

% fp is the variable containing each fibers score for matching the
% atlas. We will set each fibers score to 0 if it does not pass through
% the necessary ROIs
fp(~(keep1'&keep2'&~keep3'))=0;

%Also note: Tracts that cross through slf_t rois should be automatically
%classified as slf_t, without considering their probs.
fp(19, (keep1(:, 19)'&keep2(:, 19)'&~keep3(:, 19)'))=max(fp(:));
fp(20, (keep1(:, 20)'&keep2(:, 20)'&~keep3(:, 20)'))=max(fp(:));


%% Eliminate fibers that don't match any of the atlases very well
% We have a set of atlas scores for each each fiber. To categorize the
% fibers, we will find the atlas with the highest score (using 'sort').
[atlasScore,atlasInd] = sort(fp,1,'descend');
% Eliminate fibers that don't match any of the atlases very well:
goodEnough = atlasScore(1,:) ~= 0;
curAtlasFibers = cell(1,sz(4));
for ii=1:sz(4)
    curAtlasFibers{ii} = find(atlasInd(1,:)==ii & goodEnough);
end

% Populate the structure denoting the fiber group number that each fiber in
% the origional wholebrain group was assigned to
fiberIndex = zeros(length(wbfg.fibers),1);

for ii = 1 : length(curAtlasFibers)
    fiberIndex(curAtlasFibers{ii}) = ii;
end

% Create a structure with fiber indices and group names.
classification.index = fiberIndex;
classification.names = labels;

return;