function [fgsegment_lmm_gradient, fgsegment_rmm_gradient, pass1L, pass1R,  L_pArc, R_pArc, L_pArc_fibersIndices, R_pArc_fibersIndices, L_pArc_vot, R_pArc_vot, L_pArc_vot_Indexes, R_pArc_vot_Indexes] =bsc_segmentVOF(fg, fsROIdir, outdir, dtFile, L_arcuateFile, R_arcuateFile)

%% This script identify the VOF from whole-brain tractography using AFQ-VOF toolbox.
%%
%% Prerequisite:
% freesurferROIs in mat format:
% whole-brain connectome file, in fg format
%
% (C) Hiromasa Takemura, CiNet HHS 2016
% mod: DNB 2017, Bloomington

disp('loading the relevant files if necessary...')
if ~and(isstruct (L_arcuateFile), isstruct (R_arcuateFile))
    L_arcuate = fgRead(L_arcuateFile);
    R_arcuate = fgRead(R_arcuateFile);
else
 L_arcuate=L_arcuateFile;  
 R_arcuate=R_arcuateFile;  
end
if ~isstruct (dtFile)
    dt = dtiLoadDt6(dtFile);
else
    dt=dtFile;
end
fprintf('Segmenting the VOF from Connectome ...\n')
% Segment the VOF using AFQ-VOF extension
[L_VOF, R_VOF, L_VOF_Indexes, R_VOF_Indexes,  L_pArc, R_pArc, L_pArc_fibersIndices, R_pArc_fibersIndices, L_pArc_vot, R_pArc_vot, L_pArc_vot_Indexes, R_pArc_vot_Indexes] = AFQ_FindVOF(fg,L_arcuate,R_arcuate,fsROIdir,outdir, [],[],dt);


% Exclude streamlines with outlier gradient
[fgsegment_lmm_gradient, L_Indexes] = vofe_gradient_removeoutlier(L_VOF,2);
[fgsegment_rmm_gradient, R_Indexes] = vofe_gradient_removeoutlier(R_VOF,2);
BoolIndex1L=logical(L_Indexes);
BoolIndex1R=logical(R_Indexes);
pass1L=L_VOF_Indexes(BoolIndex1L);
pass1R=R_VOF_Indexes(BoolIndex1R);

% Exclude outlier streamline in terms of position and length
% IN THE MODED VERSION WE ARE LEAVING THIS OUT BECAUSE NONE OF THE OTHER
% TRACTS HAVE HAD mbaComputeFibersOutliers done on them.
% [L_VOF_clean, keepFascicleL] = mbaComputeFibersOutliers(fgsegment_lmm_gradient,3,3,40);
% [R_VOF_clean, keepFascicleR] = mbaComputeFibersOutliers(fgsegment_rmm_gradient,3,3,40);
% pass2L=pass1L(keepFascicleL);
% pass2R=pass1R(keepFascicleR);





