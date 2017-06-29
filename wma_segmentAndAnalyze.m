function [results, classificationRAW]= wma_segmentAndAnalyze(fe,dt6,fsDIR)
%
% [results, classification]= wma_segmentAndAnalyze(fe,dt6,fsDIR)
%
% Overview:  This code takes in an fe structure (or FG) and then performs
% segmentation and analysis on both (the fe structure and the segmentation)
%
% Inputs:
%
% -fe: either a LiFE object or a path to one.
%
% -dt6: either a path to a dt6 file or a dt6 object
%
% -fsDIR: path  to THIS SUBJECT'S freesurfer directory
%
% Outputs:
%
% results: a summary structure containing information about both the
% connectome and the segmented tracts.  
%%  preliminary loading
%
tic
if ischar(fe)
    load(fe);
    %haven't been able to test this with (path to wbFG passed in) yet.
    %if it is a fe structure, get the wbFG out of it
end

[classificationRAW]=bsc_AFQseg_and_bloomtracks_v2(fe,dt6,fsDIR);

% only pruning the origional afq tracts
classificationCut= removeOutliersClassification(classificationRAW,fe, 4, 4, 1:20);

[tractStats] = wma_multiTractAnalysis(classificationCut,fe,dt6);

[~, results]= bsc_feAndAFQqualityCheck(fe, classificationCut);

results.AFQstats.tractStats=tractStats;
results.AFQstats.classification=classificationCut;

segTime=toc;
fprintf('\n Segmentation and analysis for %s has taken %4.2f hours.', fe.name, segTime/(60*60))

end

