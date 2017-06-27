function [results, classification]= wma_segmentAndAnalyze(fe,dt6,fsDIR)
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
if ischar(fe)
    load(fe);
    %haven't been able to test this with (path to wbFG passed in) yet.
    %if it is a fe structure, get the wbFG out of it
end

[classification]=bsc_AFQseg_and_bloomtracks_v2(fe,dt6,fsDIR);


[tractStats] = wma_multiTractAnalysis(classification,fe,dt6);


[~, results]= bsc_feAndAFQqualityCheck(fe, classification);

results.AFQstats.tractStats=tractStats;
results.AFQstats.classification=classification;

end

