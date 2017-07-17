function bsc_makeVolumeClassifiedStreams(wbFG, classification, saveDir,subSelect)
%
%   bsc_plotClassifiedStreams(wbFG, classification, t1, view, saveDir,subSelect,colors)
%
%   PURPOSE: This function creates volumetric reperesentations for multiple
%   fiber tracts.  It will do this for all classified fiber tracts or those
%   that are subselected. Also, if you pass in an fe structure it will do
%   this for only validated streamlines.  If you'd like to prune your
%   fibers as well, feel free to use removeOutliersClassification.
%
%  -wbFG:  a structure containing the streamlines referenced in the
%    classification structure.  Can be either a fe structure or a whole
%    brain fiber group.  Will load paths.
%
%  -classification: Either the path to structure or the structure itself.
%   The strucure has a field "names" with (N) names of the tracts classified
%   while the field "indexes" has a j long vector (where  j = the nubmer of
%   streamlines in wbFG (i.e. length(wbFG.fibers)).  This j long vector has
%   a 0 for to indicate a streamline has gone unclassified, or a number 1:N
%   indicatate that the streamline has been classified as a member of tract
%   (N).
%
%  -t1: the t1 image for this subject.  Either a path or an object will
%   sufffice
%
%  -saveDir:  the directory you would like these figures saved to.  If not
%   defined, saves to current directory.
%
%  -subSelect: a vector corresponding to the indexes of the tracts (in the
%   classification.names structure) which you would like to plot.  If this
%   is not defined, then the function will plot all classified fiber
%   tracts.
%
% (C) Lindsey Kitchell and Daniel Bullock, 2017, Indiana University


%% preliminaries

% loads requisite structures from input
[wbFG, fe] = bsc_LoadAndParseFiberStructure(wbFG);

%loads classificaiton file if a path is passed
if ischar(classification)
    load(classification);
end

% if an fe structure is detected, alters classificaiton index to only
% include positively weighted fibers
if ~isempty(fe)
    if length(fe.life.fit.weights)==length(classification)
        classification=wma_clearNonvalidClassifications(classification,fe);
    else
        warning('mismatch between classification structure and fe weights')
    end
end


% if user does not pass in a subselection
if notDefined('subSelect')
    subSelect=1:length(classification.names);
end


%% plotting preliminaries

% creates a structure containing an fg object for each fiber tract
% classification
tractStruc = bsc_makeFGsFromClassification(classification, wbFG);

% loops over the selected fibers and performs the volume generation
for ifg=1:length(subSelect)

    fiberBoolNifti=bsc_singleTractVolume(fg,20,1,[5 5 5],.5);
    
    %sets file names appropriately
    boolSaveName=strcat(saveDir',fg.name, '_Vol.nii.gz');
    boolSaveName = strrep(boolSaveName, ' ', '_');

    %saves volume as nifti
    niftiWrite(fiberBoolNifti,boolSaveName);
    fprintf(strcat('\n Done with ', tractStruc(ifg).name, '\n'))
end

end




