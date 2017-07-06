function bsc_plotClassifiedStreams(wbFG, classification, t1, view, saveDir,subSelect,colors)
%
%   bsc_plotClassifiedStreams(wbFG, classification, t1, view, saveDir,subSelect,colors)
%
%   PURPOSE: This function plots classified fibers using
%   mbaDisplayConnectome.  Will either plot all classified fiber tracts or
%   those that are subselected. Also, if you pass in an fe structure it
%   will only plot the validated fibers.  If you'd like to prune your
%   fibers as well, feel free to use removeOutliersClassification
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
%  -view: either 'coronal', 'axial', 'transverse', or 'saggital'.  Case
%   does not matter.  Creates views from both sides (i.e., top and bottom,
%   left and right, front and back).
%
%  -saveDir:  the directory you would like these figures saved to.  If not
%   defined, saves to current directory.
%
%  -subSelect: a vector corresponding to the indexes of the tracts (in the
%   classification.names structure) which you would like to plot.  If this
%   is not defined, then the function will plot all classified fiber
%   tracts.
%
%  -colors: an Nx3 matrix of values between 0 and 1 corresponding to the
%   RGB mapping desired of the corresponding fiber tracts (i.e. N=1
%   corresponds to fiber tract 1 in the classificaiton.names structure).
%   If this is not defined, the function will attempt determine if there is
%   a left/right pattern in the indexing and generate colors accordingly.
%   Barring all this, it will simply select random colors for each fiber
%   tract.
% 
% mba plotting component of function based on code origionally written by
% Franco Pestilli
%
% (C) Daniel Bullock, 2017, Indiana University


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
classification=wma_clearNonvalidClassifications(classification,fe);
end

%loads t1 if a path is passed
if ischar(t1)
    t1= niftiRead(t1);
end

% if user does not pass in a subselection
if notDefined('subSelect')
    subSelect=1:length(classification.names);
end

%defines colors if necessary
if notDefined('colors')
    if iseven (length(subSelect))
        if sum(subSelect(2:2:end)==subSelect(1:2:end)+1)==length(subSelect)/2
            % if there appears to be a pattern to the index choices, colors
            % will be assigned on the assumption of L/R pairings 
            for iTracts=2:2:length(subSelect)
                colors(iTracts,1)=rand;
                colors(iTracts,2)=rand;
                colors(iTracts,3)=rand;
                colors(iTracts-1,1)=colors(iTracts,1);
                colors(iTracts-1,2)=colors;
                colors(iTracts-1,3)=colors;
            end
        end
    else
        % otherwise, each fiber group gets a random color
        for iTracts=1:length(subSelect)
            colors(iTracts,1)=rand;
            colors(iTracts,2)=rand;
            colors(iTracts,3)=rand;
        end
    end
end

% save path is set to current directory if not defined
if notDefined('saveDir')
    saveDir=pwd;
end

%% plotting preliminaries

% creates a structure containing an fg object for each fiber tract
% classification
tractStruc = bsc_makeFGsFromClassification(classification, wbFG);

% starts a parpool
[pool, poolDir]=startUniqueParpool(8);

% defines the t1 slices that will serve as the background
slices      = {[-1 0 0],[0 -1 0],[0 0 -1]};

% figure setting
fh = figure('name','FiberClassificationPlot','color','k','units','normalized','position',[.5 .5 .5 .5]);
axis square
fhNum = fh.Number;
hold on

%set background according to input choice
switch lower(view)
    case 'saggital'
        h  = mbaDisplayBrainSlice(t1, slices{1});
    case 'coronal'
        h  = mbaDisplayBrainSlice(t1, slices{2});
    case {'axial', 'transverse'}
        h  = mbaDisplayBrainSlice(t1, slices{3});
end

%% plotting

%plot each of the fiber tracts iteratively
for itract = 1:length(subSelect)
    tractIndex=subSelect(itract);
    if exist('lh','var'), delete(lh); end
    if ~isempty(tractStruc(tractIndex).fg.fibers)
        [fh, lh] = mbaDisplayConnectome(tractStruc{itract}.fg.fibers,fhNum, colors{itract}, 'single');
        delete(lh)
        fprintf('\n %s \n',classification.names)
    end
end

% set views to the desired azimuth/elevation
switch lower(view)
    case 'saggital'
        fig.views = {[90,0],[-90,0]};
        light.angle = {[90,45],[-90,45]};
        fig.names = {strcat(saveHeaderData,'_ClassifiedTracks_saggital')};
    case 'coronal'
        fig.views = {[0,0],[180,0]};
        light.angle = {[0,45],[180,45]};
        fig.names = {strcat(saveHeaderData,'_ClassifiedTracks_Coronal')};
    case {'axial', 'transverse'}
        fig.views = {[0,90],[0,-90]};
        light.angle = {[0,45],[0,-45]};
        fig.names = {strcat(saveHeaderData,'_ClassifiedTracks_Axial')};
end

% iterate through view settings and save a picture for each one.
for iview = 1:length(fig.views)
    view(fig.views{iview,:})
    lh = camlight('left'); 
    lh.Position=light.angle{iview};
    axis square
    feSavefig(fhNum,'verbose','yes', ...
        'figName',strcat(fig.names,num2str(iview)), ...
        'figDir',saveDir, ...
        'figType','jpg');
    delete(lh)
end

close all

%deletes pool and temporary directory if it exists.
delete(pool);
if exists('poolDir','var')
    if ~isempty(poolDir)
        rmdir(poolDir)
    end
    fprintf('\n parpool directory deleted \n')
end


end