function [L_fg_vert, R_fg_vert, L_vertical_fascicles_identities, R_vertical_fascicles_identities] = wma_find_vertical_fibers(fgPath,fsDir,outdir,thresh, v_crit, minLength)
% Find vertically oriented fibers projecting to VOT
%
% [L_fg_vert, R_fg_vert, L_vertical_fascicles_identities, R_vertical_fascicles_identities] = AFQ_FindVerticalFibers(fgPath,fsROIdir,outdir,thresh,v_crit, minLength)
%
% Inputs:
%
% fgPath    - Patht to a whole brain connectome (or fg structure).
%             This will be segmented
% fsROIdir  - path to THIS SUBJECT'S fs directory
% outdir    - Directory to save output
% thresh    - 
% v_crit    - 
% minLength - minimum distance between ROIs and fascicles
% 
% Intiial version of code by Jason D. Yeatman, September 2014. 
%
% Code released with:
% Yeatman J.D., Weiner K.S., Pestilli F., Rokem A., Mezer A., Wandell B.A.
% (2014). The vertical occipital fasciculus: A forgotten highway. PNAS.
%
% Modified by Daniel Bullock 2017, Indiana University

%% Set parameters

% Remove any fiber that doesn't go vertical (positive z) for thresh% of its
% coordinates
if notDefined('thresh')
    thresh = [.95 .6];
end
% Fibers must travel this much farther vertically than other directions
if notDefined('v_crit')
    v_crit = 1.3;
end
% Minumum length
if notDefined('minLength')
    minLength = 20;
end
%% Set paths and load fibers
if notDefined('antBoundary')
    antBoundary = -15;
end
roi_names = {'fusiform.mat' 'inferiortemporal.mat'...
    'lateraloccipital.mat'}; %'middletemporal.mat'

% Path to ROIs
if notDefined('fsDir')
    fsROIdir = uigetdir([],'Select the appropriate FreeSurfer subject directory');
end

% Path to fibers
if notDefined('fgPath')
    [fname,pname]=uigetfile('*','Select Fiber Group',[],'MultiSelect','on');
    if iscell(fname)
        for ii=1:length(fname)
            fgPath{ii} = fullfile(pname{ii},fname{ii});
        end
    else
        fgPath = fullfile(pname,fname);
    end
end
% output directory
if notDefined('outdir')
    outdir = uigetdir([],'Select an output directory');
end

% Load fibers. If there are multiple fiber groups merge them into 1
if iscell(fgPath)
    fg = fgRead(fgPath{1});
    for ii = 2:length(fgPath)
        % Load and merge each fiber group
        fgtmp = fgRead(fgPath{ii});
        fg    = dtiMergeFiberGroups(fg,fgtmp);
    end
elseif ischar(fgPath)
    fg = fgRead(fgPath);
elseif isstruct(fgPath);
    fg = fgPath;
end

% Remove fibers that are less than 2cm
% CHECK CHECK CHECK
L  = cellfun(@(x) length(x),fg.fibers);
ac = cellfun(@(x) max(x(2,:)),fg.fibers);
fasciles_long_indices      = L > minLength; %#ok<*NASGU>
fascicles_anterior_indices = ac < antBoundary;%was done in terms of length instead of antBoundary
fascicles_indices          = and(fasciles_long_indices, fascicles_anterior_indices );
fascicles_identities       = find(fascicles_indices);
fg.fibers                  = fg.fibers(fascicles_indices);

%% Create VOT roi from freesurfer ROIs
% We select a few FreeSurfer ROIs from the segmentation.
% Left Hemisphere Rois -> [1007,1009,1011]
% Right Hemisphere Rois -> [2007,2009,2011]
[L_roi_all] = wma_roiFromFSnums(fsDir,[1007,1009,1011]);
[R_roi_all] = wma_roiFromFSnums(fsDir,[2007,2009,2011]);

% set names
L_roi_all.name = 'LVOT';
R_roi_all.name = 'RVOT';

%% Find all vertical fibers projecting to VOT
% Intersect fibers with ROI
[L_fg, ~, L_vertical_fascicles_indices] = dtiIntersectFibersWithRoi([],{'and' 'endpoints'},4,L_roi_all,fg);
L_vertical_fascicles_identities         = fascicles_identities(L_vertical_fascicles_indices);

[R_fg, ~, R_vertical_fascicles_indices] = dtiIntersectFibersWithRoi([],{'and' 'endpoints'},4,R_roi_all,fg);
R_vertical_fascicles_identities         = fascicles_identities(R_vertical_fascicles_indices);

% Flip each fiber so that the first point is the most inferior one
for ii = 1:length(L_fg.fibers)
   if L_fg.fibers{ii}(3,1) >  L_fg.fibers{ii}(3,end)
      L_fg.fibers{ii} = L_fg.fibers{ii}(:,end:-1:1); 
   end
end
for ii = 1:length(R_fg.fibers)
   if R_fg.fibers{ii}(3,1) >  R_fg.fibers{ii}(3,end)
       R_fg.fibers{ii} = R_fg.fibers{ii}(:,end:-1:1);
   end
end

% Compute the difference between each fiber coordinate and the previous one
f_diffL = cellfun(@(x) diff(x,[],2), L_fg.fibers, 'uniformoutput',0);
f_diffR = cellfun(@(x) diff(x,[],2), R_fg.fibers, 'uniformoutput',0);

% Compute the total distance traveled in each direction
f_distL = cellfun(@(x) sum(abs(x),2), f_diffL, 'uniformoutput',0);
f_distR = cellfun(@(x) sum(abs(x),2), f_diffR, 'uniformoutput',0);

% For each fiber coordinate see if it is superior/inferior, left/right, 
% ant/post to the previous one
f_dirL = cellfun(@(x) x > 0,f_diffL, 'uniformoutput',0);
f_dirR = cellfun(@(x) x > 0,f_diffR, 'uniformoutput',0);

%% Loop over fibers and remove ones that are not vertical

L_fg_vert = dtiNewFiberGroup('L_vertical');
R_fg_vert = dtiNewFiberGroup('R_vertical');


for ii = 1:length(L_fg.fibers)
    % Compute the % of the fiber length that is consistant in a cardinal
    % direction
    f_dir_percL = sum(f_dirL{ii},2)./size(L_fg.fibers{ii},2); 
    
    % Compute the % of the fiber length that travels more in the z
    % direction than other directions
    f_z_percL = sum(abs(f_diffL{ii}(3,:)) > abs(f_diffL{ii}(2,:)) & abs(f_diffL{ii}(3,:)) > abs(f_diffL{ii}(1,:)))./size(L_fg.fibers{ii},2); 
   
    % Remove fibers that don't:
    % (1) go vertical for long enough
    % (2) travel farther in the z direction than y direction
    % (3) travel farther in the z direction than x direction
    if f_dir_percL(3) > thresh(1) && ...
       f_z_percL      > thresh(2) && ...
       f_distL{ii}(3) > v_crit*f_distL{ii}(2) && ...
       f_distL{ii}(3) > v_crit*f_distL{ii}(1)
       L_fg_vert.fibers = vertcat(L_fg_vert.fibers,L_fg.fibers{ii});
    else
       L_vertical_fascicles_identities(ii) = 0;
    end
end

for ii = 1:length(R_fg.fibers)
    f_dir_percR = sum(f_dirR{ii},2)./size(R_fg.fibers{ii},2);
    f_z_percR   = sum(abs(f_diffR{ii}(3,:)) > abs(f_diffR{ii}(2,:)) & abs(f_diffR{ii}(3,:)) > abs(f_diffR{ii}(1,:)))./size(R_fg.fibers{ii},2); 

    if f_dir_percR(3) > thresh(1) && ...
       f_z_percR      > thresh(2) && ...
       f_distR{ii}(3) > v_crit*f_distR{ii}(2) && ...
       f_distR{ii}(3) > v_crit*f_distR{ii}(1)
       R_fg_vert.fibers = vertcat(R_fg_vert.fibers,R_fg.fibers{ii});
    else
       R_vertical_fascicles_identities(ii) = 0;
    end
end

L_vertical_fascicles_identities = L_vertical_fascicles_identities( (L_vertical_fascicles_identities~=0) );
R_vertical_fascicles_identities = R_vertical_fascicles_identities( (R_vertical_fascicles_identities~=0) );

%% Save fibers
if length(L_fg_vert.fibers) > 10
    dtiWriteFiberGroup(L_fg_vert,fullfile(outdir,'Left_VerticalFG')); L = 1;
else
    L = 0;
end
if length(R_fg_vert.fibers) > 10
    dtiWriteFiberGroup(R_fg_vert,fullfile(outdir,'Right_VerticalFG')); R = 1;
else
    R = 0;
end

return