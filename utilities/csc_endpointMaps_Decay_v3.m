function  [nii, nii2, nii_normalized, nii2_normalized] = csc_endpointMaps_Decay_v3(fg, fsDir,saveDir, saveHeader, thresholdinmm, outlierreject, decayFunc)
% function  csc_endpointMaps(path2t1,path2graynii,fg,saveDir)
%
% OVERVIEW: generate .nifti files for (both) of the tract endpoint density
% mappings for a given fibergroup and save them down in the designated
% output directory.
%
%
% INPUTS:
% -fg: fiber group structure, in an ACPC coordinate scheme, as a transform
% to an img space coordinate scheme is later performed
% -fsDir: path to THIS SUBJECT'S freesurfer directory
% -saveDir: The designared output directory
% -thresholdinmm: A distance threshold between tract endpoint and gray
% matter voxels (default: 2 mm).
% -outlierreject: True, perform outlier rejection; false, do not perform
% outlier rejection. Default is false.
% -decayFunc: the decay function to use to calculate weight to throw in
% voxels 
%          ->uniform: no distance loss until threshold, then zero
%          ->linear: normalized linear cost to distance, 0 at dist > threshold
%          ->exponential: exponential cost to distance, 0 at dist > threshold
%          ->exact:  only counts exact endpoint voxel, uses floor()
%
%
% REQUIREMENT: Parallel computing toolbox and MATLAB newer than 2013b.
%
% OUTPUTS:
% -nii, nii2: Endpoint density file for two distinct endpoint clouds of the tract.
% In this file, we simply sum the values across streamline endpoints.
% -nii_normalized, nii2_normalized: Endpoint density with normalization, 
% based on max endpoint value.
%
% 

% (C) Daniel Bullock and Hiromasa Takemura, 2016, CiNet HHS
%  Decay functions added by DNB 05/2017

%% Parameter settings & preliminaries 
if notDefined('thresholdinmm'), thresholdinmm=2;end
% Define a distance threshold. Default is 2 mm.
if notDefined('outlierreject'), outlierreject=0;end

if notDefined('decayFunc'), decayFunc='uniform';end

%generate the gm mask if it doesn't exist
    if ~exist(strcat(fsDir,'/gm_mask.nii.gz'))
        %apaprently necessary for matlab?
        spaceChar={' '};
        %why are these structures?
        string1=strcat('mri_binarize --i',spaceChar,fsDir,'mri/aparc+aseg.mgz --gm --o',spaceChar, fsDir, 'gm_mask.mgz');
        string2=strcat('mri_convert',spaceChar,fsDir,'gm_mask.mgz',spaceChar, fsDir, 'gm_mask.nii.gz');
        [status1 result] = system(string1{1});
        [status2 result] = system(string2{1});
        if or(status1,status2)
            fprintf('/n Error generating gm nifti file.  There may be a problem finding the aaparc+aseg.mgz or the gm_mask.mgz file.')
            keyboard
        end
    end


%% Perform outlier rejection if desired
if outlierreject == 1,
    %perform fiber pruning if desired
    [fg, ~]=mbaComputeFibersOutliers(fg, 4, 4);
else
end

%% Load files

% DONT DO THIS
% % Read T1 file and use it to get the size of voxels in the image
% t1 = niftiRead(path2t1);
% xform = t1.qto_ijk;
% thresholdVoxel = thresholdinmm/t1.pixdim(1); % Define the threshold in voxel space of T1

% Load Gray matter mask; previously this was either a left or right mask,
% but now we make a whole brain grey matter mask.  As such we can probably
% handle cross-hemispheric fibers.
graynii = niftiRead(strcat(fsDir,'/gm_mask.nii.gz'));

%% Estimate the orientation of the tract

%convert fg to imgspace coordinate scheme

%fg = dtiXformFiberCoords(fg, xform, 'img');



%calculate the net traversal of the fiber group in each dimension,
%which can then be used to determine the primary dimension of
%orientation for the fiber group
for ifibers=1:length(fg.fibers)
    xDisplacement(ifibers)=abs(fg.fibers{ifibers}(1,1)-fg.fibers{ifibers}(1,end));
    yDisplacement(ifibers)=abs(fg.fibers{ifibers}(2,1)-fg.fibers{ifibers}(2,end));
    zDisplacement(ifibers)=abs(fg.fibers{ifibers}(3,1)-fg.fibers{ifibers}(3,end));
end
displaceVec=[mean(xDisplacement) mean(yDisplacement) mean(zDisplacement)];
tractorientation=find(displaceVec==max(displaceVec)); % Find out a maximum dimenstion from x, y, z to check a predominant orientation of the tract

% Generate file names (tractorientation = 1, left-right; 2,
% anterior-posterior; 3, superior-inferior)
switch tractorientation
    case 1
        orientation1='right';
        orientation2='left';
    case 2
        orientation1='anterior';
        orientation2='posterior';
    case 3
        orientation1='superior';
        orientation2='inferior';
end


%% Set output file names
fname{1} = strcat(saveHeader,'_',fg.name,'_',decayFunc,'_',num2str(thresholdinmm),'mm_',orientation1,'FiberEndpoint.nii.gz'); % Tract endpoint density without normalization
fname{2} = strcat(saveHeader,'_',fg.name,'_',decayFunc,'_',num2str(thresholdinmm),'mm_',orientation2,'FiberEndpoint.nii.gz');
saveFileName{1} =  strcat(saveHeader,'_',fg.name,'_',decayFunc,'_',num2str(thresholdinmm),'mm_',orientation1,'FiberEndpoint_normalized.nii.gz'); % Tract endpoint density with normalization
saveFileName{2} =  strcat(saveHeader,'_',fg.name,'_',decayFunc,'_',num2str(thresholdinmm),'mm_',orientation2,'FiberEndpoint_normalized.nii.gz');

%% Extract fiber endpoint coordinates
fbsize = size(fg.fibers);

% Determine number of 1000 streamlines sub-groups the fg can be divided to improve
% the computational speed. If the total number of streamlines are less than
% thousand, all streamlines are going into a single group.
fiberBins=0:1000:fbsize(1);
if fiberBins(end)~=fbsize(1)
    fiberBins=horzcat(fiberBins,fbsize(1)); % Set bin size, just in case number of streamlines are less than 1000
end

%% generate blank nii structures for the endpoint mappings
nii2 = graynii;
nii2.fname = fullfile(saveDir,fname{2});
nii2.data = zeros(size(graynii.data));

nii = graynii;
nii.fname = fullfile(saveDir,fname{1});
nii.data = zeros(size(graynii.data));

fprintf('\n Fiber group %s contains %i fibers and will therefore require %i iteration(s) ',fg.name,fbsize(1) , length(fiberBins)-1)
fprintf('\n Iteration ')
tic

    %% Get Gray coordinates apropriate coordinate scheme.
    % Here we return a linear indexing of each voxel of the grey matter mask
    graynii.data(1,1,1)=1;
    [gray_coords(1,:), gray_coords(2,:), gray_coords(3,:)] = ind2sub(size(graynii.data),find(graynii.data));
    
    [CenterOfMass_gray_coords, outMat, imScale] = mrAnatXformCoords(graynii.qto_xyz, gray_coords);
%     scatter3(CenterOfMass_gray_coords(:,1), CenterOfMass_gray_coords(:,2), CenterOfMass_gray_coords(:,3))
%     view([0,90])
    
    CenterOfMass_gray_coords=CenterOfMass_gray_coords';
    
    %% parpool
        %create a parpool object if one does not already exist
    poolobj = gcp('nocreate'); % If no pool, do not create new one.
    if isempty(poolobj)
        parpool
    else
    end

%% Iteratively compute the density maping for each 1000 fiber sub-group
for iFibBins=1:length(fiberBins)-1
    fprintf('%i _',iFibBins)
    %determine the lower and upper bin bounds for this iteration
    
    lowBin=fiberBins(iFibBins);
    highBin=fiberBins(iFibBins+1);
    
    for ph = 1:fiberBins(iFibBins+1)-fiberBins(iFibBins)
        %convert the cell object fgimg.fibers into a standard matrix (I
        %guess?)
        f_coords = cell2mat(fg.fibers(fiberBins(iFibBins)+ph));
        %determine the endpoint of the fiber
        fcoordsize = size(f_coords);
        
        %for each fiber, determine how the fiber is oriented (i.e. if
        %the first indexed coordinate is the most
        %left/anterior/superior), some preprocessing steps (i.e.
        %bsc_reorientfibers) do this automatically.
        if f_coords(tractorientation,1)>f_coords(tractorientation,fcoordsize(2)) % Classify two endpoints based on its positions
            fendpoints(1, ph,1) = f_coords(1,1);
            fendpoints(2, ph,1) = f_coords(2,1);
            fendpoints(3, ph,1) = f_coords(3,1);
            fendpoints(1, ph,2) = f_coords(1,fcoordsize(2));
            fendpoints(2, ph,2) = f_coords(2,fcoordsize(2));
            fendpoints(3, ph,2) = f_coords(3,fcoordsize(2));
        else
            fendpoints(1, ph,1) = f_coords(1,fcoordsize(2));
            fendpoints(2, ph,1) = f_coords(2,fcoordsize(2));
            fendpoints(3, ph,1) = f_coords(3,fcoordsize(2));
            fendpoints(1, ph,2) = f_coords(1,1);
            fendpoints(2, ph,2) = f_coords(2,1);
            fendpoints(3, ph,2) = f_coords(3,1);
        end
    end
    
    %next generate a bounding matrix that represents the minimum and
    %maxium node coordinate value for each dimension.  Basically
    %this creates a 3 dimensional window that we look for grey matter
    %in, greatly reducing our computational overhead in the subsequent
    %steps.
    for iDimensions=1:3
        for iClusters=[1 2]
            %interpretation: 1st = dimension (i.e. 1 = x, 2 =y, 3=z) , 2nd = cluster
            %(i.e. left/anterior/superior vs right/posterior/inferior endpoints), 3rd = min/max (i.e. 1 = min, 2 =
            %max)
            endMinMax(iDimensions,iClusters,1)=min(fendpoints(iDimensions,:,iClusters));
            endMinMax(iDimensions,iClusters,2)=max(fendpoints(iDimensions,:,iClusters));
            
        end
    end
    
    %expand the window so that you include the maximally distanced gray
    %matter, prior to this step, we were only considering exactly where the fibers
    %ended and went no further.  Also, change min
    endMinMax(:,:,1)=endMinMax(:,:,1)-thresholdinmm;
    endMinMax(:,:,2)=endMinMax(:,:,2)+thresholdinmm;
        
    % Restrict the gray matter coordinate to chose a voxel within a certain
    % distance from streamline endpoint to improve the computational
    % efficiency.
    criteriaD1=CenterOfMass_gray_coords(1,:) < endMinMax (1,1,2);
    criteriaD2=CenterOfMass_gray_coords(1,:) > endMinMax (1,1,1);
    criteriaD3=CenterOfMass_gray_coords(2,:) < endMinMax (2,1,2);
    criteriaD4=CenterOfMass_gray_coords(2,:) > endMinMax (2,1,1);
    criteriaD5=CenterOfMass_gray_coords(3,:) < endMinMax (3,1,2);
    criteriaD6=CenterOfMass_gray_coords(3,:) > endMinMax (3,1,1);
    
    criteriaV1=CenterOfMass_gray_coords(1,:) < endMinMax (1,2,2);
    criteriaV2=CenterOfMass_gray_coords(1,:) > endMinMax (1,2,1);
    criteriaV3=CenterOfMass_gray_coords(2,:) < endMinMax (2,2,2);
    criteriaV4=CenterOfMass_gray_coords(2,:) > endMinMax (2,2,1);
    criteriaV5=CenterOfMass_gray_coords(3,:) < endMinMax (3,2,2);
    criteriaV6=CenterOfMass_gray_coords(3,:) > endMinMax (3,2,1);
    
    % Find those gray matter voxel linear indicies that meet *all* three criterion set forth by the endMinMax vector
    validDGrayCoords=find(criteriaD1 & criteriaD2 & criteriaD3 & criteriaD4 & criteriaD5 & criteriaD6);
    validVGrayCoords=find(criteriaV1 & criteriaV2 & criteriaV3 & criteriaV4 & criteriaV5 & criteriaV6);
    
    %% Get distance between streamline endpoint and gray matter voxels
    %Initialize data structures for the outputs for the next several steps
    bestSqDist_D=zeros(fiberBins(iFibBins+1)-fiberBins(iFibBins),length(validDGrayCoords));
    bestSqDist_V=zeros(fiberBins(iFibBins+1)-fiberBins(iFibBins),length(validVGrayCoords));
    indices_D=zeros(fiberBins(iFibBins+1)-fiberBins(iFibBins),length(validDGrayCoords));
    indices_V=zeros(fiberBins(iFibBins+1)-fiberBins(iFibBins),length(validVGrayCoords));
    

    
    % Within a parfor process compute the distance between each streamline
    % endpoint and each grey matter voxel within the 3D window. Previously,
    % without the windowing operation, this computation would occur for
    % every voxel of the white matter mask.  Further computational speed
    % up was obtained by the binning process, which reduced the number
    % of fibers this was computed for *and* resulted in a smaller 3D
    % window
    
    
    parfor kk = 1:fiberBins(iFibBins+1)-fiberBins(iFibBins)
        [indices_D(kk,:), bestSqDist_D(kk,:)] = nearpoints(CenterOfMass_gray_coords(:,validDGrayCoords), fendpoints(:,kk,1));
        [indices_V(kk,:), bestSqDist_V(kk,:)] = nearpoints(CenterOfMass_gray_coords(:,validVGrayCoords), fendpoints(:,kk,2));
    end
    
    %create containers for the distances between each fiber endpoint
    %and valid grey matter voxel.  Note: the dorsal/ventral naming
    %conventions are just a holdover from a previous iteration of the
    %code.  They do not refer to the actual orientation of the fiber or
    %the... something
    
    exact_D=floor(fendpoints(:,:,1));
    
    exact_V=floor(fendpoints(:,:,2));
    
    %structure to hold fiber distance from grey matter voxel (2nd dimension) 
    % which is within the threshold distance from a given endpoint (1st Dim)
    dorsalvox_indiceDist = zeros(fiberBins(iFibBins+1)-fiberBins(iFibBins),length(gray_coords));
    ventralvox_indiceDist = zeros(fiberBins(iFibBins+1)-fiberBins(iFibBins),length(gray_coords));
    
    %structure to hold a 1 (true) in those grey matter voxels (2nd dimension)
    % which are within the threshold distance of a given fiber endpoint (1st Dim)
    dorsalvox_indiceBool = zeros(fiberBins(iFibBins+1)-fiberBins(iFibBins),length(gray_coords));
    ventralvox_indiceBool = zeros(fiberBins(iFibBins+1)-fiberBins(iFibBins),length(gray_coords));
    
    %same for other coresponding structures
    dorsalvox_lin_indexed = zeros(fiberBins(iFibBins+1)-fiberBins(iFibBins),length(gray_coords));
    ventralvox_lin_indexed = zeros(fiberBins(iFibBins+1)-fiberBins(iFibBins),length(gray_coords));
    dorsalvox_exp_indexed = zeros(fiberBins(iFibBins+1)-fiberBins(iFibBins),length(gray_coords));
    ventralvox_exp_indexed = zeros(fiberBins(iFibBins+1)-fiberBins(iFibBins),length(gray_coords));
    
    %an x by y array of integers i, where x=the number of fibers, y =
    %number of within window grey matter voxels and i = the radius
    %threshold established by the input thresholdinmm parameter.
    blankDistVec_D=indices_D*thresholdinmm;
    blankDistVec_V=indices_V*thresholdinmm;
    
    % an an x by y (as established above) boolean vector containing 1 where
    % the endpoint to grey matter voxel distance is less than the
    % threshold, and 0 everywhere else
    validEntries_D=bestSqDist_D<thresholdinmm;
    validEntries_V=bestSqDist_V<thresholdinmm;
    
    %linearly normalized array where each nonzero entry equals the grey
    % matter voxel's distance from the associated endpoint, such that a
    % distance of zero =1 and a distance approaching the threshold
    % approaches zero in a linear fashion
    LinearNormDist_D=validEntries_D.*((blankDistVec_D-bestSqDist_D)./blankDistVec_D);
    LinearNormDist_V=validEntries_V.*((blankDistVec_V-bestSqDist_V)./blankDistVec_V);
    dorsalvox_lin_indexed(:,validDGrayCoords)=LinearNormDist_D(:,:);
    ventralvox_lin_indexed(:,validVGrayCoords)=LinearNormDist_V(:,:);
    
    %exponential decay array where each nonzero entry equals the grey
    % matter voxel's distance from the associated endpoint, such that a
    % distance of zero =1 and a distance approaching the threshold
    % approaches zero in an exponential fashion
    exponentialDecay_D=validEntries_D.*((exp(blankDistVec_D-bestSqDist_D)));
    exponentialDecay_D(exponentialDecay_D~=0)=1./exponentialDecay_D(exponentialDecay_D~=0);
    exponentialDecay_V=validEntries_V.*((exp(blankDistVec_V-bestSqDist_V)));
    exponentialDecay_V(exponentialDecay_V~=0)=1./exponentialDecay_V(exponentialDecay_V~=0);
    dorsalvox_exp_indexed(:,validDGrayCoords)=exponentialDecay_D(:,:);
    ventralvox_exp_indexed(:,validVGrayCoords)=exponentialDecay_V(:,:);
    
    % Match the matrix size to original gray matter coordinate structure
    dorsalvox_indiceDist(:,validDGrayCoords)=bestSqDist_D(:,:);
    ventralvox_indiceDist(:,validVGrayCoords)=bestSqDist_V(:,:);
    
    % same for exact
    %actually, turns out you dont need any code for this.
    %     dorsalvox_exact_indexed(:,validDGrayCoords)=
    %     ventralvox_exact_indexed(:,validVGrayCoords)
    
    %generates the structure (dorsal/ventral)vox_indiceBool using the
    %indiceDist matrix from the previous step.  The first part of the
    %conditional ensures that the relevant endpoints are within the set
    %threshold distance.  The second part ensures that all zero
    %distance endpoints (i.e. all of the distances which were not
    %computed, but are zero due to the initial generatin of the
    %matricies) are ignored.  In theory, this could cause a problem if
    %an endpoint was *exactly* the value of a mesh vertex, but because
    %the mesh points are created at regular distances with a limited
    %number of decimal points (relative to fiber node coordinates, at
    %least) this probably will not ever happen.
    dorsalvox_indiceBool(dorsalvox_indiceDist<thresholdinmm & dorsalvox_indiceDist > 0) = 1;
    ventralvox_indiceBool(ventralvox_indiceDist<thresholdinmm & ventralvox_indiceDist > 0) = 1;
    
    %sums across the fibers in order to compute the total number of endpoints
    %that are within the threshold value of each grey matter coordinate
    dorsalvox_fbdensity_Uniform = sum(dorsalvox_indiceBool,1);
    ventralvox_fbdensity_Uniform = sum(ventralvox_indiceBool,1);
    
    dorsalvox_fbdensity_LinearDecay = sum(dorsalvox_lin_indexed,1);
    ventralvox_fbdensity_LinearDecay = sum(ventralvox_lin_indexed,1);

    dorsalvox_fbdensity_ExponentialDecay = sum(dorsalvox_exp_indexed,1);
    ventralvox_fbdensity_ExponentialDecay = sum(ventralvox_exp_indexed,1);

    %end decay section
    
    
    %% add current counts/sums to nii data
    %add the current density data to the existing nii.data structure.
    %This can be acheieved thorugh simple summation as it is a density
    %count.
    switch decayFunc
        case 'uniform'
            for i = 1:length(dorsalvox_fbdensity_Uniform)
                nii.data(gray_coords(1,i), gray_coords(2,i), gray_coords(3,i)) =nii.data(gray_coords(1,i), gray_coords(2,i), gray_coords(3,i))+ dorsalvox_fbdensity_Uniform(:,i);
            end
            
            for i = 1:length(ventralvox_fbdensity_Uniform)
                nii2.data(gray_coords(1,i), gray_coords(2,i), gray_coords(3,i)) =nii2.data(gray_coords(1,i), gray_coords(2,i), gray_coords(3,i))+  ventralvox_fbdensity_Uniform(:,i);
            end
        case 'linear'
            for i = 1:length(dorsalvox_fbdensity_LinearDecay)
                nii.data(gray_coords(1,i), gray_coords(2,i), gray_coords(3,i)) =nii.data(gray_coords(1,i), gray_coords(2,i), gray_coords(3,i))+ dorsalvox_fbdensity_LinearDecay(:,i);
            end
            
            for i = 1:length(ventralvox_fbdensity_LinearDecay)
                nii2.data(gray_coords(1,i), gray_coords(2,i), gray_coords(3,i)) =nii2.data(gray_coords(1,i), gray_coords(2,i), gray_coords(3,i))+  ventralvox_fbdensity_LinearDecay(:,i);
            end
        case 'exponential'
            for i = 1:length(dorsalvox_fbdensity_ExponentialDecay)
                nii.data(gray_coords(1,i), gray_coords(2,i), gray_coords(3,i)) =nii.data(gray_coords(1,i), gray_coords(2,i), gray_coords(3,i))+ dorsalvox_fbdensity_ExponentialDecay(:,i);
            end
            
            for i = 1:length(ventralvox_fbdensity_ExponentialDecay)
                nii2.data(gray_coords(1,i), gray_coords(2,i), gray_coords(3,i)) =nii2.data(gray_coords(1,i), gray_coords(2,i), gray_coords(3,i))+  ventralvox_fbdensity_ExponentialDecay(:,i);
            end
        case 'exact'
            for i = 1:length(exact_D)
                nii.data(exact_D(1,i),exact_D(2,i),exact_D(3,i)) =nii.data(exact_D(1,i),exact_D(2,i),exact_D(3,i)) +1;
            end
            
            for i = 1:length(exact_V)
                nii.data(exact_V(1,i),exact_V(2,i),exact_V(3,i)) =nii.data(exact_V(1,i),exact_V(2,i),exact_V(3,i)) +1;
            end
    end
    
    %clearing variables that are generated by indexing
    clear validDGrayCoords
    clear validVGrayCoords
    clear bestSqDist_D
    clear bestSqDist_V
    clear fendpoints
    
end
computeTime=toc;
fprintf('\n Density mesh generation complete in %i seconds',computeTime)

%% Save file
%writing both the raw and normalized density meshes.
niftiWrite(nii);
niftiWrite(nii2);

[nii_normalized] = csc_normalizedensitynifti(nii,fullfile(saveDir,saveFileName{1}));
[nii2_normalized] = csc_normalizedensitynifti(nii2,fullfile(saveDir,saveFileName{2}));
end

