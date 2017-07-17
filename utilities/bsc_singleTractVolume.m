function fiberBoolNifti=bsc_singleTractVolume(fg,t1)


thresholdPercent = 20;
%islandFlag       = false;
smoothKernel     = [3 3 3];
voxelResize = config.voxelResize;

for ifibers=1:length(fg.fibers)
    fg.fibers{ifibers}=mrAnatXformCoords(t1.qto_ijk, fg.fibers{ifibers})';
    
end

    fiberBoolNifti = t1;
    imgDim         = t1.dim(1:3);
    imgRes         = t1.pixdim(1,1);
    imgResize      = imgDim*imgRes;
    imgBins        = ceil(imgResize/voxelResize);
    emptyMatrix    = zeros(imgBins);



    % adjusts nifti object field information (probably not necessary)
    fiberBoolNifti.dim       = imgBins;
    fiberBoolNifti.pixdim    = [voxelResize voxelResize voxelResize];
    fiberBoolNifti.qoffset_x = fiberBoolNifti.qoffset_x*(imgRes/voxelResize);
    fiberBoolNifti.qoffset_y = fiberBoolNifti.qoffset_y*(imgRes/voxelResize);
    fiberBoolNifti.qoffset_z = fiberBoolNifti.qoffset_z*(imgRes/voxelResize);


% Here we begin counting nodes per (resized) voxel
    for ifibers = 1:length(fg.fibers)
        roundedExpandedTract = round(fg.fibers{ifibers}*(imgRes/voxelResize));
        for inodes = 1:length(fg.fibers{ifibers})
            emptyMatrix(  roundedExpandedTract(1,inodes), ...
                          roundedExpandedTract(2,inodes), ...
                          roundedExpandedTract(3,inodes))= ...
            emptyMatrix(  roundedExpandedTract(1,inodes), ...
                          roundedExpandedTract(2,inodes), ...
                          roundedExpandedTract(3,inodes)) +1;
        end
    end


    %% actually do smoothing if necessary
    % smooths if necessary.  In theory de-islands as well, though I've never
    % seen it in action.  It is set to go to keyboard if it detects islands.
    % It could be that the function to detect islands is broken, so this may
    % simply be pointless.
    if config.smooth == 0
        boolMatrixVersion   = emptyMatrix>config.threshold;
        fiberBoolNifti.data = boolMatrixVersion;
    else
        smoothData = smooth3(emptyMatrix,'gaussian',smoothKernel);
        % auto threshold computation
        % computes the appropriate threshold for the given percentile value.
        % Probably not computationally efficient, but principled (aside from the
        % arbitrary threshold) and adaptive to the specific case.
        uniqueVals = unique(smoothData);

        % calculate bin interval according to uniquevals
        binz       = 0:(uniqueVals(end)/10000):uniqueVals(end);

        % is there something wrong here with the calculation of the second bin?
        hisDist    = histcounts(smoothData,binz);

        % count of nonzero entries, used for percentile calculation later
        nonZeroTotal = sum(hisDist(2:end));

        % loop to calculate cumulative sum of distribution
        for ibins = 1:length(hisDist)
            if ibins == 1
                cumulativeSums(ibins) = 0;
            else
                cumulativeSums(ibins) = cumulativeSums(ibins-1)+hisDist(ibins);
            end
        end
        % calculate percentiles
        percentiles = (cumulativeSums / nonZeroTotal)*100;
        failureVec  = (percentiles > thresholdPercent);

        % find corresponding bin value
        thresholdBinVal = binz(find(failureVec, 1 ));

        %the bwislands is currently returning an empty output for me, at least
        %on this test fiber.   I'm not sure if that means there are no islands
        %to be found for this fiber, or if the function doesn't work.  This
        %part may cause a problem at some point.
%         islands = bwislands(smoothData);
%         if ~isempty(islands)
%             islandFlag = true;
%             smoothData = deislands3d(smoothData,3);
%         end

        boolMatrixVersion=smoothData>thresholdBinVal;
        fiberBoolNifti.data=boolMatrixVersion;
    end

    % convert to unit8 from bool (necessary for niftiSave function)
    fiberBoolNifti.data = uint8(fiberBoolNifti.data);
    
end
