function fiberBoolNifti=bsc_singleTractVolume(fg,thresholdPercent,smoothBool,smoothKernel,voxelResize)

if notDefined('thresholdPercent')
thresholdPercent = 20;
end

if smoothBool && notDefined ('smoothKernel')
smoothKernel     = [3 3 3];
end

if notDefined('voxelResize')
voxelResize = 1;
end

%finding the minimum node coordinate for each streamline for each dimension
for ifibers=1:length(fg.fibers)
    minX(ifibers)=min(fg.fibers{ifibers},1);
    minY(ifibers)=min(fg.fibers{ifibers},2);
    minZ(ifibers)=min(fg.fibers{ifibers},3); 
end

%find absolute minimum value for entire FG for reach dimensiuon
absolutexMin=min(minX);
absoluteyMin=min(minY);
absolutezMin=min(minZ);

%establish offset necessary to make all coordinates positive and thus indexable (with a 5 mm window)
if absolutexMin<0 
xOffset=abs(absolutexMin)+floor(smoothKernel(1)/2)+5;
else
xOffset=-abs(absolutexMin)+floor(smoothKernel(1)/2)+5;
end

if absoluteyMin<0 
yOffset=abs(absoluteyMin)+floor(smoothKernel(2)/2)+5;
else
yOffset=-abs(absoluteyMin)+floor(smoothKernel(2)/2)+5;
end

if absolutezMin<0 
zOffset=abs(absolutezMin)+floor(smoothKernel(3)/2)+5;
else
zOffset=-abs(absolutezMin)+floor(smoothKernel(3)/2)+5;
end

%change the units of the smoothingKernel to the resized voxel units
smoothKernelResize=smoothKernel/voxelResize;

%round it to the nearest odd number
smoothKernelResize=2.*round((smoothKernelResize+1)./2)-1;

% chagne the coordinate scheme of the fg streamlines to positive only.
for ifibers=1:length(fg.fibers)
    fg.fibers{ifibers}(1,:)=fg.fibers{ifibers}(1,:)+xOffset;
    fg.fibers{ifibers}(2,:)=fg.fibers{ifibers}(2,:)+yOffset;
    fg.fibers{ifibers}(3,:)=fg.fibers{ifibers}(3,:)+zOffset;
end

%resize the streamlines to the coordinate scheme desired and find the bounding box coordinates
for ifibers=1:length(fg.fibers)
    fgResize.fibers{ifibers}=fg.fibers{ifibers}./voxelResize;
    
    resizeMinX(ifibers)=min(fgResize.fibers{ifibers},1);
    resizeMinY(ifibers)=min(fgResize.fibers{ifibers},2);
    resizeMinZ(ifibers)=min(fgResize.fibers{ifibers},3); 
    
    resizeMaxX(ifibers)=max(fgResize.fibers{ifibers},1);
    resizeMaxY(ifibers)=max(fgResize.fibers{ifibers},2);
    resizeMaxZ(ifibers)=max(fgResize.fibers{ifibers},3); 
end

%not useful
absolutexResizeMin=min(resizeMinX);
absoluteyResizeMin=min(resizeMinY);
absolutezResizeMin=min(resizeMinZ);

%find the max coords of the bounding box
absolutexResizeMax=max(resizeMaxX);
absoluteyResizeMax=max(resizeMaxY);
absolutezResizeMax=max(resizeMaxZ);

%create a matrix to hold a count of the endpoints in particular resized voxels
FiberVolume=zeros(floor(absolutexResizeMax),floor(absoluteyResizeMax),floor(absolutezResizeMax))

%go through each node for each streamline and count the number of nodes in each of the resized voxels
for ifibers=1:length(fgResize.fibers)
    for iNodes=1:length(fgResize.fibers{ifibers})
    %make sure this is going by columns
    coordinates=floor(fgResize.fibers{ifibers}(:,iNodes));
    FiberVolume(coordinates)=FiberVolume(coordinates)+1;
    end   
end


 if smoothBool == 0
        boolMatrixVersion   = FiberVolume>0;
        fiberBoolNifti.data = boolMatrixVersion;
    else
        smoothData = smooth3(FiberVolume,'gaussian',smoothKernelResize);
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

fiberBoolNifti.dim=size(fiberBoolNifti.data)
fiberBoolNifti.pixdim=voxelResize

