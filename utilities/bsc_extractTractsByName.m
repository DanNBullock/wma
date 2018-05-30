function subClassification = bsc_extractTractsByName(classification,names)
subClassification.names=[];
subClassification.index=zeros(size(classification.index));
classificationGrouping = wma_classificationStrucGrouping(classification);


for iNames=1:length(names)
    
    singleNameIndex=find(strcmpi(classificationGrouping.names,names{iNames}));
    
    singleStreamsIndexs=find(classificationGrouping.index==singleNameIndex);
    multiStreamIndexs=unique(classification.index(singleStreamsIndexs));
    
    for iTracts=1:length(multiStreamIndexs)
        curIndex=max(unique(subClassification.index))+1;
        subClassification.names{curIndex}= classification.names{multiStreamIndexs(iTracts)};
        subClassification.index(classification.index==multiStreamIndexs(iTracts))=curIndex;
    end
end

end
    