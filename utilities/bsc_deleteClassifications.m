function [classification] = bsc_deleteClassifications(baseClassification,removeIndexes)





removeIndexes=sort(removeIndexes,'descend');



baseNames=baseClassification.names;
baseNameNum=length(baseClassification.names);

for iRemove=1:length(removeIndexes)
    indx2remove=find(baseClassification.index==removeIndexes(iRemove));
    baseClassification.index(indx2remove)=0;
    baseClassification.names(removeIndexes(iRemove))=[];
end

classification=baseClassification;

end
    