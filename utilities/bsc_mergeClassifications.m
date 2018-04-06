function [classification] = bsc_mergeClassifications(baseClassification,classificationAdd)









baseNames=baseClassification.names;
baseNameNum=length(baseClassification.names);

addNames=classificationAdd.names;
addNameNum=length(classificationAdd.names);

presumeNameNum=baseNameNum+addNameNum;

uniqueNamesTotal=unique(horzcat(baseNames,addNames));
uniqueNamesLength=length(uniqueNamesTotal);

if uniqueNamesLength==presumeNameNum
    %hyper inelegant, guarenteed to cause problems.
    classification.names=horzcat(baseNames,addNames);
    addIndex=[classificationAdd.index]+[classificationAdd.index>0] *baseNameNum;
    classification.index(classificationAdd.index>0)=0;
    if sum(size(addIndex)==size(baseClassification.index))==2
    classification.index=baseClassification.index+addIndex;
    else
        classification.index=baseClassification.index+addIndex';
    end
else
    keyboard
    %WRITE THIS
end
end
    