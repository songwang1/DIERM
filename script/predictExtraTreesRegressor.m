function prediction = predictExtraTreesRegressor(model, X)
%PREDICTEXTRATREESREGRESSOR Predict responses from a fitted Extra-Trees model.

n = size(X,1); k = model.numberOfResponses;
if isfield(model,'jointOutputs') && model.jointOutputs
    prediction=zeros(n,k);
    for treeNumber=1:numel(model.trees)
        prediction=prediction+predictTree(model.trees{treeNumber},X,k);
    end
    prediction=prediction/numel(model.trees);
    return
end
prediction = zeros(n,k);
for response = 1:k
    accumulator = zeros(n,1);
    for treeNumber = 1:size(model.trees,2)
        accumulator = accumulator + predictTree(model.trees{response,treeNumber},X,1);
    end
    prediction(:,response) = accumulator/size(model.trees,2);
end
end

function value = predictTree(tree,X,numberOfResponses)
n = size(X,1); value = nan(n,numberOfResponses);
nodeStack={tree}; rowStack={(1:n)'};
while ~isempty(nodeStack)
    node=nodeStack{end}; rows=rowStack{end};
    nodeStack(end)=[]; rowStack(end)=[];
    if isempty(rows), continue, end
    if node.feature==0
        value(rows,:)=repmat(node.value,numel(rows),1);
    else
        isLeft=X(rows,node.feature)<=node.threshold;
        nodeStack{end+1}=node.left; rowStack{end+1}=rows(isLeft); %#ok<AGROW>
        nodeStack{end+1}=node.right; rowStack{end+1}=rows(~isLeft); %#ok<AGROW>
    end
end
end
