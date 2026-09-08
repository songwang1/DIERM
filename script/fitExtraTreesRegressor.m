function model = fitExtraTreesRegressor(X, Y, options)
%FITEXTRATREESREGRESSOR Fit a multi-output extremely randomized-tree model.
%   This is a self-contained MATLAB implementation of the essential
%   Extra-Trees algorithm used by the Python workflow: every tree uses the
%   full training sample, a random subset of predictors is considered at
%   each node, and one uniformly random threshold is tested per predictor.
%
%   Required inputs:
%     X       - N-by-P finite predictor matrix.
%     Y       - N-by-K finite response matrix.
%     options - structure; supported fields are NumTrees (600),
%               MinLeafSize (4), MaxFeatures (0.8), and Seed (42).

arguments
    X double
    Y double
    options.NumTrees (1,1) double = 600
    options.MinLeafSize (1,1) double = 4
    options.MaxFeatures (1,1) double = 0.8
    options.Seed (1,1) double = 42
    options.JointOutputs (1,1) logical = false
end
if any(~isfinite(X),'all') || any(~isfinite(Y),'all')
    error('X and Y must be finite. Apply training-derived imputation first.');
end
[~,p] = size(X);
k = size(Y,2);
mtry = max(1,min(p,ceil(options.MaxFeatures*p)));
model.options = options;
model.numberOfPredictors = p;
model.numberOfResponses = k;
model.jointOutputs = options.JointOutputs;
if options.JointOutputs
    model.trees=cell(1,options.NumTrees);
    for treeNumber=1:options.NumTrees
        stream=RandStream('mt19937ar','Seed',options.Seed+treeNumber);
        model.trees{treeNumber}=growNode((1:size(X,1))',0,stream,X,Y,mtry,options.MinLeafSize);
    end
    return
end
model.trees = cell(k,options.NumTrees);
for response = 1:k
    for treeNumber = 1:options.NumTrees
        stream = RandStream('mt19937ar','Seed', ...
            options.Seed + 100000*response + treeNumber);
        model.trees{response,treeNumber} = growNode((1:size(X,1))',0,stream,X,Y(:,response),mtry,options.MinLeafSize);
    end
end
end

function node = growNode(index, depth, stream, X, y, mtry, minLeaf)
% Recursively grow a single extremely randomized regression tree.
node.value = mean(y(index,:),1);
node.feature = int32(0); node.threshold = NaN;
node.left = []; node.right = [];
if numel(index) < 2*minLeaf || depth >= 128 || sum(var(y(index,:),1,1)) <= eps
    return
end
p = size(X,2);
features = randperm(stream,p,mtry);
bestLoss = inf; bestFeature = 0; bestThreshold = NaN; bestLeft = []; bestRight = [];
for feature = features
    values = X(index,feature);
    lo = min(values); hi = max(values);
    if ~(hi > lo), continue, end
    threshold = lo + rand(stream)*(hi-lo);
    isLeft = values <= threshold;
    if sum(isLeft) < minLeaf || sum(~isLeft) < minLeaf, continue, end
    leftIndex = index(isLeft); rightIndex = index(~isLeft);
    leftResidual = y(leftIndex,:)-mean(y(leftIndex,:),1);
    rightResidual = y(rightIndex,:)-mean(y(rightIndex,:),1);
    loss = sum(leftResidual.^2,'all')+sum(rightResidual.^2,'all');
    if loss < bestLoss
        bestLoss=loss; bestFeature=feature; bestThreshold=threshold;
        bestLeft=leftIndex; bestRight=rightIndex;
    end
end
parentResidual = y(index,:)-mean(y(index,:),1);
if bestFeature == 0 || bestLoss >= sum(parentResidual.^2,'all')-eps, return, end
node.feature = int32(bestFeature); node.threshold = bestThreshold;
node.left = growNode(bestLeft,depth+1,stream,X,y,mtry,minLeaf);
node.right = growNode(bestRight,depth+1,stream,X,y,mtry,minLeaf);
end
