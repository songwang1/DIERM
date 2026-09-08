function fold = makeSiteGroupedFolds(siteID, numberOfFolds)
%MAKESITEGROUPEDFOLDS Reproduce unshuffled scikit-learn GroupKFold logic.
%   Complete sites are assigned greedily to the currently smallest fold,
%   starting with the sites that contain the most observations. This keeps
%   all years from a site together while balancing fold sizes.

if nargin < 2, numberOfFolds = 5; end
siteID = string(siteID(:));
[sites,~,siteIndex]=unique(siteID,'sorted');
assert(numel(sites)>=numberOfFolds, ...
    'The number of unique sites must be at least the number of folds.');
counts=accumarray(siteIndex,1,[numel(sites),1]);
% NumPy's descending argsort reverses the sorted group index for ties.
[~,order]=sortrows([-counts,-(1:numel(sites))']);
foldLoad=zeros(numberOfFolds,1); siteFold=zeros(numel(sites),1);
for k=1:numel(order)
    [~,targetFold]=min(foldLoad);
    siteFold(order(k))=targetFold;
    foldLoad(targetFold)=foldLoad(targetFold)+counts(order(k));
end
fold=siteFold(siteIndex);
end
