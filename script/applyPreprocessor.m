function X = applyPreprocessor(frame, prep)
%APPLYPREPROCESSOR Apply fixed imputation and one-hot encoding rules.

n = height(frame);
X = nan(n,numel(prep.numericNames));
for j = 1:numel(prep.numericNames)
    x = double(frame.(prep.numericNames{j}));
    x(~isfinite(x)) = prep.numericMedian(j);
    X(:,j) = x;
end
for j = 1:numel(prep.categoricalNames)
    x = string(frame.(prep.categoricalNames{j}));
    x(ismissing(x)) = prep.categoryMode(j);
    levels = prep.levels{j};
    block = zeros(n,numel(levels));
    for k = 1:numel(levels)
        block(:,k) = x == levels(k);
    end
    X = [X block]; %#ok<AGROW>
end
end
