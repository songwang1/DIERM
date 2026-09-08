function prep = fitPreprocessor(frame, numericNames, categoricalNames)
%FITPREPROCESSOR Learn training-only imputations and category levels.
%   The fitted structure must be reused without modification on validation
%   data. This mirrors a leakage-safe preprocessing pipeline.

prep.numericNames = numericNames;
prep.categoricalNames = categoricalNames;
prep.numericMedian = nan(1,numel(numericNames));
for j = 1:numel(numericNames)
    prep.numericMedian(j) = median(double(frame.(numericNames{j})),'omitnan');
    if ~isfinite(prep.numericMedian(j)), prep.numericMedian(j)=0; end
end
prep.levels = cell(1,numel(categoricalNames));
prep.categoryMode = strings(1,numel(categoricalNames));
for j = 1:numel(categoricalNames)
    x = string(frame.(categoricalNames{j}));
    x(ismissing(x)) = "missing";
    prep.levels{j} = unique(x,'stable');
    prep.categoryMode(j) = string(mode(categorical(x)));
end
end
