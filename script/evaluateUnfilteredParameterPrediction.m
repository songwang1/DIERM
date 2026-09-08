function results = evaluateUnfilteredParameterPrediction(dailyFile,fitsFile,outputDirectory,options)
%EVALUATEUNFILTEREDPARAMETERPREDICTION Run grouped-site parameter diagnostics.
%   This sensitivity workflow evaluates legacy a/b/c parameterizations when
%   those columns are present in FITSFILE. It is not required for the final
%   beta0/beta1/beta2 global projection.

arguments
    dailyFile (1,1) string
    fitsFile (1,1) string
    outputDirectory (1,1) string
    options.NumTrees (1,1) double = 500
    options.Seed (1,1) double = 42
end
if ~isfolder(outputDirectory), mkdir(outputDirectory); end
daily=readtable(dailyFile,'TextType','string'); fits=readtable(fitsFile,'TextType','string');
predictors=buildAnnualPredictors(daily,true);
frame=innerjoin(fits,predictors,'Keys',{'site_id','year'});
candidateTargets={'a','b_topt_c','c','beta0','beta1','beta2'};
targets=intersect(candidateTargets,frame.Properties.VariableNames,'stable');
% Restrict the feature set to independently constructed environmental
% predictors. Fit diagnostics and alternative response parameters must not
% leak into the prediction matrix.
numericNames={'latitude','longitude','ta_mean','ta_min','ta_max','ta_std', ...
    'vpd_mean','vpd_min','vpd_max','vpd_std','precip_sum','precip_mean', ...
    'precip_max','gpp_mean','gpp_min','gpp_max','gpp_std','ta_p05','ta_p95','wet_days'};
categoricalNames=intersect({'igbp','data_hub'},frame.Properties.VariableNames,'stable');
fold=makeSiteGroupedFolds(frame.site_id,5);
metricRows={}; predictionRows={};
for j=1:numel(targets)
    target=targets{j}; y=double(frame.(target)); eligible=isfinite(y); prediction=nan(height(frame),1);
    for f=1:5
        train=eligible&fold~=f; test=eligible&fold==f;
        prep=fitPreprocessor(frame(train,:),numericNames,categoricalNames);
        Xtrain=applyPreprocessor(frame(train,:),prep); Xtest=applyPreprocessor(frame(test,:),prep);
        model=fitExtraTreesRegressor(Xtrain,y(train),NumTrees=options.NumTrees, ...
            MinLeafSize=4,MaxFeatures=0.8,Seed=options.Seed);
        prediction(test)=predictExtraTreesRegressor(model,Xtest);
    end
    use=eligible&isfinite(prediction); e=y(use)-prediction(use);
    r2=1-sum(e.^2)/sum((y(use)-mean(y(use))).^2);
    metricRows{end+1}=table(string(target),sum(use),numel(unique(frame.site_id(use))),r2, ...
        sqrt(mean(e.^2)),mean(abs(e)),'VariableNames',{'target','site_years','sites','grouped_oof_r2','rmse','mae'}); %#ok<AGROW>
    predictionRows{end+1}=table(frame.site_id(use),frame.year(use),repmat(string(target),sum(use),1), ...
        y(use),prediction(use),'VariableNames',{'site_id','year','target','observed','predicted'}); %#ok<AGROW>
end
metrics=vertcat(metricRows{:}); predictions=vertcat(predictionRows{:});
writetable(metrics,fullfile(outputDirectory,'grouped_site_parameter_prediction.csv'));
writetable(predictions,fullfile(outputDirectory,'grouped_site_parameter_oof_predictions.csv'));
results=struct('metrics',metrics,'predictions',predictions);
save(fullfile(outputDirectory,'parameter_prediction_results.mat'),'results');
end
