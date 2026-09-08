function bundle = trainGlobalBetaModel(dailyFile, fitsFile, outputDirectory, options)
%TRAINGLOBALBETAMODEL Train the final globally applicable DIERM beta model.
%   Cross-validation is grouped by site. The final model is fitted only
%   after OOF predictions have been generated for model evaluation.

arguments
    dailyFile (1,1) string
    fitsFile (1,1) string
    outputDirectory (1,1) string
    options.NumTrees (1,1) double = 600
    options.Seed (1,1) double = 20260826
end
if ~isfolder(outputDirectory), mkdir(outputDirectory); end
daily=readtable(dailyFile,'TextType','string');
fits=readtable(fitsFile,'TextType','string');
predictors=buildAnnualPredictors(daily,false);
frame=innerjoin(fits,predictors,'Keys',{'site_id','year'});
features={'latitude','longitude','ta_mean','ta_min','ta_max','ta_std','ta_p05','ta_p95', ...
    'vpd_mean','vpd_min','vpd_max','vpd_std','precip_sum','precip_mean','gpp_mean', ...
    'precip_month_min','precip_month_max','precip_month_std','gpp_month_min', ...
    'gpp_month_max','gpp_month_std'};
targets={'beta0','beta1','beta2'};
Y=frame{:,targets};
valid=all(isfinite(Y),2); frame=frame(valid,:); Y=Y(valid,:);
fold=makeSiteGroupedFolds(frame.site_id,5);
oof=nan(size(Y));
for f=1:5
    train=fold~=f; test=fold==f;
    prep=fitPreprocessor(frame(train,:),features,{});
    Xtrain=applyPreprocessor(frame(train,:),prep);
    Xtest=applyPreprocessor(frame(test,:),prep);
    model=fitExtraTreesRegressor(Xtrain,Y(train,:),NumTrees=options.NumTrees, ...
        MinLeafSize=4,MaxFeatures=0.8,Seed=42+f);
    oof(test,:)=predictExtraTreesRegressor(model,Xtest);
end
oof(:,3)=min(oof(:,3),-1e-10);
rows=cell(3,1);
for j=1:3
    e=Y(:,j)-oof(:,j); r2=1-sum(e.^2)/sum((Y(:,j)-mean(Y(:,j))).^2);
    rows{j}=table("extra_trees",string(targets{j}),height(frame),numel(unique(frame.site_id)), ...
        r2,sqrt(mean(e.^2)),mean(abs(e)),'VariableNames', ...
        {'model','target','n','sites','oof_r2','rmse','mae'});
end
metrics=vertcat(rows{:});
writetable(metrics,fullfile(outputDirectory,'global_beta_parameter_oof_metrics.csv'));
oofTable=table(frame.site_id,frame.year,Y(:,1),Y(:,2),Y(:,3),oof(:,1),oof(:,2),oof(:,3), ...
    'VariableNames',{'site_id','year','observed_beta0','observed_beta1','observed_beta2', ...
    'predicted_beta0','predicted_beta1','predicted_beta2'});
writetable(oofTable,fullfile(outputDirectory,'global_beta_parameter_oof_predictions.csv'));

prep=fitPreprocessor(frame,features,{});
X=applyPreprocessor(frame,prep);
model=fitExtraTreesRegressor(X,Y,NumTrees=options.NumTrees,MinLeafSize=4, ...
    MaxFeatures=0.8,Seed=options.Seed);
bundle=struct('model',model,'preprocessor',prep,'features',{features}, ...
    'targets',{targets},'tref_c',12,'equation', ...
    'log(ER)=beta0+beta1*(Ta-12)+beta2*(Ta-12)^2; beta2<0');
rangeTable=table(string(features(:)),'VariableNames',{'feature'});
rangeTable.p01=nan(numel(features),1); rangeTable.median=nan(numel(features),1); rangeTable.p99=nan(numel(features),1);
for j=1:numel(features)
    value=double(frame.(features{j})); value=value(isfinite(value));
    rangeTable{j,{'p01','median','p99'}}=prctile(value,[1 50 99]);
end
writetable(rangeTable,fullfile(outputDirectory,'training_feature_ranges.csv'));
bundle.feature_p01=rangeTable.p01';
bundle.feature_p99=rangeTable.p99';
save(fullfile(outputDirectory,'global_beta_model.mat'),'bundle','-v7.3');
end
