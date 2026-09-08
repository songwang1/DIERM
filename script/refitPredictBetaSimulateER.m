function results = refitPredictBetaSimulateER(dailyFile, outputDirectory, options)
%REFITPREDICTBETASIMULATEER Fit DIERM curves and run site-grouped OOF tests.
%   RESULTS = REFITPREDICTBETASIMULATEER(DAILYFILE, OUTPUTDIRECTORY) is the
%   MATLAB entry point corresponding to the current site-year DIERM fitting
%   workflow. DAILYFILE is the merged daily FLUXNET CSV file.

arguments
    dailyFile (1,1) string
    outputDirectory (1,1) string
    options.Tref (1,1) double = 12
    options.NumTrees (1,1) double = 500
    options.Seed (1,1) double = 42
end
if ~isfolder(outputDirectory), mkdir(outputDirectory); end
daily = readtable(dailyFile,'TextType','string');
required = {'site_id','TIMESTAMP','TA_F','VPD_F','P_F','GPP_DT_VUT_REF', ...
    'RECO_DT_VUT_REF','location_lat','location_long','igbp','data_hub'};
assert(all(ismember(required,daily.Properties.VariableNames)), ...
    'The daily input file is missing one or more required variables.');

stamp = numericColumn(daily.TIMESTAMP);
daily.year = floor(stamp/10000);
daily.site_id = string(daily.site_id);
daily.TA_F = numericColumn(daily.TA_F);
daily.RECO_DT_VUT_REF = numericColumn(daily.RECO_DT_VUT_REF);
[group,site,year] = findgroups(daily.site_id,daily.year);
fitRows = cell(max(group),1);
for g = 1:max(group)
    one = group==g;
    ta = daily.TA_F(one); er = daily.RECO_DT_VUT_REF(one);
    valid = isfinite(ta)&isfinite(er)&ta>=-55&ta<=60&er>0;
    if sum(valid)<300, continue, end
    ta = ta(valid); er = er(valid);
    bin = floor(ta);
    [bg,~] = findgroups(bin);
    binnedTa = splitapply(@mean,ta,bg);
    binnedER = splitapply(@mean,er,bg);
    days = splitapply(@numel,er,bg);
    keep = days>=3;
    if sum(keep)<6, continue, end
    fitted = fitSiteYearBeta(binnedTa(keep),binnedER(keep),options.Tref);
    if isempty(fitted), continue, end
    tmin = min(ta); tmax = max(ta); vertex=fitted.implied_topt_c;
    if ~isfinite(vertex), optimumClass="unidentified_near_linear";
    elseif vertex<tmin, optimumClass="below_observed_range";
    elseif vertex>tmax, optimumClass="above_observed_range";
    else, optimumClass="within_observed_range";
    end
    prefix = table(site(g),year(g),options.Tref,sum(valid),sum(keep),tmin,tmax,mean(ta),optimumClass, ...
        'VariableNames',{'site_id','year','tref_c','valid_days','temperature_bins', ...
        'tmin_c','tmax_c','mat_c','optimum_class'});
    fitRows{g}=[prefix fitted];
end
fits = vertcat(fitRows{~cellfun(@isempty,fitRows)});
writetable(fits,fullfile(outputDirectory,'all_site_year_beta_fits.csv'));

predictors = buildAnnualPredictors(daily,true);
frame = innerjoin(fits,predictors,'Keys',{'site_id','year'});
numericNames = {'latitude','longitude','ta_mean','ta_min','ta_max','ta_std', ...
    'vpd_mean','vpd_min','vpd_max','vpd_std','precip_sum','precip_mean', ...
    'precip_max','gpp_mean','gpp_min','gpp_max','gpp_std','ta_p05','ta_p95','wet_days'};
categoricalNames = {'igbp','data_hub'};
Y = frame{:,{'beta0','beta1','beta2'}};
fold = makeSiteGroupedFolds(frame.site_id,5);
oof = nan(size(Y));
for f = 1:5
    train = fold~=f; test=fold==f;
    prep = fitPreprocessor(frame(train,:),numericNames,categoricalNames);
    Xtrain = applyPreprocessor(frame(train,:),prep);
    Xtest = applyPreprocessor(frame(test,:),prep);
    model = fitExtraTreesRegressor(Xtrain,Y(train,:), ...
        NumTrees=options.NumTrees,MinLeafSize=4,MaxFeatures=0.8, ...
        Seed=options.Seed,JointOutputs=true);
    oof(test,:) = predictExtraTreesRegressor(model,Xtest);
end
oof(:,3)=min(oof(:,3),-1e-10);
oofTable=table(frame.site_id,frame.year,Y(:,1),Y(:,2),Y(:,3), ...
    oof(:,1),oof(:,2),oof(:,3),repmat("extra_trees_joint",height(frame),1), ...
    'VariableNames',{'site_id','year','observed_beta0','observed_beta1','observed_beta2', ...
    'predicted_beta0','predicted_beta1','predicted_beta2','model'});
writetable(oofTable,fullfile(outputDirectory,'beta_parameter_oof_predictions.csv'));
target = ["beta0" "beta1" "beta2"];
metricRows=cell(3,1);
for j=1:3
    metricRows{j}=regressionMetrics(Y(:,j),oof(:,j),target(j),numel(unique(frame.site_id)));
end
metrics=vertcat(metricRows{:});
writetable(metrics,fullfile(outputDirectory,'beta_parameter_oof_metrics.csv'));
[erMetrics,annualER]=evaluateOofER(daily,oofTable,options.Tref);
writetable(erMetrics,fullfile(outputDirectory,'beta_er_oof_metrics.csv'));
writetable(annualER,fullfile(outputDirectory,'beta_er_site_year_predictions.csv'));
save(fullfile(outputDirectory,'site_year_oof_results.mat'),'fits','predictors','frame', ...
    'oofTable','metrics','erMetrics','annualER','-v7.3');
results=struct('fits',fits,'oofPredictions',oofTable,'metrics',metrics, ...
    'erMetrics',erMetrics,'annualER',annualER);
end

function x=numericColumn(x)
if ~isnumeric(x), x=str2double(string(x)); else, x=double(x); end
end
function row=regressionMetrics(y,pred,target,numberOfSites)
residual=y-pred; sst=sum((y-mean(y)).^2);
r2=1-sum(residual.^2)/sst;
row=table("extra_trees_joint",target,numel(y),numberOfSites,r2, ...
    sqrt(mean(residual.^2)),mean(abs(residual)), ...
    'VariableNames',{'model','target','n','sites','oof_r2','rmse','mae'});
end
