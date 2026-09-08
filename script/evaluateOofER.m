function [metrics,annual] = evaluateOofER(daily,oofTable,trefC)
%EVALUATEOOFER Reconstruct ER from OOF beta estimates and score predictions.
%   Metrics are reported for pooled daily records, retained temperature-bin
%   means, site-year means, within-site site-year anomalies, and site means.

params=oofTable(:,{'site_id','year','predicted_beta0','predicted_beta1','predicted_beta2'});
params.Properties.VariableNames(3:5)={'beta0','beta1','beta2'};
if ~ismember('year',daily.Properties.VariableNames)
    stamp=daily.TIMESTAMP; if ~isnumeric(stamp), stamp=str2double(string(stamp)); end
    daily.year=floor(double(stamp)/10000);
end
daily.site_id=string(daily.site_id);
ta=daily.TA_F; if ~isnumeric(ta), ta=str2double(string(ta)); end
er=daily.RECO_DT_VUT_REF; if ~isnumeric(er), er=str2double(string(er)); end
daily.TA_F=double(ta); daily.RECO_DT_VUT_REF=double(er);
valid=isfinite(daily.TA_F)&isfinite(daily.RECO_DT_VUT_REF)& ...
    daily.TA_F>=-55&daily.TA_F<=60&daily.RECO_DT_VUT_REF>0;
x=innerjoin(daily(valid,{'site_id','year','TA_F','RECO_DT_VUT_REF'}),params,'Keys',{'site_id','year'});
x.predicted_er=diermER(x{:,{'beta0','beta1','beta2'}},x.TA_F,trefC);
rows=cell(5,1); rows{1}=metricRow(x.RECO_DT_VUT_REF,x.predicted_er,"daily_pooled",NaN);

x.temperature_bin=floor(x.TA_F);
[g,site,year,bin]=findgroups(x.site_id,x.year,x.temperature_bin);
observed=splitapply(@mean,x.RECO_DT_VUT_REF,g); predicted=splitapply(@mean,x.predicted_er,g); days=splitapply(@numel,x.TA_F,g);
binned=table(site,year,bin,observed,predicted,days); binned=binned(binned.days>=3,:);
rows{2}=metricRow(binned.observed,binned.predicted,"temperature_bin_pooled",NaN);

[g,site,year]=findgroups(x.site_id,x.year);
observed_er=splitapply(@mean,x.RECO_DT_VUT_REF,g); predicted_er=splitapply(@mean,x.predicted_er,g); valid_days=splitapply(@numel,x.TA_F,g);
annual=table(site,year,observed_er,predicted_er,valid_days, ...
    'VariableNames',{'site_id','year','observed_er','predicted_er','valid_days'});
rows{3}=metricRow(annual.observed_er,annual.predicted_er,"site_year_mean",NaN);
[siteGroup,~]=findgroups(annual.site_id);
siteObservedMean=splitapply(@mean,annual.observed_er,siteGroup);
sitePredictedMean=splitapply(@mean,annual.predicted_er,siteGroup);
annual.obs_anom=annual.observed_er-siteObservedMean(siteGroup);
annual.pred_anom=annual.predicted_er-sitePredictedMean(siteGroup);
siteCount=splitapply(@numel,annual.year,siteGroup); multi=siteCount(siteGroup)>=2;
rows{4}=metricRow(annual.obs_anom(multi),annual.pred_anom(multi), ...
    "site_year_within_site_anomaly",numel(unique(annual.site_id(multi))));
siteObserved=siteObservedMean; sitePredicted=sitePredictedMean;
rows{5}=metricRow(siteObserved,sitePredicted,"site_mean_between_sites",numel(siteObserved));
metrics=vertcat(rows{:});
end

function row=metricRow(y,pred,level,sites)
ok=isfinite(y)&isfinite(pred); y=y(ok); pred=pred(ok); e=y-pred;
r2=1-sum(e.^2)/sum((y-mean(y)).^2);
row=table("extra_trees_joint",level,numel(y),sites,r2,sqrt(mean(e.^2)),mean(abs(e)), ...
    'VariableNames',{'scenario','level','n','sites','r2','rmse','mae'});
end
