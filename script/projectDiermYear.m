function result = projectDiermYear(forcing,bundle,latitude,longitude)
%PROJECTDIERMYEAR Predict beta, Topt, HOT days, and ER for one gridded year.
%   Required forcing fields use time-by-latitude-by-longitude order:
%     taDaily, vpdDaily, gppMonthly
%   Supply either precipDaily or precipMonthly. Monthly rates must represent
%   daily means for their respective months.

monthDays=[31 28 31 30 31 30 31 31 30 31 30 31];
ta=double(forcing.taDaily); vpd=double(forcing.vpdDaily);
gpp=double(forcing.gppMonthly);
[nDay,nLat,nLon]=size(ta);
assert(nDay==365,'The current workflow uses a 365-day calendar.');
if isfield(forcing,'precipDaily')
    precipDaily=double(forcing.precipDaily);
    edge=[0 cumsum(monthDays)]; precipMonthly=nan(12,nLat,nLon);
    for m=1:12
        precipMonthly(m,:,:)=mean(precipDaily(edge(m)+1:edge(m+1),:,:),1,'omitnan');
    end
    precipSum=squeeze(sum(precipDaily,1,'omitnan'));
    precipMean=squeeze(mean(precipDaily,1,'omitnan'));
else
    precipMonthly=double(forcing.precipMonthly);
    precipSum=squeeze(sum(precipMonthly.*reshape(monthDays,12,1,1),1,'omitnan'));
    precipMean=precipSum/365;
end
[lonGrid,latGrid]=meshgrid(longitude,latitude);
q=prctile(ta,[5 95],1);
arrays=struct;
arrays.latitude=latGrid; arrays.longitude=lonGrid;
arrays.ta_mean=squeeze(mean(ta,1,'omitnan')); arrays.ta_min=squeeze(min(ta,[],1,'omitnan'));
arrays.ta_max=squeeze(max(ta,[],1,'omitnan')); arrays.ta_std=squeeze(std(ta,0,1,'omitnan'));
arrays.ta_p05=squeeze(q(1,:,:)); arrays.ta_p95=squeeze(q(2,:,:));
arrays.vpd_mean=squeeze(mean(vpd,1,'omitnan')); arrays.vpd_min=squeeze(min(vpd,[],1,'omitnan'));
arrays.vpd_max=squeeze(max(vpd,[],1,'omitnan')); arrays.vpd_std=squeeze(std(vpd,0,1,'omitnan'));
arrays.precip_sum=precipSum; arrays.precip_mean=precipMean;
arrays.gpp_mean=squeeze(sum(gpp.*reshape(monthDays,12,1,1),1,'omitnan'))/365;
arrays.precip_month_min=squeeze(min(precipMonthly,[],1,'omitnan'));
arrays.precip_month_max=squeeze(max(precipMonthly,[],1,'omitnan'));
arrays.precip_month_std=squeeze(std(precipMonthly,0,1,'omitnan'));
arrays.gpp_month_min=squeeze(min(gpp,[],1,'omitnan'));
arrays.gpp_month_max=squeeze(max(gpp,[],1,'omitnan'));
arrays.gpp_month_std=squeeze(std(gpp,0,1,'omitnan'));

frame=table;
rawX=nan(nLat*nLon,numel(bundle.features));
for j=1:numel(bundle.features)
    name=bundle.features{j}; value=arrays.(name);
    frame.(name)=value(:); rawX(:,j)=value(:);
end
X=applyPreprocessor(frame,bundle.preprocessor);
% Match the production projection: training-derived imputation is part of
% the fitted model, but gridded cells with any missing raw predictor are not
% projected and must remain missing in the output.
valid=all(isfinite(rawX),2)&arrays.gpp_mean(:)>1e-6;
beta=nan(nLat*nLon,3);
batch=20000; index=find(valid);
for first=1:batch:numel(index)
    take=index(first:min(first+batch-1,numel(index)));
    beta(take,:)=predictExtraTreesRegressor(bundle.model,X(take,:));
end
beta(:,3)=min(beta(:,3),-1e-10);
annual=diermAnnualMetrics(reshape(ta,nDay,[]),beta,bundle.tref_c);
result=struct;
result.beta0=reshape(beta(:,1),nLat,nLon); result.beta1=reshape(beta(:,2),nLat,nLon);
result.beta2=reshape(beta(:,3),nLat,nLon); result.topt_c=reshape(annual.topt_c,nLat,nLon);
result.hot_days=reshape(annual.hot_days,nLat,nLon);
result.annual_er_g_c_m2=reshape(annual.annual_er_g_c_m2,nLat,nLon);
result.mat_c=arrays.ta_mean; result.valid=reshape(valid,nLat,nLon);
if isfield(bundle,'feature_p01') && isfield(bundle,'feature_p99')
    outsideEach=(rawX<bundle.feature_p01 | rawX>bundle.feature_p99) & valid;
    outside=any(outsideEach,2);
    result.outside_training_range=reshape(outside,nLat,nLon);
    result.outside_training_by_feature=reshape(outsideEach,nLat,nLon,[]);
end
end
