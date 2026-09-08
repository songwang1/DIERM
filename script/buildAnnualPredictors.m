function out = buildAnnualPredictors(daily, includeCategories)
%BUILDANNUALPREDICTORS Aggregate daily FLUXNET data to site-year predictors.
%   OUT = BUILDANNUALPREDICTORS(DAILY) reproduces the annual and monthly
%   summaries used by the DIERM beta-parameter model. Missing-value codes at
%   or below -9990 are converted to NaN. Air temperature is restricted to
%   [-55, 60] degrees C, and non-positive GPP is treated as missing.

if nargin < 2, includeCategories = false; end
daily.site_id = string(daily.site_id);
stamp = double(daily.TIMESTAMP);
daily.year = floor(stamp/10000);
daily.month = mod(floor(stamp/100),100);
names = {'TA_F','VPD_F','P_F','GPP_DT_VUT_REF'};
for j = 1:numel(names)
    v = double(daily.(names{j}));
    v(v <= -9990) = NaN;
    daily.(names{j}) = v;
end
daily.TA_F(daily.TA_F < -55 | daily.TA_F > 60) = NaN;
daily.GPP_DT_VUT_REF(daily.GPP_DT_VUT_REF <= 0) = NaN;

[group,site,year] = findgroups(daily.site_id,daily.year);
n = max(group);
latitude = splitapply(@firstFinite,double(daily.location_lat),group);
longitude = splitapply(@firstFinite,double(daily.location_long),group);
ta_mean = splitapply(@nanMean,daily.TA_F,group);
ta_min = splitapply(@nanMin,daily.TA_F,group);
ta_max = splitapply(@nanMax,daily.TA_F,group);
ta_std = splitapply(@nanStd,daily.TA_F,group);
ta_p05 = splitapply(@(x) nanPct(x,5),daily.TA_F,group);
ta_p95 = splitapply(@(x) nanPct(x,95),daily.TA_F,group);
vpd_mean = splitapply(@nanMean,daily.VPD_F,group);
vpd_min = splitapply(@nanMin,daily.VPD_F,group);
vpd_max = splitapply(@nanMax,daily.VPD_F,group);
vpd_std = splitapply(@nanStd,daily.VPD_F,group);
precip_sum = splitapply(@(x) sum(x,'omitnan'),daily.P_F,group);
precip_mean = splitapply(@nanMean,daily.P_F,group);
precip_max = splitapply(@nanMax,daily.P_F,group);
gpp_mean = splitapply(@nanMean,daily.GPP_DT_VUT_REF,group);
gpp_min = splitapply(@nanMin,daily.GPP_DT_VUT_REF,group);
gpp_max = splitapply(@nanMax,daily.GPP_DT_VUT_REF,group);
gpp_std = splitapply(@nanStd,daily.GPP_DT_VUT_REF,group);
wet_days = splitapply(@(x) sum(x>0 & isfinite(x)),daily.P_F,group);

% Monthly means are summarized across the 12 calendar months. Missing
% months remain missing and do not contribute to the summary statistics.
precipMonthly = nan(n,12); gppMonthly = nan(n,12);
for m = 1:12
    use = daily.month == m;
    if any(use)
        precipMonthly(:,m) = accumarray(group(use),daily.P_F(use),[n 1],@(x) mean(x,'omitnan'),NaN);
        gppMonthly(:,m) = accumarray(group(use),daily.GPP_DT_VUT_REF(use),[n 1],@(x) mean(x,'omitnan'),NaN);
    end
end
precip_month_min = min(precipMonthly,[],2,'omitnan');
precip_month_max = max(precipMonthly,[],2,'omitnan');
precip_month_std = std(precipMonthly,0,2,'omitnan');
gpp_month_min = min(gppMonthly,[],2,'omitnan');
gpp_month_max = max(gppMonthly,[],2,'omitnan');
gpp_month_std = std(gppMonthly,0,2,'omitnan');

out = table(site,year,latitude,longitude,ta_mean,ta_min,ta_max,ta_std, ...
    ta_p05,ta_p95,vpd_mean,vpd_min,vpd_max,vpd_std,precip_sum,precip_mean,precip_max, ...
    gpp_mean,gpp_min,gpp_max,gpp_std,wet_days,precip_month_min,precip_month_max,precip_month_std, ...
    gpp_month_min,gpp_month_max,gpp_month_std, ...
    'VariableNames',{'site_id','year','latitude','longitude','ta_mean','ta_min', ...
    'ta_max','ta_std','ta_p05','ta_p95','vpd_mean','vpd_min','vpd_max', ...
    'vpd_std','precip_sum','precip_mean','precip_max','gpp_mean','gpp_min','gpp_max', ...
    'gpp_std','wet_days','precip_month_min', ...
    'precip_month_max','precip_month_std','gpp_month_min','gpp_month_max', ...
    'gpp_month_std'});
if includeCategories
    out.igbp = splitapply(@(x) firstString(x),string(daily.igbp),group);
    out.data_hub = splitapply(@(x) firstString(x),string(daily.data_hub),group);
end
end

function y = firstFinite(x)
x = x(isfinite(x)); if isempty(x), y = NaN; else, y = x(1); end
end
function y = nanMean(x), y = mean(x,'omitnan'); end
function y = nanMin(x), if all(isnan(x)), y=NaN; else, y=min(x,[],'omitnan'); end, end
function y = nanMax(x), if all(isnan(x)), y=NaN; else, y=max(x,[],'omitnan'); end, end
function y = nanStd(x), y = std(x,0,'omitnan'); end
function y = nanPct(x,p)
x=x(isfinite(x)); if isempty(x), y=NaN; else, y=prctile(x,p); end
end
function y = firstString(x)
x=x(~ismissing(x)); if isempty(x), y="missing"; else, y=x(1); end
end
