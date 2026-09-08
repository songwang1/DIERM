function output = summarizeGlobalBetaSSP(summaryDirectory,outputFile)
%SUMMARIZEGLOBALBETASSP Summarize baseline, future change, and HAC trends.
scenarios=["ssp126" "ssp245" "ssp370" "ssp585"];
sources={'beta0_mean','beta1_mean','beta2_mean','topt_c_mean','hot_days_mean','global_er_gt_c_yr'};
labels={'beta0','beta1','beta2','topt_c','hot_days','global_er_gt_c_yr'};
rows=cell(numel(scenarios),1);
for s=1:numel(scenarios)
    file=fullfile(summaryDirectory,"global_beta_"+scenarios(s)+"_annual_summary.csv");
    data=readtable(file);
    assert(height(data)==86&&min(data.year)==2015&&max(data.year)==2100,'Incomplete SSP series.');
    row=struct('scenario',scenarios(s)); early=data.year>=2015&data.year<=2024; late=data.year>=2091&data.year<=2100;
    for j=1:numel(sources)
        value=data.(sources{j}); trend=hacRegression(value,[ones(height(data),1),data.year-mean(data.year)],3);
        row.(labels{j}+"_2015_2024")=mean(value(early),'omitnan');
        row.(labels{j}+"_2091_2100")=mean(value(late),'omitnan');
        row.(labels{j}+"_change")=row.(labels{j}+"_2091_2100")-row.(labels{j}+"_2015_2024");
        row.(labels{j}+"_trend_per_year")=trend.beta(2); row.(labels{j}+"_trend_hac_p")=trend.p(2);
    end
    hot=data.hot_days_mean; er=data.global_er_gt_c_yr;
    fit=hacRegression(er,[ones(height(data),1),hot,hot.^2],3);
    row.er_hot_quadratic_p=fit.p(3); row.er_hot_peak_days=-fit.beta(2)/(2*fit.beta(3)); rows{s}=row;
end
output=struct2table(vertcat(rows{:}));
[parent,~,~]=fileparts(outputFile); if strlength(parent)>0 && ~isfolder(parent), mkdir(parent); end
writetable(output,outputFile);
end
