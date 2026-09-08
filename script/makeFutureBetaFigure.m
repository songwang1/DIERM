function statistics = makeFutureBetaFigure(summaryDirectory,outputBase)
%MAKEFUTUREBETAFIGURE Create the three-panel future HOT-day/ER figure.
[parent,~,~]=fileparts(outputBase); if strlength(parent)>0 && ~isfolder(parent), mkdir(parent); end
scenarios=["ssp126" "ssp245" "ssp370" "ssp585"];
labels=["SSP1-2.6" "SSP2-4.5" "SSP3-7.0" "SSP5-8.5"];
colors=[35 137 233;108 182 255;242 142 43;216 59 93]/255;
figure('Color','w','Position',[100 100 1050 480]);
layout=tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
axA=nexttile(1); hold(axA,'on'); axB=nexttile(3); hold(axB,'on'); axC=nexttile(2,[2 1]); hold(axC,'on');
statistics=struct;
for s=1:numel(scenarios)
    data=readtable(fullfile(summaryDirectory,"global_beta_"+scenarios(s)+"_annual_summary.csv"));
    plot(axA,data.year,data.hot_days_mean,'LineWidth',1.45,'Color',colors(s,:),'DisplayName',labels(s));
    plot(axB,data.year,data.global_er_gt_c_yr,'LineWidth',1.45,'Color',colors(s,:));
    scatter(axC,data.hot_days_mean,data.global_er_gt_c_yr,25,'MarkerEdgeColor',colors(s,:), ...
        'MarkerFaceColor','w','DisplayName',labels(s));
    x=data.hot_days_mean; y=data.global_er_gt_c_yr; fit=hacRegression(y,[ones(height(data),1),x,x.^2],3);
    xx=linspace(min(x),max(x),180)'; design=[ones(size(xx)),xx,xx.^2]; yy=design*fit.beta;
    predictionSE=sqrt(max(0,sum((design*fit.covariance).*design,2)));
    critical=tinv(0.975,max(1,fit.n-size(design,2)));
    lower=yy-critical.*predictionSE; upper=yy+critical.*predictionSE;
    fill(axC,[xx;flipud(xx)],[lower;flipud(upper)],colors(s,:), ...
        'FaceAlpha',0.12,'EdgeColor','none','HandleVisibility','off');
    plot(axC,xx,yy,'LineWidth',1.8,'Color',colors(s,:),'HandleVisibility','off');
    field=char(scenarios(s)); statistics.(field)=struct( ...
        'hot_2015_2024_mean',mean(x(data.year<=2024)), ...
        'hot_2091_2100_mean',mean(x(data.year>=2091)), ...
        'er_2015_2024_mean',mean(y(data.year<=2024)), ...
        'er_2091_2100_mean',mean(y(data.year>=2091)), ...
        'quadratic_hac_p_linear',fit.p(2),'quadratic_hac_p_quadratic',fit.p(3));
end
xlabel(axB,'Year'); ylabel(axA,'HOT days (days yr^{-1})'); ylabel(axB,'ER (Gt C yr^{-1})');
xlabel(axC,'HOT days (days yr^{-1})'); ylabel(axC,'ER (Gt C yr^{-1})');
legend(axA,'Location','northwest','Box','off'); legend(axC,'Location','best','Box','off');
title(axA,'A','HorizontalAlignment','left'); title(axB,'B','HorizontalAlignment','left'); title(axC,'C','HorizontalAlignment','left');
set(findall(gcf,'Type','axes'),'FontName','Arial','FontSize',11,'Box','on');
exportgraphics(layout,outputBase+".png",'Resolution',600); exportgraphics(layout,outputBase+".pdf",'ContentType','vector');
fid=fopen(outputBase+"_statistics.json",'w'); fprintf(fid,'%s',jsonencode(statistics,'PrettyPrint',true)); fclose(fid);
end
