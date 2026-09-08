function statistics = makeFigure2CD(summaryFile,h5File,outputBase,options)
%MAKEFIGURE2CD Create the historical DIERM panels corresponding to Fig. 2C-D.
%   Annual confidence bands are 95% intervals from a 5-degree spatial-block
%   bootstrap. The panel-D test uses a 5-degree cluster-robust covariance.

arguments
    summaryFile (1,1) string
    h5File (1,1) string
    outputBase (1,1) string
    options.BootstrapDraws (1,1) double = 2000
    options.Seed (1,1) double = 20260826
end
annual=readtable(summaryFile); year=annual.year;
[parent,~,~]=fileparts(outputBase); if strlength(parent)>0 && ~isfolder(parent), mkdir(parent); end
names=["hot_days" "topt_c" "mat_c"];
for j=1:numel(names)
    y=annual.(names(j)+"_mean"); fits.(names(j))=hacRegression(y,[ones(numel(year),1),year],2);
end
ci=spatialBlockBootstrap(h5File,names,options.BootstrapDraws,options.Seed);
red=[199 54 47]/255; blue=[21 155 211]/255; black=[22 22 22]/255;
figure('Color','w','Position',[100 100 1040 435]); layout=tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
ax=nexttile; hold(ax,'on'); yyaxis(ax,'left');
fillBand(ax,year,ci.hot_days,0.85*[1 1 1]);
plot(ax,year,annual.hot_days_mean,'Color',black,'LineWidth',1.55,'DisplayName','Global HOT days');
plot(ax,year,[ones(numel(year),1),year]*fits.hot_days.beta,'--','Color',black,'LineWidth',1.2, ...
    'DisplayName',sprintf('HOT days trend: %.2f d yr^{-1}, %s',fits.hot_days.beta(2),pText(fits.hot_days.p(2))));
ylabel(ax,'HOT days (days yr^{-1})'); ylim(ax,[0 160]);
yyaxis(ax,'right'); fillBand(ax,year,ci.mat_c,[0.95 0.81 0.80]);
plot(ax,year,annual.mat_c_mean,'Color',red,'LineWidth',1.15,'DisplayName','Global MAT');
plot(ax,year,[ones(numel(year),1),year]*fits.mat_c.beta,'--','Color',red,'LineWidth',1.2, ...
    'DisplayName',sprintf('MAT trend: %.3f deg C yr^{-1}, %s',fits.mat_c.beta(2),pText(fits.mat_c.p(2))));
fillBand(ax,year,ci.topt_c,[0.80 0.92 0.97]);
plot(ax,year,annual.topt_c_mean,'Color',blue,'LineWidth',1.15,'DisplayName','Global T_{opt}');
plot(ax,year,[ones(numel(year),1),year]*fits.topt_c.beta,'--','Color',blue,'LineWidth',1.2, ...
    'DisplayName',sprintf('T_{opt} trend: %.3f deg C yr^{-1}, %s',fits.topt_c.beta(2),pText(fits.topt_c.p(2))));
ylabel(ax,'Temperature (degrees C)'); xlabel(ax,'Year'); title(ax,'C','HorizontalAlignment','left');
legend(ax,'Location','best','Box','off','FontSize',8);

ax=nexttile; hold(ax,'on'); years=double(h5read(h5File,'/year')); yi=find(years==2020,1);
lat=double(h5read(h5File,'/lat')); lon=double(h5read(h5File,'/lon')); area=gridCellArea(lat,lon);
x=double(h5read(h5File,'/mat_c',[1 1 yi],[numel(lat) numel(lon) 1]));
y=double(h5read(h5File,'/hot_days',[1 1 yi],[numel(lat) numel(lon) 1]));
[row,col]=ndgrid(1:numel(lat),1:numel(lon)); cluster=floor((row-1)/10)*ceil(numel(lon)/10)+floor((col-1)/10)+1;
ok=isfinite(x)&isfinite(y); x=x(ok); y=y(ok); w=area(ok); cluster=cluster(ok);
if numel(x)>1500000
    rng(options.Seed); take=randperm(numel(x),1500000); x=x(take); y=y(take); w=w(take); cluster=cluster(take);
end
[count,xEdge,yEdge]=histcounts2(x,y,105); imageData=log1p(count'); imageData=imageData/max(imageData,[],'all');
xc=(xEdge(1:end-1)+xEdge(2:end))/2; yc=(yEdge(1:end-1)+yEdge(2:end))/2;
imagesc(ax,xc,yc,imageData,'AlphaData',count'>0); axis(ax,'xy'); colormap(ax,'parula'); cb=colorbar(ax); ylabel(cb,'Normalized pixel density');
scaled=w/mean(w,'omitnan'); curve=@(b,z) b(3)+b(1).*exp(b(2).*z);
objective=@(b) sqrt(scaled).*(curve(b,x)-y);
b=lsqnonlin(objective,[12 0.08 5],[0 0 0],[1000 1 365],optimoptions('lsqnonlin','Display','off'));
xx=linspace(max(-10,min(x)),min(32,max(x)),300); plot(ax,xx,curve(b,xx),'Color',red,'LineWidth',2);
robust=clusterRobustRegression(y,[ones(numel(x),1),x,x.^2],cluster,scaled);
R=[0 1 0;0 0 1]; statistic=(R*robust.beta)'/(R*robust.covariance*R')*(R*robust.beta);
p=1-chi2cdf(statistic,2); text(ax,0.03,0.94,pText(p),'Units','normalized','VerticalAlignment','top');
xlabel(ax,'Mean annual air temperature (degrees C)'); ylabel(ax,'HOT days (days yr^{-1})');
xlim(ax,[-10 32]); ylim(ax,[0 365]); title(ax,'D','HorizontalAlignment','left');
set(findall(gcf,'Type','axes'),'FontName','Arial','FontSize',11,'Box','on');
exportgraphics(layout,outputBase+".png",'Resolution',600); exportgraphics(layout,outputBase+".pdf",'ContentType','vector');
statistics=struct('trends_hac',fits,'panel_d_year',2020,'panel_d_points',numel(x), ...
    'panel_d_exponential_parameters',b,'panel_d_cluster_robust_joint_p',p);
fid=fopen(outputBase+"_statistics.json",'w'); fprintf(fid,'%s',jsonencode(statistics,'PrettyPrint',true)); fclose(fid);
end

function ci=spatialBlockBootstrap(h5File,names,draws,seed)
lat=double(h5read(h5File,'/lat')); lon=double(h5read(h5File,'/lon')); years=h5read(h5File,'/year');
area=gridCellArea(lat,lon); nLat=numel(lat); nLon=numel(lon); nYear=numel(years);
assert(mod(nLat,10)==0&&mod(nLon,10)==0,'The 5-degree bootstrap requires a 0.5-degree grid.');
valid=isfinite(h5read(h5File,'/topt_c')); blockWeight=zeros(nLat/10*nLon/10,nYear);
sums=struct; for name=names, sums.(name)=zeros(size(blockWeight)); end
for yi=1:nYear
    w=area.*valid(:,:,yi); blockWeight(:,yi)=blockSum(w);
    for name=names
        value=double(h5read(h5File,"/"+name,[1 1 yi],[nLat nLon 1]));
        current=sums.(name); current(:,yi)=blockSum(value.*w); sums.(name)=current;
    end
end
keep=any(blockWeight>0,2); blockWeight=blockWeight(keep,:);
for name=names, current=sums.(name); sums.(name)=current(keep,:); end
rng(seed); nBlock=size(blockWeight,1);
for name=names
    boot=zeros(draws,nYear);
    for d=1:draws
        take=randi(nBlock,nBlock,1); current=sums.(name);
        boot(d,:)=sum(current(take,:),1)./sum(blockWeight(take,:),1);
    end
    ci.(name)=prctile(boot,[2.5 97.5],1);
end
end
function s=blockSum(x)
[nLat,nLon]=size(x); block=reshape(x,10,nLat/10,10,nLon/10); s=reshape(sum(sum(block,1),3),[],1);
end
function fillBand(ax,x,interval,color)
fill(ax,[x;flipud(x)],[interval(1,:)';flipud(interval(2,:)')],color,'EdgeColor','none','FaceAlpha',0.55,'HandleVisibility','off');
end
function textValue=pText(p)
if p<0.001, textValue='P < 0.001'; else, textValue=sprintf('P = %.3f',p); end
end
