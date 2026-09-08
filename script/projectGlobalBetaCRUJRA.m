function summary = projectGlobalBetaCRUJRA(config)
%PROJECTGLOBALBETACRUJRA Run the 1991-2020 CRU-JRA/CEDAR DIERM projection.
%   CONFIG requires ModelFile, CruDirectory, CedarCache, OutputH5,
%   OutputCSV, StartYear, and EndYear. The output HDF5 fields use
%   latitude-by-longitude-by-year order.

arguments
    config struct
end
required={'ModelFile','CruDirectory','CedarCache','OutputH5','OutputCSV','StartYear','EndYear'};
for k=1:numel(required), assert(isfield(config,required{k}),'Missing config.%s',required{k}); end
loaded=load(config.ModelFile,'bundle'); bundle=loaded.bundle;
[parentH5,~,~]=fileparts(config.OutputH5); if strlength(parentH5)>0 && ~isfolder(parentH5), mkdir(parentH5); end
[parentCSV,~,~]=fileparts(config.OutputCSV); if strlength(parentCSV)>0 && ~isfolder(parentCSV), mkdir(parentCSV); end
years=config.StartYear:config.EndYear;
sample=fullfile(config.CruDirectory,sprintf('crujra.v2.5.5d.tmp.%d.365d.noc.nc',years(1)));
latitude=double(ncread(sample,'lat')); longitude=double(ncread(sample,'lon'));
nLat=numel(latitude); nLon=numel(longitude); nYear=numel(years);
area=gridCellArea(latitude,longitude);
if isfile(config.OutputH5), error('Output already exists: %s',config.OutputH5); end
h5create(config.OutputH5,'/year',[nYear 1],'Datatype','int32'); h5write(config.OutputH5,'/year',int32(years(:)));
h5create(config.OutputH5,'/lat',[nLat 1]); h5write(config.OutputH5,'/lat',latitude(:));
h5create(config.OutputH5,'/lon',[nLon 1]); h5write(config.OutputH5,'/lon',longitude(:));
variables={'beta0','beta1','beta2','topt_c','hot_days','mat_c','annual_er_g_c_m2'};
for k=1:numel(variables)
    h5create(config.OutputH5,"/"+variables{k},[nLat nLon nYear], ...
        'Datatype','single','ChunkSize',[min(45,nLat) min(90,nLon) 1],'Deflate',4,'FillValue',single(NaN));
end
gppCache=matfile(config.CedarCache);
cacheYears=double(gppCache.years);
rows=cell(nYear,1);
for yi=1:nYear
    year=years(yi);
    ta=readCRU(config.CruDirectory,'tmp','tmp',year,'mean',nLat,nLon)-273.15;
    q=readCRU(config.CruDirectory,'spfh','spfh',year,'mean',nLat,nLon);
    pressure=readCRU(config.CruDirectory,'pres','pres',year,'mean',nLat,nLon);
    actualVP=q.*pressure./(0.622+0.378.*q);
    saturationVP=611.*exp(17.27.*ta./(ta+237.3));
    vpd=max((saturationVP-actualVP)./100,0);
    precip=readCRU(config.CruDirectory,'pre','pre',year,'sum',nLat,nLon);
    cacheIndex=find(cacheYears==year,1);
    assert(~isempty(cacheIndex),'CEDAR cache does not contain year %d.',year);
    gpp=squeeze(double(gppCache.gpp_g_c_m2_day(cacheIndex,:,:,:)))/(86400*12e-6);
    forcing=struct('taDaily',ta,'vpdDaily',vpd,'precipDaily',precip,'gppMonthly',gpp);
    result=projectDiermYear(forcing,bundle,latitude,longitude);
    for k=1:numel(variables)
        h5write(config.OutputH5,"/"+variables{k},single(result.(variables{k})),[1 1 yi],[nLat nLon 1]);
    end
    weights=area; weights(~result.valid)=0;
    row=struct('year',year);
    for name=["beta0" "beta1" "beta2" "topt_c" "hot_days" "mat_c"]
        one=weightedSummary(result.(name),weights,name); row=mergeStruct(row,one);
    end
    row.global_er_gt_c_yr=sum(result.annual_er_g_c_m2(:).*weights(:),'omitnan')/1e15;
    if isfield(result,'outside_training_range')
        row.outside_any_training_p01_p99_area_fraction= ...
            sum(area(result.outside_training_range),'omitnan')/sum(area(result.valid),'omitnan');
    end
    row.valid_land_area_million_km2=sum(weights(:),'omitnan')/1e12;
    rows{yi}=row; summary=struct2table(vertcat(rows{1:yi}));
    writetable(summary,config.OutputCSV);
    fprintf('%d: MAT=%.2f, Topt=%.2f, HOT=%.1f, ER=%.2f Gt C yr^-1\n', ...
        year,row.mat_c_mean,row.topt_c_mean,row.hot_days_mean,row.global_er_gt_c_yr);
end
end

function daily=readCRU(directory,stem,variable,year,reducer,nLat,nLon)
file=fullfile(directory,sprintf('crujra.v2.5.5d.%s.%d.365d.noc.nc',stem,year));
sixHourly=readNetCDFTimeSlice(file,variable,1,1460,nLat,nLon);
sixHourly=reshape(sixHourly,4,365,nLat,nLon);
if reducer=="mean", daily=squeeze(mean(sixHourly,1,'omitnan'));
else, daily=squeeze(sum(sixHourly,1,'omitnan')); end
end
function a=mergeStruct(a,b)
names=fieldnames(b); for j=1:numel(names), a.(names{j})=b.(names{j}); end
end
