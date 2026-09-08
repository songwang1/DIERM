function summary = projectGlobalBetaSSP(config)
%PROJECTGLOBALBETASSP Project DIERM parameters and ER for one CMIP6 SSP.
%   CONFIG requires Scenario, ModelFile, DailyDirectory, MonthlyDirectory,
%   OutputH5, OutputCSV, StartYear, and EndYear. Expected scenarios are
%   ssp126, ssp245, ssp370, and ssp585.

arguments
    config struct
end
required={'Scenario','ModelFile','DailyDirectory','MonthlyDirectory','OutputH5', ...
    'OutputCSV','StartYear','EndYear'};
for k=1:numel(required), assert(isfield(config,required{k}),'Missing config.%s',required{k}); end
scenario=string(config.Scenario);
assert(any(scenario==["ssp126" "ssp245" "ssp370" "ssp585"]),'Unsupported SSP.');
loaded=load(config.ModelFile,'bundle'); bundle=loaded.bundle;
[parentH5,~,~]=fileparts(config.OutputH5); if strlength(parentH5)>0 && ~isfolder(parentH5), mkdir(parentH5); end
[parentCSV,~,~]=fileparts(config.OutputCSV); if strlength(parentCSV)>0 && ~isfolder(parentCSV), mkdir(parentCSV); end
areaFile=firstFile(config.DailyDirectory,"areacella_fx_CESM2-WACCM_"+scenario+"_*.nc");
landFile=firstFile(config.DailyDirectory,"sftlf_fx_CESM2-WACCM_"+scenario+"_*.nc");
latitude=double(ncread(areaFile,'lat')); longitude=double(ncread(areaFile,'lon'));
nLat=numel(latitude); nLon=numel(longitude);
area=standardizeGridMap(ncread(areaFile,'areacella'),nLat,nLon);
landFraction=standardizeGridMap(ncread(landFile,'sftlf'),nLat,nLon)/100;
cellWeight=area.*landFraction;
years=config.StartYear:config.EndYear; nYear=numel(years);
if isfile(config.OutputH5), error('Output already exists: %s',config.OutputH5); end
h5create(config.OutputH5,'/year',[nYear 1],'Datatype','int32'); h5write(config.OutputH5,'/year',int32(years(:)));
h5create(config.OutputH5,'/lat',[nLat 1]); h5write(config.OutputH5,'/lat',latitude(:));
h5create(config.OutputH5,'/lon',[nLon 1]); h5write(config.OutputH5,'/lon',longitude(:));
variables={'beta0','beta1','beta2','topt_c','hot_days','annual_er_g_c_m2'};
for k=1:numel(variables)
    h5create(config.OutputH5,"/"+variables{k},[nLat nLon nYear],'Datatype','single', ...
        'ChunkSize',[min(48,nLat) min(72,nLon) 1],'Deflate',4,'FillValue',single(NaN));
end
rows=cell(nYear,1);
for yi=1:nYear
    year=years(yi);
    tasFile=findTimeFile(config.DailyDirectory,"tas_day_CESM2-WACCM_"+scenario+"_r1i1p1f1_gn_*.nc",year);
    hursFile=findTimeFile(config.DailyDirectory,"hurs_day_CESM2-WACCM_"+scenario+"_r1i1p1f1_gn_*.nc",year);
    gppFile=findTimeFile(config.MonthlyDirectory,"gpp_Lmon_CESM2-WACCM_"+scenario+"_r1i1p1f1_gn_*.nc",year);
    prFile=findTimeFile(config.MonthlyDirectory,"pr_Amon_CESM2-WACCM_"+scenario+"_r1i1p1f1_gn_*.nc",year);
    ta=readYear(tasFile,'tas',year,365,nLat,nLon)-273.15;
    rh=readYear(hursFile,'hurs',year,365,nLat,nLon);
    vpd=max(0.611.*exp(17.27.*ta./(ta+237.3)).*(1-rh/100).*10,0);
    gpp=readYear(gppFile,'gpp',year,12,nLat,nLon)/(12e-9);
    precip=readYear(prFile,'pr',year,12,nLat,nLon)*86400;
    forcing=struct('taDaily',ta,'vpdDaily',vpd,'gppMonthly',gpp,'precipMonthly',precip);
    result=projectDiermYear(forcing,bundle,latitude,longitude);
    for k=1:numel(variables)
        h5write(config.OutputH5,"/"+variables{k},single(result.(variables{k})),[1 1 yi],[nLat nLon 1]);
    end
    weights=cellWeight; weights(~result.valid)=0;
    row=struct('scenario',scenario,'year',year);
    for name=["beta0" "beta1" "beta2" "topt_c" "hot_days"]
        row=mergeStruct(row,weightedSummary(result.(name),weights,name));
    end
    row.global_er_gt_c_yr=sum(result.annual_er_g_c_m2(:).*weights(:),'omitnan')/1e15;
    if isfield(result,'outside_training_range')
        row.land_weight_outside_any_training_p01_p99= ...
            sum(cellWeight(result.outside_training_range),'omitnan')/sum(cellWeight(result.valid),'omitnan');
        denominator=sum(cellWeight(result.valid),'omitnan');
        for featureIndex=1:numel(bundle.features)
            mask=result.outside_training_by_feature(:,:,featureIndex);
            field="outside_"+bundle.features{featureIndex}+"_p01_p99_area_fraction";
            row.(field)=sum(cellWeight(mask),'omitnan')/denominator;
        end
    end
    finiteTopt=result.valid&isfinite(result.topt_c); denominator=sum(cellWeight(finiteTopt),'omitnan');
    row.topt_below_minus20_area_fraction=sum(cellWeight(finiteTopt&result.topt_c<-20),'omitnan')/denominator;
    row.topt_above_50_area_fraction=sum(cellWeight(finiteTopt&result.topt_c>50),'omitnan')/denominator;
    row.valid_land_area_million_km2=sum(weights(:),'omitnan')/1e12;
    rows{yi}=row; summary=struct2table(vertcat(rows{1:yi})); writetable(summary,config.OutputCSV);
    fprintf('%s %d: Topt=%.2f, HOT=%.1f, ER=%.2f Gt C yr^-1\n',scenario,year, ...
        row.topt_c_mean,row.hot_days_mean,row.global_er_gt_c_yr);
end
end

function file=firstFile(directory,pattern)
files=dir(fullfile(directory,pattern)); assert(~isempty(files),'No file matched %s',pattern);
file=fullfile(files(1).folder,files(1).name);
end
function file=findTimeFile(directory,pattern,year)
files=dir(fullfile(directory,pattern)); file="";
for k=1:numel(files)
    token=regexp(files(k).name,'(\d{6,8})-(\d{6,8})','tokens','once');
    if isempty(token), continue, end
    if year>=str2double(token{1}(1:4)) && year<=str2double(token{2}(1:4))
        file=string(fullfile(files(k).folder,files(k).name)); return
    end
end
assert(strlength(file)>0,'No file covering year %d matched %s.',year,pattern);
end
function cube=readYear(file,variable,year,stepsPerYear,nLat,nLon)
token=regexp(file,'(\d{6,8})-(\d{6,8})','tokens','once'); startYear=str2double(token{1}(1:4));
first=(year-startYear)*stepsPerYear+1;
cube=readNetCDFTimeSlice(file,variable,first,stepsPerYear,nLat,nLon);
end
function a=mergeStruct(a,b)
names=fieldnames(b); for j=1:numel(names), a.(names{j})=b.(names{j}); end
end
