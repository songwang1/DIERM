function preaggregateCedarGPP(cedarDirectory,outputFile,startYear,endYear)
%PREAGGREGATECEDARGPP Aggregate monthly CEDAR-GPP to the CRU 0.5-degree grid.
%   CEDARDIRECTORY must contain the extracted monthly NetCDF files. The
%   function writes a version-7.3 MAT file with gpp_g_c_m2_day stored as
%   year-by-month-by-latitude-by-longitude.

arguments
    cedarDirectory (1,1) string
    outputFile (1,1) string
    startYear (1,1) double = 1991
    endYear (1,1) double = 2020
end
files=dir(fullfile(cedarDirectory,'**','*.nc'));
years=(startYear:endYear)'; months=(1:12)';
if isfile(outputFile), error('Output already exists: %s',outputFile); end
[parent,~,~]=fileparts(outputFile);
if strlength(parent)>0 && ~isfolder(parent), mkdir(parent); end
save(outputFile,'years','months','-v7.3');
target=matfile(outputFile,'Writable',true);
target.gpp_g_c_m2_day(endYear-startYear+1,12,360,720)=single(NaN);
for y=startYear:endYear
    for m=1:12
        key=sprintf('%04d%02d',y,m);
        match=find(contains(string({files.name}),key),1);
        if isempty(match), error('Missing CEDAR-GPP month %s.',key); end
        file=fullfile(files(match).folder,files(match).name);
        value=double(squeeze(ncread(file,'GPP_mean')));
        info=ncinfo(file); variableNames=string({info.Variables.Name});
        if any(variableNames=="lat"), lat=double(ncread(file,'lat'));
        else, lat=double(ncread(file,'y')); end
        if any(variableNames=="lon"), lon=double(ncread(file,'lon'));
        else, lon=double(ncread(file,'x')); end
        if isequal(size(value),[7200 3600]), value=value'; end
        if ~isequal(size(value),[3600 7200])
            error('Unexpected CEDAR grid in %s.',file);
        end
        value(value==-9999)=NaN; value=value*0.01;
        value(value<0|value>100)=NaN;
        % Match the ascending CRU latitude grid and the [-180, 180)
        % longitude convention before block aggregation.
        if lat(1)>lat(end), value=flipud(value); lat=flipud(lat); end %#ok<NASGU>
        wrappedLon=mod(lon+180,360)-180; [~,order]=sort(wrappedLon);
        value=value(:,order);
        blocks=reshape(value,10,360,10,720);
        coarse=squeeze(mean(mean(blocks,1,'omitnan'),3,'omitnan'));
        target.gpp_g_c_m2_day(y-startYear+1,m,:,:)=single(coarse);
    end
    fprintf('Aggregated CEDAR-GPP %d\n',y);
end
end
