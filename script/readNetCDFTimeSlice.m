function cube = readNetCDFTimeSlice(file,variable,firstStep,numberOfSteps,nLat,nLon)
%READNETCDFTIMESLICE Read a time interval and standardize its grid order.
info=ncinfo(file,variable); names=string({info.Dimensions.Name});
timeDimension=find(contains(lower(names),'time'),1);
if isempty(timeDimension)
    sizes=[info.Dimensions.Length];
    timeDimension=find(sizes>=firstStep+numberOfSteps-1,1,'last');
end
start=ones(1,numel(info.Dimensions)); count=[info.Dimensions.Length];
start(timeDimension)=firstStep; count(timeDimension)=numberOfSteps;
raw=double(ncread(file,variable,start,count));
cube=standardizeGridCube(raw,numberOfSteps,nLat,nLon);
cube(abs(cube)>1e10)=NaN;
end
