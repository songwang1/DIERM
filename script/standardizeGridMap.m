function map = standardizeGridMap(raw,nLat,nLon)
%STANDARDIZEGRIDMAP Return a two-dimensional map as latitude-by-longitude.
raw=squeeze(double(raw));
if isequal(size(raw),[nLat nLon])
    map=raw;
elseif isequal(size(raw),[nLon nLat])
    map=raw';
else
    error('Unexpected spatial dimensions [%s].',num2str(size(raw)));
end
end
