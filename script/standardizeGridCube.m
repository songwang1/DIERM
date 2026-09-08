function cube = standardizeGridCube(raw,nTime,nLat,nLon)
%STANDARDIZEGRIDCUBE Return data in time-by-latitude-by-longitude order.
%   NetCDF dimension order differs among products and MATLAB preserves file
%   order. This helper identifies the required permutation from dimensions.
s=size(raw); s(end+1:3)=1;
target=[nTime nLat nLon];
perms3=perms(1:3); match=[];
for k=1:size(perms3,1)
    if isequal(s(perms3(k,:)),target), match=perms3(k,:); break, end
end
if isempty(match)
    error('Cannot map array size [%s] to [%s].',num2str(s),num2str(target));
end
cube=permute(raw,match);
end
