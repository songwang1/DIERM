function area = gridCellArea(latitude,longitude)
%GRIDCELLAREA Calculate spherical grid-cell area in square metres.
radius=6371000;
dlat=median(diff(latitude)); dlon=deg2rad(median(diff(longitude)));
lower=deg2rad(max(-90,latitude(:)-dlat/2));
upper=deg2rad(min(90,latitude(:)+dlat/2));
row=radius^2*dlon.*(sin(upper)-sin(lower));
area=repmat(row,1,numel(longitude));
end
