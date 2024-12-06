function [cdata_M,mindata,maxdata] = get_color(data,colormapp,clim,c,ndim,bed,ice,mindata,maxdata)

caxis(clim); cmap = colormapp;
mintopo = min(data(:)); %mintopo = log(1);
maxtopo = max(data(:)); 
if c,
    cdata = 0.99*(data-mintopo)/(maxtopo-mintopo);
    I = find(cdata(:) > 0.99); cdata(I) = 0.99;
    icemap = ice>5;
    cdata(icemap==1) = 1;
else
    indx = data > 0;
    I = find(abs(data(:)) > 0.000001*(maxdata-mindata));
    cdata =  0.95*(data-mindata)/(maxdata-mindata+1e-16); 
    I = cdata(:) < 0; cdata(I) = 0;
    I = cdata(:) > 1; cdata(I) = 1;
    cdata = cdata +4;
    cdata(indx == 0) = 0.99*(data(indx==0)-mindata)/(maxdata -mindata);
    I = find(cdata(indx==0) > 0.99); cdata(I) = 0.99;
end
%make Cdata_M
%Elevation
xi = linspace(clim(1),clim(2),length(cmap(:,1)));
rc = interp1(xi,cmap(:,1),cdata);
gc = interp1(xi,cmap(:,2),cdata);
bc = interp1(xi,cmap(:,3),cdata);
if c
    rc(icemap) = 0.8;
    bc(icemap)=0.8;
    gc(icemap) =0.8;
end
if ndim,
    cdata_M(:,:,1) = rc;
    cdata_M(:,:,2) = gc;
    cdata_M(:,:,3) = bc;
else
    cdata_M(:,1) = rc;
    cdata_M(:,2) = gc;
    cdata_M(:,3) = bc;
end
    
   
end