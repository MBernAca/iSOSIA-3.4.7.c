% Swath profile
function [d,z,x,y,mind,maxd,meand,xd] = get_swath_profile(xi, yi, n, fnr, limy, fig)
    SPM = SPMload;
    nx = SPM.mesh.nx; ny = SPM.mesh.ny;
    loaddata
    include = SPM.data.include;
    DEM = GRIDobj(reshape(xc(include==1),236,108)./1e3,reshape(yc(include==1),236,108)./1e3,reshape(bed(include==1),236,108)./1e3);
    yl = reshape(yc(include==1),236,108);
    % Profile
    [d,z,x,y] = demprofile(DEM,n,xi,yi);
    
    meand = mean(DEM.Z');
    mind = min(DEM.Z');
    maxd = max(DEM.Z');
    xd = yl(:,1)./1e3;
    % plot profile
    if fig
        figure
        plot(d,z)
        xlabel('Distance along profile (km)')
        ylabel('Elevation (km)')
        ylim([limy(1) limy(2)])
    end

end