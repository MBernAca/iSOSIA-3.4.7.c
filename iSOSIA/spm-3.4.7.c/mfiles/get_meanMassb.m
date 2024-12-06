function [meanMb, stdMb] = get_meanMassb(fnr)

SPM = SPMload;
ny=SPM.mesh.ny;
nx = SPM.mesh.nx;
include = SPM.data.include;
loaddata;

index = find(massb(SPM.data.include==1) < 0.0001 & massb(SPM.data.include==1) > -0.0001 & ice(SPM.data.include==1) > 0);
bed_include = bed(SPM.data.include==1); massb_include = massb(SPM.data.include==1); ice_include = ice(SPM.data.include==1);
meanELA = mean(bed_include(index)+ice_include(index))
stdELA = std(bed_include(index)+ice_include(index))


end