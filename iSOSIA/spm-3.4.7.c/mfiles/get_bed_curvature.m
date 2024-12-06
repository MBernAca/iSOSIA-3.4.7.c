function [bcurv] = get_bed_curvature(bed,dx,dy)

ny = size(bed,1);
nx = size(bed,2);
bcurv= zeros(size(bed));

for i=2:ny-1
    for j=2:nx-1
        bcurv(i,j) = (bed(i,j+1) - 2.0*bed(i,j) + bed(i,j-1)) / (dx*dx) + (bed(i+1,j) - 2.0*bed(i,j) + bed(i-1,j)) / (dy*dy);
    end
end

end