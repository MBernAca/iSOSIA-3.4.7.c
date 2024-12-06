function [ELA_array, maxglacierExtent,area_above_ELA] = get_ELATime(fnr_array,interpFactor,Massb_threshold,IceThick_threshold)

    SPM = SPMload;
    ny=SPM.mesh.ny;
    nx = SPM.mesh.nx;
	dx = SPM.mesh.dx;
	dy = SPM.mesh.dy;
	cell_area = dx*dy;
    include = SPM.data.include;
    nfiles = max(fnr_array);
    
    % Get ELA along the stream lines for each time step
    ELA_array = zeros(1,nfiles);
    std_array = zeros(1,nfiles);
	area_above_ELA = zeros(1,nfiles);
    maxglacierExtent = 9999;
    maxglacierExtentgot = zeros(1,nfiles);

    for i = 1:nfiles
        fnr = i;
        loaddata
        index_ELA = find(Ta(include==1) < Massb_threshold & Ta(include==1) > -Massb_threshold);
        bed_include = bed(include==1); ice_include = ice(include==1);
        % Find min elevation glacier extent
        % We apply threshold on ice thickness to focus on larger glaciers
        index_extent = find(got(include==1) > 0.3);
        if index_extent
            maxglacierExtent = min(maxglacierExtent, min(bed_include(index_extent)+ice_include(index_extent)));
        else
            maxglacierExtent = max(bed_include(:)+ice_include(:));
        end

        % get GOT
        index_extent = find(got(include==1) > 2.5e5);
        if index_extent
            maxglacierExtentgot(i) = min(bed_include(index_extent)+ice_include(index_extent));
        else
            maxglacierExtentgot(i) = max(bed_include(:)+ice_include(:));
        end

        ELA_array(i) = Ta(1,1)/SPM.mprop.lrate;
		index_above_ELA = find(bed_include >= ELA_array(i));
		area_above_ELA(i) = length(index_above_ELA) .*cell_area;
    end


end