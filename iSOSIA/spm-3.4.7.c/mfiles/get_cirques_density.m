function [cirques_density,abrasion_mean,fluvial_mean,indices_cirques,cirques_density_initial] = get_cirques_density(curv,relief,bed,cc,rc,abrasion,fluvial,abrasion_threshold,xc,yc,fig)

    bin_curvature = linspace(-0.0020,0.0020,50);
    %get channels
    DEM = GRIDobj(xc,yc,bed);
    curv_DEM = DEM;
    curv_DEM.Z = flipud(curv);
    grid_size = DEM.size;
    curv = flipud(curv);
    relief = flipud(relief);
    bed = flipud(bed);
    Meanbed = localtopography(DEM,10000,'type','mean');
    xc = flipud(xc); yc = flipud(yc);

    % Flow accumulation
    FD = FLOWobj(DEM,"preprocess",'carve');
    A = flowacc(FD);

    % Stream
    S = STREAMobj(FD,A>0.001e2);
    S = removeshortstreams(S,1500);
    S = modify(S,'maxdsdistance',1000);

    % get channel heads
    xy = streampoi(S,'channelheads');
    channel_heads = streampoi(S,'channelheads','ix');
    indexes_ch = channel_heads;
    indices_channelhead = find(bed(indexes_ch) > Meanbed.Z(indexes_ch));
    indexes_ch = indexes_ch(bed(indexes_ch) > Meanbed.Z(indexes_ch));
    
	Meanbed = Meanbed.Z;
    if fig
        figure
        subplot(1,2,1)
        imageschs(DEM,dilate(sqrt(A),ones(2)),'colormap','flowcolor','caxis',[0 10])
        subplot(1,2,2)
        imageschs(DEM,[],'colormap',[1 1 1]);
        hold on
        plotc(S,DEM);
        plot(xy(indices_channelhead,1),xy(indices_channelhead,2),'s')
        colormap(jet);
        hold off
    end
    % For each channel heads, check if this is a glacial cirque
    % Count cirques above cc
    abrasion_DEM = abrasion;
    indices_cirques = [];
    indices_initialTopo = [];
    nn = 1;
    nb_cell_abrasion = 0;
    abrasion_mean = 0;
    for k=1:length(indexes_ch)
        [row,col] = ind2sub(grid_size,indexes_ch(k));
        row_min = max(row - nn, 1);
        row_max = min(row + nn, grid_size(1));
        col_min = max(col - nn, 1);
        col_max = min(col + nn, grid_size(2));
        num_neighbors = (row_max - row_min + 1) * (col_max - col_min + 1) - 1; % minus central point
        neighbor_indices = zeros(num_neighbors+1, 1);
        
        idx = 1;
        % get neighbors
        for i = row_min:row_max
            for j = col_min:col_max
                if i == row && j == col
                    continue; % Ignore central point
                end
                
                % Convert coordinates of neighbors to linear index
                neighbor_indices(idx) = sub2ind(grid_size, i, j);
                idx = idx + 1;
            end
        end
        % 
        neighbor_indices(end) = indexes_ch(k);
        % find cirque
        bcurv = curv(neighbor_indices); relief_temp = relief(neighbor_indices); bed_temp = bed(neighbor_indices); abrasion_temp = abrasion(neighbor_indices); MeanBed_temp = Meanbed(neighbor_indices);
        fluvial_temp = fluvial(neighbor_indices);
        index = find(bcurv>= cc & abrasion_temp>=abrasion_threshold & relief_temp > rc);
        index2 = find(bcurv>= cc & relief_temp > rc);
        if index % save channel head
            indices_cirques = [indices_cirques,indexes_ch(k)];
            abrasion_mean = abrasion_mean + sum(abrasion_temp(index));
            nb_cell_abrasion = nb_cell_abrasion + length(index);
        end
        if index2
            indices_initialTopo = [indices_initialTopo,indexes_ch(k)];
        end
    end
    fluvial_mean = mean(fluvial_temp(index));
    abrasion_mean = abrasion_mean / nb_cell_abrasion;
    cirques(indices_cirques) = 1;
    cirques_density = length(indices_cirques);
    if isempty(indices_initialTopo)
        cirques_density_initial = 0;
    else
        cirques_density_initial = length(indices_initialTopo);
    end

    if fig
        figure 
        subplot(2,2,1)
        hold on
        imageschs(DEM,[],'colormap',[1 1 1]);view([90 90])
        plot(xc(indices_cirques),yc(indices_cirques),'o','MarkerFaceColor','b','MarkerSize',3,'MarkerEdgeColor','none')
        hold off
        subplot(2,2,2)
        hold on
        imageschs(DEM,[],'colormap',[1 1 1]);view([90 90])
        scatter(xy(indices_channelhead,1),xy(indices_channelhead,2),5,curv(indexes_ch),'fill');caxis([0 cc]);colormap(jet);
        hold off
        subplot(2,2,3)
        imageschs(DEM,flipud(abrasion_DEM),'caxis',[0 abrasion_threshold]);view([90 90])
    end

end