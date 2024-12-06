
% This code extract the outputs from 9 models : Colder (T0 - 2°C),
% Reference, Warmer (T0 + 2°C), Climatic cycles amplitude x2, No climatic
% cycles, 41 kyr cycles only, 100 kyr cycles only, Precipitation rate x
% 1.5, and Precipitation rate x 0.75.

% The results are saved in .mat files to be used by the matlab codes
% 'plot_figures_part1.m' and 'Merge_results.m'

%%%%%% Extract information %%%%%%

% Plotting paramaters 
path_to_models = '../Models/';
path_to_results = '../Models_output';
path_to_figures = '../Figures/';
path_to_save_analyses = '../Model_analyses/';
path_to_codes = '../Codes/';
Models = {'Temp_-2C','Reference','Temp_+2C','Ampx2','41_kyr_cycles','100_kyr_cycles',...
    'No_cycles','Prec_x1-5','Prec_x075'};
number_of_Models = length(Models);
Colors_Models = {'c';[0,0,0];'b';'r';[0,255,0]./255;[0,153,0]./255;
    [255,153,153]./255;[1,0.5,1];[255,0,255]./255};
Legend_Models = {'Colder (-2°C)';'Reference';'Warmer (+2°C)';'Amp. x2';'41 kyrs';'100 kyrs';'No cycles';'Prx1.5';...
    'Prx0.75'};

% Load model characteritics
load([path_to_models,Models{1},'/input/SPM.mat'])
number_of_output_files = 146;
nx = SPM.mesh.nx; Model_length_x = SPM.mesh.L;
ny = SPM.mesh.ny; Model_length_y = SPM.mesh.H;
dx = Model_length_x / nx;
dy = Model_length_y / ny;
cell_area = dx * dy;
x = dx/2:dx:dx*128;
y = dy/2:dy:dy*256;
fid = fopen([path_to_models,'Reference/output/output1.dat']); %file handle
xc = fread(fid,[ny,nx],'double'); %cell x coordinate (m)
yc = fread(fid,[ny,nx],'double'); %cell y coordinate (m)
fclose(fid);

% Parameters for model analyses
Bed_slope_threshold = 0.177; % [gradient]
Relief_threshold = 200;
curvature_threshold_for_cirques = 0.001;
abrasion_threshold_for_cirques = 5;
Relief_threshold_for_cirques = 500;
minimum_elevation = -1000;
maximum_elevation = 2000;
number_of_elevation_bin = 50;
elevation_bins = linspace(minimum_elevation,maximum_elevation,number_of_elevation_bin+1);
region_of_interest = SPM.data.include;
Initial_topography = SPM.data.bed;
number_of_cells = 108*236; % number of elements in include zone
time1 = [1:85].*0.0205;
time2 = [89:149].*0.02;
Model_time = [time1,time2];


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Main model analyses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section extract main results and perform some analyses from the
% model results. This section can be run once. All results for plotting
% are save in .mat file at the end of this section.


% Get bed gradients for initial topography
Initial_bed_slope = get_bed_curvature(SPM.data.bed,dx,dy);

Initial_LRS_area_per_elevation_bin = zeros(size(elevation_bins));
Initial_hypsometry_area_per_elevation_bin = zeros(length(elevation_bins),1); %for elevation

Initial_DEM = GRIDobj(xc,yc,Initial_topography);
Initial_bed = flipud(Initial_DEM.Z);
Initial_bed = Initial_bed(region_of_interest ==1);

[~,~,initial_elev_bin] = histcounts(Initial_bed,elevation_bins);
Initial_relief = localtopography(Initial_DEM,500);
Initial_relief = flipud(Initial_relief.Z);
Initial_relief = Initial_relief(region_of_interest==1);

%Compute Mean local topography (10 km) to apply as a filter for LRS
Initial_Mean_DEM = localtopography(Initial_DEM,10000,'type','mean');
Initial_Mean_DEM = flipud(Initial_Mean_DEM.Z);
Initial_Mean_DEM = Initial_Mean_DEM(region_of_interest==1);

% Calculate on initial topography
for i=1:length(elevation_bins)
	%Initial topo
	index = find(initial_elev_bin==i);
	% count slopes < 10° for the bin elevation for initial topo
	slope_found = Initial_bed_slope(index);
	relief_found = Initial_relief(index);
	index_LRS = find(slope_found < Bed_slope_threshold & relief_found < Relief_threshold & Initial_bed(index) > Initial_Mean_DEM(index));
	Initial_LRS_area_per_elevation_bin(i) = length(index_LRS) * cell_area * 1e-6; % Area of LRS (km²)
	Initial_hypsometry_area_per_elevation_bin(i) = length(index) * cell_area * 1e-6; % Area of surfaces (km²)
end

% Get data value for all models
Total_LRS_area = zeros(number_of_Models,1);
ELA_array = zeros(number_of_output_files,number_of_Models);
std_ELA_array = zeros(number_of_output_files,number_of_Models);
maxglacierExtent = zeros(number_of_output_files,number_of_Models);
LRS_area_per_elevation_bin = zeros(length(elevation_bins),number_of_Models);
Hypsometry_area_per_elevation_bin = zeros(length(elevation_bins),number_of_Models); %for elevation
Erosion_per_elevation_bin = zeros(length(elevation_bins),number_of_Models);

Erosion_models = zeros(number_of_Models,1);
Abrasion_models = zeros(number_of_Models,1);
Fluvial_models = zeros(number_of_Models,1);
Hillslope_models = zeros(number_of_Models,1);
isostasy_models = zeros(number_of_Models,1);
Mean_GOT = zeros(number_of_Models,1);
std_GOT = zeros(number_of_Models,1);
Area_above_MeanElev = zeros(number_of_Models,1);
Erosion_Spatial = zeros(number_of_cells,number_of_Models);
Fluvial_Spatial = zeros(number_of_cells,number_of_Models);
Glacial_Spatial = zeros(number_of_cells,number_of_Models);
GOT_Spatial = zeros(number_of_cells,number_of_Models);
Mean_bed_Spatial = zeros(number_of_cells,number_of_Models);
bed_Spatial = zeros(number_of_cells,number_of_Models);
isostasy_evolution_models = zeros(number_of_output_files,number_of_Models);
maxTopo_evolution_models = zeros(number_of_output_files,number_of_Models);
meanTopo_evolution_models = zeros(number_of_output_files,number_of_Models);
minTopo_evolution_models = zeros(number_of_output_files,number_of_Models);
area_above_ELA = zeros(number_of_output_files,number_of_Models);
ice_vol_mean_interglacials = zeros(number_of_output_files,number_of_Models);
Mean_water_discharge_evolution = zeros(number_of_output_files,number_of_Models);
cirques_density_time = zeros(number_of_output_files,number_of_Models);

% Loop through models
cd(path_to_models)
for i = 1:number_of_Models
    disp(['Doing model ',num2str(i),'/',num2str(number_of_Models)])
    fnr = number_of_output_files;
    loaddata_models % Load model outcomes from iSOSIA (internal iSOSIA function)

    if abs(abs((xc(1,2)-xc(1,1)))-abs(yc(2,1)-yc(1,1)))>1e-9
	    x = dx/2:dx:dx*nx;
	    y = dx/2:dx:dx*ny;
	    [xc,yc] = meshgrid(x,y);
    end

    DEM = GRIDobj(xc, yc, bed);
    bed = bed(region_of_interest == 1);
    bed_Spatial(:,i) = bed';
    slopes = bslope(region_of_interest==1);
    Abrasion_models(i) = sum(abrasion(region_of_interest==1));
    Fluvial_models(i) = sum(fluvial(region_of_interest==1));
    Hillslope_models(i) = sum(hillslope(region_of_interest==1));
    isostasy_models(i) = mean(isostasy(region_of_interest==1));
    Erosion_models(i) = sum(abrasion(region_of_interest == 1) + hillslope(region_of_interest==1) + fluvial(region_of_interest==1) + periglacial(region_of_interest==1));
    Erosion_Spatial(:,i) = (abrasion(region_of_interest == 1) + hillslope(region_of_interest==1) + fluvial(region_of_interest==1) + periglacial(region_of_interest==1))';
    Fluvial_Spatial(:,i) = fluvial(region_of_interest==1)';
    Glacial_Spatial(:,i) = abrasion(region_of_interest==1)';
    GOT_Spatial(:,i) = got(region_of_interest==1)';

    [~,~,elev_bin] = histcounts(bed(:),elevation_bins);

    Mean_GOT(i) = mean(got(region_of_interest==1));
    std_GOT(i) = std(got(region_of_interest==1));

    relief = localtopography(DEM,500);
    relief = flipud(relief.Z);
    relief = relief(region_of_interest == 1);

    Meanbed = localtopography(DEM,10000,'type','mean');
    Meanbed = flipud(Meanbed.Z);
    Meanbed = Meanbed(region_of_interest == 1);

    Mean_bed_Spatial(:,i) = Meanbed';
    bedTrue = bed;
    fluvial = fluvial(region_of_interest == 1);
    glacial = abrasion(region_of_interest == 1);
    hillslope = hillslope(region_of_interest == 1);

    Area_above_MeanElev(i) = length(bed > Meanbed) * cell_area;
    for j = 1:length(elevation_bins)
        index = find(elev_bin==j);
        % count slopes < 10° for the bin elevation
        slope_found = slopes(index);
        relief_found = relief(index);
        index_LRS = find(slope_found < Bed_slope_threshold & relief_found < Relief_threshold & bedTrue(index) > Meanbed(index));
        LRS_area_per_elevation_bin(j,i) = length(index_LRS)*dx*dy*1e-6; % Area of LRS (km²)
	    Hypsometry_area_per_elevation_bin(j,i) = length(index)*dx*dy*1e-6; % Area of surfaces (km²)

        fluvial_temp = fluvial(index); glacial_temp = glacial(index); hillslope_temp = hillslope(index);
        total_erosion = (sum(fluvial_temp(index_LRS)) + sum(glacial_temp(index_LRS)) + sum(hillslope_temp(index_LRS)));
        Erosion_per_elevation_bin(i,(j)) = mean(glacial_temp(index_LRS) + fluvial_temp(index_LRS) + hillslope_temp(index_LRS));
    end
    Total_LRS_area(i) = sum(LRS_area_per_elevation_bin(:,i));

    % loop through time
    for k = 1:number_of_output_files
        fnr = k;
        loaddata_models
        bed_include = bed(region_of_interest == 1); massb_include = massb(region_of_interest == 1); ice_include = ice(region_of_interest == 1);

        % Find min elevation glacier extent
        % We apply threshold on ice thickness to focus on larger glaciers
        index_extent = find(ice(region_of_interest == 1) > 100);
        if isempty(index_extent)
            maxglacierExtent(k,i) = max(bed_include);
        else
            maxglacierExtent(k,i) = min(bed_include(index_extent)+ice_include(index_extent));
        end
        ELA_array(k,i) = Ta(1,1)/SPM.mprop.lrate;
        maxTopo_evolution_models(k,i) = max(bed_include(:));
        meanTopo_evolution_models(k,i) = mean(bed_include(:));
        minTopo_evolution_models(k,i) = min(bed_include(:));
        isostasy_evolution_models(k,i) = mean(isostasy(region_of_interest==1));
        index_above_ELA = find(bed_include >= ELA_array(k,i));
        area_above_ELA(k,i) = length(index_above_ELA) .* cell_area;
        ice_vol_mean_interglacials(k,i) = mean(ice_include);
        Mean_water_discharge_evolution(k,i) = mean(qw(region_of_interest==1));

        DEM = GRIDobj(xc,yc,bed);

        relief = localtopography(DEM,1000);
        relief = flipud(relief.Z);
        bcurv = get_bed_curvature(bed,dx,dy);
        % get cirques density
        [cirques_density_time(k,i),~,~] = get_cirques_density(bcurv,relief,bed,curvature_threshold_for_cirques,Relief_threshold_for_cirques,abrasion,fluvial,abrasion_threshold_for_cirques,xc,yc,0);
    end
   
end
cd(path_to_codes)
cirques_density_time = cirques_density_time .* cell_area ./ 1e6;

save([path_to_save_analyses,'Models_results_part1.mat'],'Total_LRS_area','Erosion_models', 'ELA_array','std_ELA_array', ...
    'elevation_bins','LRS_area_per_elevation_bin','maxglacierExtent','Abrasion_models',...
    'Fluvial_models','Hillslope_models','isostasy_models','Area_above_MeanElev',...
    'Erosion_Spatial','Fluvial_Spatial','Glacial_Spatial','Mean_bed_Spatial','isostasy_evolution_models',...
    'minTopo_evolution_models','maxTopo_evolution_models','meanTopo_evolution_models','GOT_Spatial',...
    'Mean_GOT','std_GOT','ice_vol_mean_interglacials','area_above_ELA','Erosion_per_elevation_bin','bed_Spatial',...
    'Mean_water_discharge_evolution','Hypsometry_area_per_elevation_bin','cirques_density_time','Initial_hypsometry_area_per_elevation_bin',...
    'Initial_LRS_area_per_elevation_bin');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Extract results from time series %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time_from_time_series = 100:100:2.990e6;
ntime = length(time_from_time_series);
Mean_erosion_rate_evolution = zeros(ntime,number_of_Models);
Mean_fluvial_rate_evolution = zeros(ntime,number_of_Models);
Mean_glacial_rate_evolution = zeros(ntime,number_of_Models);
Temperature_evolution = zeros(ntime,number_of_Models);
ice_volume_evolution = zeros(ntime,number_of_Models);
Ice_vol_mean_stage1 = zeros(number_of_Models,1);
Ice_vol_mean_stage2 = zeros(number_of_Models,1);
for i = 1:number_of_Models
    disp(['Doing model ',num2str(i),'/',num2str(number_of_Models)])
    tdata = load([path_to_models,Models{i},'/output/tseries.dat']);
    time = tdata(:,1);
    [time_unique,index,ic] = unique(time);
    Temp = tdata(:,10);
    Temperature_evolution(:,i) = interp1(time_unique,Temp(index),time_from_time_series);
    ero_flu = tdata(:,36) .* 1e3; % mm/yr
    ero_gla = tdata(:,21) .* 1e3; % mm/yr
    ero_hill = tdata(:,19) .* 1e3; % mm/yr
    ero_peri = tdata(:,6) .* 1e3; % mm/yr
    ice_temp =  tdata(:,3) .* cell_area .* 1e-6;
    ice_volume_evolution(:,i) =  interp1(time_unique,ice_temp(index),time_from_time_series);
    Ice_vol_mean_stage1(i) = mean(ice_volume_evolution(1:17400,i));
    Ice_vol_mean_stage2(i) = mean(ice_volume_evolution(:,i));
    Mean_fluvial_rate_evolution(:,i) = interp1(time_unique,ero_flu(index),time_from_time_series);
    Mean_glacial_rate_evolution(:,i) = interp1(time_unique,ero_gla(index),time_from_time_series);
    Mean_erosion_rate_evolution(:,i) = interp1(time_unique,ero_gla(index) + ero_flu(index) + ero_hill(index) + ero_peri(index), time_from_time_series);
end

save([path_to_save_analyses,'Erosion_timeSeries_part1.mat'],'Mean_erosion_rate_evolution','Mean_glacial_rate_evolution','Mean_fluvial_rate_evolution', ...
    'Ice_vol_mean_stage1','Ice_vol_mean_stage2','time_from_time_series','Temperature_evolution','ice_volume_evolution')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LRS area %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize arrays
LRS_area_evolution = zeros(number_of_output_files,number_of_Models);
LRS_area_distribution_stage1 = zeros(number_of_Models,length(elevation_bins));

% Loop through models
cd(path_to_models)
for i = 1:number_of_Models
    disp(['Doing model ',num2str(i),'/',num2str(number_of_Models)])

    % Loop through model output files
    for t = 1:number_of_output_files
        fnr = t;
        Area_dist = zeros(1,length(elevation_bins));
        loaddata_models
        if abs(abs((xc(1,2)-xc(1,1)))-abs(yc(2,1)-yc(1,1)))>1e-9
		    x = dx/2:dx:dx*nx;
		    y = dx/2:dx:dx*ny;
		    [xc,yc] = meshgrid(x,y);
        end

        DEM = GRIDobj(xc,yc,bed);
        bed = bed(region_of_interest==1);
        slopes = bslope(region_of_interest==1);
    
        [~,~,elev_bin] = histcounts(bed(:),elevation_bins);
        % Get local relief 
        relief = localtopography(DEM,500);
        relief = flipud(relief.Z);
        relief = relief(region_of_interest==1);
        % Get mean topography within 10 km
        Meanbed = localtopography(DEM,10000,'type','mean');
        Meanbed = flipud(Meanbed.Z);
        Meanbed = Meanbed(region_of_interest==1);
        % Get low-relief surfaces
        bedTrue = bed;
        for j=1:length(elevation_bins)
	        index = find(elev_bin==j);
	        % Found slopes < 10° for the bin elevation
	        slope_found = slopes(index);
	        relief_found = relief(index);
            % LRS are surfaces with bed slope < 10°, with local relief <
            % 200 m and wich elevation is > than local mean topography
            % (filter valley floor)
	        index_slope = find(slope_found < 0.177 & relief_found< Relief_threshold & bedTrue(index) > Meanbed(index));
	        Area_dist(j) = length(index_slope) * cell_area * 1e-6; % LRS area in km2
            if t==85
                LRS_area_distribution_stage1(i,j) = length(index_slope)*dx*dy*1e-6;
            end
        end
        LRS_area_evolution(t,i) = sum(Area_dist);
    end
end
cd(path_to_codes)
save([path_to_save_analyses,'LRS_area_evolution_part1.mat'],'LRS_area_evolution','LRS_area_distribution_stage1')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% glacial/interglacial maximum  %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Find glacial/interglacial maximum %%%%%
Models_ID = 1:number_of_Models;
Mean_Erate_interglacials = zeros(size(Models));
Mean_Temp_interglacials = zeros(size(Models));
delta_max_Temp_array =  zeros(size(Models));
min_Sep = 70000; % 50 kyrs
min_Sep2 = 30000;
for i=1:length(Models_ID)
    % find local maxima
    if Models_ID(i) == 5
        ind = islocalmax(Temperature_evolution(:,Models_ID(i)),'MinSeparation',min_Sep2,'SamplePoints',time_from_time_series);
    else
        ind = islocalmax(Temperature_evolution(:,Models_ID(i)),'MinSeparation',min_Sep,'SamplePoints',time_from_time_series);
    end

    % find local minima
    if Models_ID(i) == 5
        ind_gla = islocalmin(Temperature_evolution(:,Models_ID(i)),'MinSeparation',min_Sep2,'SamplePoints',time_from_time_series);
    else
        ind_gla = islocalmin(Temperature_evolution(:,Models_ID(i)),'MinSeparation',min_Sep,'SamplePoints',time_from_time_series);
    end
    % keep only time > 1.7 Myr
    if Models_ID(i) == 7 % No cycle model, take interglacial from reference model
        ind = ind_previous;
        ind_gla = ind_gla_previous;
    else
        ind(time_from_time_series<=1.7.*1e6) = 0;
        ind = find(ind==1);
        ind = ind(1:end);
        ind_previous = ind;
        ind_gla(time_from_time_series<=1.7.*1e6) = 0;
        ind_gla = find(ind_gla==1);
        ind_gla = ind_gla(1:end);
        ind_gla_previous = ind_gla;
    end
    Mean_Erate_interglacials(Models_ID(i)) = mean(Mean_erosion_rate_evolution(ind,Models_ID(i)));
    Mean_Temp_interglacials(Models_ID(i)) = mean(Temperature_evolution(ind,Models_ID(i)));
    delta_max_Temp_array(Models_ID(i)) = abs(mean(Temperature_evolution(ind,Models_ID(i)))-mean(Temperature_evolution(ind_gla,Models_ID(i))));
%     figure
%     plot(time_from_time_series,Temperature_evolution(:,Models_ID(i)))
%     hold on
%     scatter(time_from_time_series(ind),Temperature_evolution(ind,Models_ID(i)),20,'r','fill')
%     scatter(time_from_time_series(ind_gla),Temperature_evolution(ind_gla,Models_ID(i)),20,'b','fill')
%     hold off
end

save([path_to_save_analyses,'Glacial_interglacials_part1.mat'],'delta_max_Temp_array','Mean_Temp_interglacials','Mean_Erate_interglacials')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Erosion vs elevation  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
minele = 0;
maxele = 1;
Nh = 30;
Norm_elevation_bins = linspace(minele,maxele,Nh+1);
Erosion_dist = zeros(length(Norm_elevation_bins),number_of_Models);
bed_dist = zeros(length(Norm_elevation_bins),number_of_Models);
for m=1:length(Models)
    bed_norm = bed_Spatial(:,m) - min(bed_Spatial(:,m));
    bed_norm = bed_norm ./ max(bed_norm);
    [~,~,elev_bin] = histcounts(bed_norm,Norm_elevation_bins);
    for i=1:length(Norm_elevation_bins)
        index = find(elev_bin==i);
        Erosion_dist(i,m) = sum(Erosion_Spatial(index,m))./length(index);
        bed_dist(i,m) = sum(bed_Spatial(index,m))./length(index);
    end
end

% Save data
save([path_to_save_analyses,'Erosion_vs_elevation_part1.mat'],'Erosion_dist','bed_dist')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% local vs large scale relief  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% plot local relief distribution %%%%%
range_local = 2000; % 2 km window
range_large = 20000; % 20 km window
fnr = 146; % Last output file

% Inititalize arrays
mean_local_relief = zeros(size(Models));
mean_local_relief_init = zeros(size(Models));
mean_large_relief = zeros(size(Models));
mean_large_relief_init = zeros(size(Models));
mean_erosion = zeros(size(Models));

% Loop through models
for m=1:length(Models)
    cd([path_to_models,Models{m}])

    loaddata
    xc = xc(region_of_interest==1); xc = reshape(xc,236,108);
    yc = yc(region_of_interest==1); yc = reshape(yc,236,108);

    %%%  on Initial topo %%%
    % local relief
    bed0 = SPM.data.bed; bed0 = bed0(region_of_interest == 1); bed0 = reshape(bed0,236,108); 
    DEM = GRIDobj(xc, yc, bed0);
    relief = localtopography(DEM, range_local);
    relief = flipud(relief.Z); 
    mean_local_relief_init(m) = median(relief(:));
    % large-scale relief
    relief = localtopography(DEM, range_large);
    relief = flipud(relief.Z); 
    mean_large_relief_init(m) = median(relief(:));

    %%% on Final topo %%%
    % local relief
    bed = bed(region_of_interest == 1); bed = reshape(bed, 236, 108);
    DEM = GRIDobj(xc, yc, bed);
    relief = localtopography(DEM, range_local);
    relief = flipud(relief.Z); 
	N = histc(relief(:),elevation_bins) / number_of_cells;
    mean_local_relief(m) = median(relief(:));
    % large-scale relief
    relief = localtopography(DEM,range_large);
    relief = flipud(relief.Z); 
    mean_large_relief(m) = median(relief(:));
    mean_erosion(m) = mean(abrasion(region_of_interest == 1) + fluvial(region_of_interest == 1) + hillslope(region_of_interest == 1));

    cd(path_to_codes)
end

save([path_to_save_analyses,'Relief_scales_part1.mat'],'mean_erosion','mean_local_relief','mean_large_relief_init','mean_local_relief_init','mean_large_relief')




