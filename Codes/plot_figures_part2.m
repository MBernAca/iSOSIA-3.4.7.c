
% This code plot main figures from model results extracted by the code
% 'Extract_model_results_part2.m'.


% Plotting paramaters
path_to_models = '../Models/';
path_to_results = '../Models_output';
path_to_figures = '../Figures/';
path_to_save_analyses = '../Model_analyses/';
path_to_codes = '../Codes/';
Models = {'Kg_x50','Kg_x8','Kf_01','Kf_06','Exponent_2','No_ice','Temp_+1C',...
    'Temp_-1C','Prec_x09','Prec_x1-25'};
number_of_Models = length(Models);
Colors_Models = {[0,0,0];'b';'r';[0,255,0]./255;[0,153,0]./255;[0.5,0.92,1];
    [255,153,153]./255;[1,0.5,1];[255,0,255]./255};
Legend_Models = {'Reference';'K_Gx50';'K_Gx8';'K_Fx0.5';'K_Fx0.08';'l = 2';'No ice';...
    '+1°C';'-1°C';'Pr x0.9';'Pr x1.25'};

% Load model characteritics
load([path_to_models,Models{1},'/input/SPM.mat'])
number_of_output_files = 149;
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
Model_time = [1:149].*0.02;


%% Load results analyses
load([path_to_save_analyses,'Models_results_part2.mat'])
load([path_to_save_analyses, 'Erosion_timeSeries_part2.mat'])
load([path_to_save_analyses, 'LRS_area_evolution_part2.mat'])
load([path_to_save_analyses, 'reference_model_results.mat'])


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Figures main text %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%--------------------------------------------------------------------------
%-------------------------- Section 3.1. ----------------------------------
%--------------------------------------------------------------------------

%% Figure 2. Formation of low-relief surfaces.

Models_to_plot = [1, 2, 3, 4];
kg = [5e-4, 8e-5, 1e-5, 1e-5];
kf = [1.2, 1.2, 0.1, 0.6];
size_font = 10;

figure
set(gcf,'units','normalized','outerposition',[0 0 1 1])
cmap = colormap('cool');
% KG vs LRS area
subplot(2,2,1)
hold on, grid on, box on
scatter(kg(1:2),Total_LRS_area(Models_to_plot(1:2)),100,Total_LRS_area(Models_to_plot(1:2)),'o','fill'),cb = colorbar;
scatter(kg(3),LRS_area_reference(end),100,LRS_area_reference(end),'o','fill'),
hold off
caxis([0 766])
ylim([0 900]), xlim([0 6e-4])
xlabel('K_G')

% Kf vs LRS area
subplot(2,2,2)
hold on, grid on, box on
scatter(kf(3:4),Total_LRS_area(Models_to_plot(3:4)),100,Total_LRS_area(Models_to_plot(3:4)),'s','fill'),
scatter(kf(1),LRS_area_reference(end),100,LRS_area_reference(end),'s','fill'),
hold off
cb = colorbar;
ylabel('LRS area (km²)'), xlabel('K_F (m^{0.5})')
axis square
caxis([0 766])
ylim([0 900]), xlim([0 1.5])
set(gca,'fontname','arial','fontsize',size_font)

% Cumulative erosion vs normalized elevation
subplot(2,2,3)
hold on, box on, grid on,
minele = 0;
maxele = 1;
Nh = 30;
Norm_elevation_bins = linspace(minele,maxele,Nh+1);
Erosion_dist = zeros(length(Norm_elevation_bins),length(Models_to_plot)+1);
bed_dist = zeros(length(Norm_elevation_bins),length(Models_to_plot)+1);
for m=1:length(Models_to_plot)
    bed_norm = bed_Spatial(:,Models_to_plot(m)) - min(bed_Spatial(:,Models_to_plot(m)));
    bed_norm = bed_norm ./ max(bed_norm);
    [~,~,elev_bin] = histcounts(bed_norm,Norm_elevation_bins);
    for i=1:length(Norm_elevation_bins)
        index = find(elev_bin==i);
        Erosion_dist(i,m) = sum(Erosion_Spatial(index,Models_to_plot(m)))./length(index);
        bed_dist(i,m) = sum(bed_Spatial(index,Models_to_plot(m)))./length(index);
    end
end
% Compute for reference model
bed_norm = bed_spatial_reference - min(bed_spatial_reference);
bed_norm = bed_norm ./ max(bed_norm);
[~,~,elev_bin] = histcounts(bed_norm,Norm_elevation_bins);
for i=1:length(Norm_elevation_bins)
    index = find(elev_bin==i);
    Erosion_dist(i,m+1) = sum(Erosion_Spatial_reference(index))./length(index);
    bed_dist(i,m+1) = sum(bed_Spatial(index))./length(index);
end

% plot figure
cmap = colormap('cool');
normalized2cmap = floor(Total_LRS_area(Models_to_plot)./max(Total_LRS_area(Models_to_plot)).*256);
legendentry = cell(size(Models_to_plot));
hold on, grid on, box on
for m=1:length(Models_to_plot)
    p(m) = plot(Erosion_dist(:,Models_to_plot(m)),Norm_elevation_bins,'-','Color',cmap(normalized2cmap(m),:));
    legendentry{m} = num2str(round(Total_LRS_area(Models_to_plot(m)),1));
end
p(m+1) = plot(Erosion_dist(:,m+1),Norm_elevation_bins,'-','Color',cmap(normalized2cmap(m),:));
hold off
axis square
xlabel('Erosion (m)'),ylabel('Normalized elevation')
set(gca,'fontname','arial','fontsize',10)

% Cirques area vs LRS area
subplot(2,2,4)
hold on, grid on, box on,
scatter(cirques_density_time(end,Models_to_plot(1:2)),Total_LRS_area(Models_to_plot(1:2)),100,Total_LRS_area(Models_to_plot(1:2)),'o','fill'),
scatter(cirque_density_reference(end),LRS_area_reference(end),100,LRS_area_reference(end),'o','fill')
scatter(cirques_density_time(end,Models_to_plot(3:4)),Total_LRS_area(Models_to_plot(3:4)),100,Total_LRS_area(Models_to_plot(3:4)),'s','fill'),
scatter(cirques_density_time(end,5),Total_LRS_area(5),100,Total_LRS_area(5),'s','fill'),
hold off
cb = colorbar;
ylabel('LRS area (km²)'), xlabel('Cirque area (km²)')
axis square
caxis([0 766])
ylim([0 900])
set(gca,'fontname','arial','fontsize',10)

% Save figure 
print([path_to_figures,'Figure2_mainText'],'-dsvg','-vector')


%% Figure 3. Formation and preservation of low-relief surfaces.

%%%%%%%%%%%%%%%%% Part 1: plot local relief distribution %%%%%%%%%%%%%%%%%%
Models_to_plot = [1, 3];
size_font = 10;
range_local =1000; % 2 km window
fnr = 149; % Last output file

figure
set(gcf,'units','normalized','outerposition',[0 0 1 1])
cd([path_to_models,'Reference'])
subplot(3,1,1)
hypsometry('file',146,'relief','radius',range_local,'dlim',[0 1500],'nofig','noclose')
cd(['../',path_to_codes])

subplot(3,1,2)
cd([path_to_models,Models{Models_to_plot(1)}])
hypsometry('file',fnr,'relief','radius',range_local,'dlim',[0 1500],'nofig','noclose')
cd(['../',path_to_codes])

subplot(3,1,3)
cd([path_to_models,Models{Models_to_plot(2)}])
hypsometry('file',fnr,'relief','radius',range_local,'dlim',[0 1500],'nofig','noclose')
cd(['../',path_to_codes])

% Save figure
print([path_to_figures,'Figure3_Part1_reliefDistributions'],'-dsvg','-vector')


%%%%% Part 2: plot hypsometry and LRS distribution with elevation %%%%%%%%%
% get erosion with elevation
minele = -1000;
maxele = 2000;
edges_ero = linspace(minele,maxele,30+1);
Erosion_dist = zeros(length(edges_ero),length(Models_to_plot));
Fluvial_dist = zeros(length(edges_ero),length(Models_to_plot));
Abrasion_dist = zeros(length(edges_ero),length(Models_to_plot));
bed_dist = zeros(length(edges_ero),length(Models_to_plot));
% Loop through models
for m=1:length(Models_to_plot)
    bed_norm = bed_Spatial(:,Models_to_plot(m));
    [~,~,elev_bin] = histcounts(bed_norm,edges_ero);
    for i=1:length(edges_ero)
        index = find(elev_bin==i);
        Erosion_dist(i,Models_to_plot(m)) = sum(Erosion_Spatial(index,Models_to_plot(m)))./length(index);
        Fluvial_dist(i,Models_to_plot(m)) = sum(Fluvial_Spatial(index,Models_to_plot(m)))./length(index);
        Abrasion_dist(i,Models_to_plot(m)) = sum(Glacial_Spatial(index,Models_to_plot(m)))./length(index);
        bed_dist(i,Models_to_plot(m)) = sum(bed_Spatial(index,Models_to_plot(m)))./length(index);
    end
end

figure
set(gcf,'units','normalized','outerposition',[0 0 1 1])
for m=1:length(Models_to_plot)
    ax = subplot(2,2,m);
    hold on, box on, grid on
    %%% hypsometry frequency %%% 
    % Initial topo
    for i=1:number_of_elevation_bin
     elevation_range = [elevation_bins(i), elevation_bins(i), elevation_bins(i+1), elevation_bins(i+1)] ./ 1e3;
     hypsometry_range = [0, -Initial_hypsometry_area_per_elevation_bin(i), -Initial_hypsometry_area_per_elevation_bin(i), 0] ./ 500 .* 100;
     hi = patch(ax,hypsometry_range,elevation_range,[0.7,0.7,0.7], 'EdgeColor',[0.3,0.3,0.3],'linewidth',0.5);
    end

    % Final topo
    for i=1:number_of_elevation_bin
      elevation_range = [elevation_bins(i), elevation_bins(i), elevation_bins(i+1), elevation_bins(i+1)] ./ 1e3;
      hypsometry_range  = [0, -Hypsometry_area_per_elevation_bin(i,Models_to_plot(m)), -Hypsometry_area_per_elevation_bin(i,Models_to_plot(m)), 0] ./500 .* 100;
      he = patch(ax,hypsometry_range,elevation_range,[255,128,0]./255, 'EdgeColor',[0.3,0.3,0.3],'linewidth',0.5,'FaceAlpha',0.4);
    end

    %%% Cumulative erosion vs elevation %%% 
    temp =  Erosion_dist(:,Models_to_plot(m))./1e3;
    he = plot(-temp./1.5.*100,edges_ero./1e3,'k-','LineWidth',1);
    for i=1:number_of_elevation_bin
      elevation_range = [elevation_bins(i), elevation_bins(i), elevation_bins(i+1), elevation_bins(i+1)] ./ 1e3;
      % stage 1
      LRS_area_range = [0, LRS_area_distribution_stage1(Models_to_plot(m), i), LRS_area_distribution_stage1(Models_to_plot(m), i), 0] ./120 .* 100;
      patch(ax, LRS_area_range, elevation_range, [0.7,0.7,0.7], 'EdgeColor', [0.3,0.3,0.3], 'LineWidth', 0.5);
      % stage 2
      LRS_area_range = [0, LRS_area_per_elevation_bin(i,Models_to_plot(m)), LRS_area_per_elevation_bin(i,Models_to_plot(m)), 0] ./120.*100;
      hf = patch(ax, LRS_area_range, elevation_range, [0,230,0]./255, 'EdgeColor', [0.3,0.3,0.3], 'LineWidth', 0.5, 'FaceAlpha', 0.4);
    end 
    disp(['Total LRS area = ',num2str(sum(LRS_area_per_elevation_bin(:,Models_to_plot(m)))),' km²'])

    %%% Get Mean ELA %%%
    cd ([path_to_models,Models{Models_to_plot(m)}])

    interpFactor = 10;
    Massb_threshold = 0.1;
    IceThick_threshold = 10;
    % Stage 1
    [ELA_array,minElevGlacierExtent] = get_ELATime(1:87,interpFactor,Massb_threshold,IceThick_threshold);
    [Mean_ELA,index_ELA] = min(ELA_array(1:87)); 
    minElevGlacier = minElevGlacierExtent/1e3; 
    plot(ax,[-100,100],[Mean_ELA,Mean_ELA]./1e3,'b--')
    plot(ax,[-100,100],[minElevGlacier,minElevGlacier]./1e3,'r--')

    % Stage 2
    [ELA_array,minElevGlacierExtent] = get_ELATime(1:149,interpFactor,Massb_threshold,IceThick_threshold);
    [Mean_ELA,index_ELA] = min(ELA_array(1:149)); 
    minElevGlacier = minElevGlacierExtent/1e3; 
    plot(ax,[-100,100],[Mean_ELA,Mean_ELA]./1e3,'b-')
    plot(ax,[-100,100],[minElevGlacier,minElevGlacier]./1e3,'r-')

    hold off
    axis square
    xlim([-100 100]), ylim([-1 2]),
    ax.XTick = [-100, -50, 0, 50, 100];
    ax.XTickLabel = {'500','250','0','60','120'};
    ylabel('Elevation (km)'),
    set(gca,'fontname','arial','fontsize',10)
    title(Models{Models_to_plot(m)})
    cd(['../',path_to_codes])
end

% Save figure
print([path_to_figure, 'Figure3_Part2_Hyspometry_LRS_vs_elevation'],'-dsvg','-vector')


%%%%%%%%%%% Part 3: Show topographies with LRS and cirques %%%%%%%%%%%%%%%%
fnr = 149;
figure 
set(gcf,'units','normalized','outerposition',[0 0 1 1])
t = tiledlayout(2,1);
nexttile
cd([path_to_models, Models{Models_to_plot(1)}])
show('file',fnr,'data','bed','tlim',[-800 1700], 'mask',0.01,'cirques','cirques','withLRS','cmap','vik','noforeland','view',[90 90],'noclose','nofig','nocbar','radius',1000)
set(gca,'Visible','off')
cd(['../',path_to_codes])

nexttile
cd([path_to_models, Models{Models_to_plot(2)}])
show('file',fnr,'data','bed','tlim',[-800 1700], 'mask',0.01,'cirques','cirques','withLRS','cmap','vik','noforeland','view',[90 90],'noclose','nofig','nocbar','radius',1000)
set(gca,'Visible','off')
cd(['../',path_to_codes])


%%%%%%%%%%% Part 4: plot zoom erosion  %%%%%%%%%%%%%%%%
figure
t = tiledlayout(2,4);
cd([path_to_models, Models{Models_to_plot(1)}])
%%% Kg x50 %%%
% Sliding rate
nexttile
show('file',1:149,'data','bed','dlim',[0 100],'mask',0.001,'sliding','sliding','dlim',[0 10],'noforeland','zoom',[1,15,50,23,65,-1,2],'view',[90 90],'nofig','nocbar','noclose');%
set(gca,'Visible','off')
% Glacial erosion
nexttile
show('file',149,'data','abrasion_rate','rate',0,'dlim',[0 100],'dlim',[0 0.1],'noforeland','zoom',[1,15,50,23,65,-1,2],'view',[90 90],'nofig','nocbar','noclose','cmap','lapaz');%
set(gca,'Visible','off')
% Fluvial + hillslope erosion
nexttile
show('file',149,'data','fluvial_rate','rate',0,'dlim',[0 100],'dlim',[0 0.4],'noforeland','zoom',[1,15,50,23,65,-1,2],'view',[90 90],'nofig','nocbar','noclose','cmap','bamako');%
set(gca,'Visible','off')
% Mean erosion
nexttile
show('file',149,'data','erosion_rate','rate',0,'dlim',[0 100],'dlim',[0 0.5],'noforeland','zoom',[1,15,50,23,65,-1,2],'view',[90 90],'nofig','nocbar','noclose','cmap','lajolla');%
set(gca,'Visible','off')
cd(['../',path_to_codes])

%%% Kf x0.08 %%%
cd([path_to_models, Models{Models_to_plot(2)}])
% Sliding rate
nexttile
show('file',1:149,'data','bed','dlim',[0 100],'mask',0.001,'sliding','sliding','dlim',[0 10],'noforeland','zoom',[1,15,50,23,65,-1,2],'view',[90 90],'nofig','nocbar','noclose');%
set(gca,'Visible','off')
% Glacial erosion
nexttile
show('file',149,'data','abrasion_rate','rate',0,'dlim',[0 100],'dlim',[0 0.1],'noforeland','zoom',[1,15,50,23,65,-1,2],'view',[90 90],'nofig','nocbar','noclose','cmap','lapaz');%
set(gca,'Visible','off')
% Fluvial + hillslope erosion
nexttile
show('file',149,'data','fluvial_rate','rate',0,'dlim',[0 100],'dlim',[0 0.4],'noforeland','zoom',[1,15,50,23,65,-1,2],'view',[90 90],'nofig','nocbar','noclose','cmap','bamako');%
set(gca,'Visible','off')
% Mean erosion
nexttile
show('file',149,'data','erosion_rate','rate',0,'dlim',[0 100],'dlim',[0 0.5],'noforeland','zoom',[1,15,50,23,65,-1,2],'view',[90 90],'nofig','nocbar','noclose','cmap','lajolla');%
set(gca,'Visible','off')
cd(['../',path_to_codes])





