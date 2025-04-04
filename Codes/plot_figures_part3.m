% This code gather results from all models and plot main summary figures

path_to_models = '../Models/';
path_to_results = '../Models_output';
path_to_figures = '../Figures/';
path_to_save_analyses = '../Model_analyses/';
path_to_codes = '../Codes/';

Models_part1 = {'Temp_-2C','Reference','Temp_+2C','Ampx2','41_kyr_cycles','100_kyr_cycles',...
    'No_cycles','Prec_x1-5','Prec_x075'};
Models_part2 = {'Kg_x50','Kg_x8','Kf_01','Kf_06','Exponent_2','No_ice','Temp_+1C',...
    'Temp_-1C','Prec_x09','Prec_x1-25'};
Models = {'Temp_-2C','Reference','Temp_+2C','Ampx2','41_kyr_cycles','100_kyr_cycles',...
    'No_cycles','Prec_x1-5','Prec_x075','Kg_x50','Kg_x8','Kf_01','Kf_06','Exponent_2','No_ice','Temp_+1C',...
    'Temp_-1C','Prec_x09','Prec_x1-25'};

number_of_Models = length(Models);
Colors_Models = {[0.5,0.92,1]; [0,0,0]; [51,51,255]./255; [0,153,76]./255; [0,255,0]./255; [0,204,0]./255;
    [153,255,153]./255; [255,204,255]./255; [178,102,255]./255; [255,178,102]./255; [255,153,51]./255; [255,128,0]./255;
    [204,102,0]./255; [153,76,0]./255; [1,0,0]; [51,153,255]./255;...
    [153,204,255]./255; [255,102,255]./255; [255,153,255]./255};

Legend_Models = {'Colder (-2°C)';'Reference';'Warmer (+2°C)';'Amp. x2';'41 kyrs';'100 kyrs';'No cycles';'Prx1.5';...
    'Prx0.75';'K_Gx50';'K_Gx8';'K_F0.1';'K_Fx0.6';'l = 2';'No ice';...
    '+1°C';'-1°C';'Pr x0.9';'Pr x1.25'};
time1 = [1:85].*0.0205;
time2 = [89:149].*0.02;
Model_time1 = [time1,time2];
Model_time2 = [1:149].*0.02;


%% Initiate arrays

LRS_area_final = zeros(number_of_Models,1);
LRS_area_stage1 = zeros(number_of_Models,1);
Mean_ice_volume = zeros(number_of_Models,1);
T0_Initial = zeros(number_of_Models,1);
delta_Temp_glacialInterglacial = zeros(number_of_Models,1);
Mean_T0_interglacial = zeros(number_of_Models,1);
Mean_Erate_interglacial = zeros(number_of_Models,1);
Total_cirque_area = zeros(number_of_Models,1);
Erosion_vs_elevation = zeros(number_of_Models,31);
local_reliefs = zeros(number_of_Models,1);
local_reliefs_1km = zeros(number_of_Models,1);
large_reliefs = zeros(number_of_Models,1);
large_reliefs_init = zeros(number_of_Models,1);
local_reliefs_init = zeros(number_of_Models,1);
local_reliefs_init_1km = zeros(number_of_Models,1);
erosion_reliefs = zeros(number_of_Models,1);
Temperature_time = zeros(29900,number_of_Models);
ice_volume_time = zeros(29900,number_of_Models);
Mean_erosion_time = zeros(29900,number_of_Models);
Mean_water_discharge_time = zeros(146,number_of_Models);
LRS_time1 = zeros(146,length(Models_part1));

% Load result models part 1
load([path_to_save_analyses, 'Models_results_part1.mat'])
load([path_to_save_analyses, 'Erosion_timeSeries_part1.mat'])
load([path_to_save_analyses, 'LRS_area_evolution_part1.mat'])
load([path_to_save_analyses, 'Glacial_interglacials_part1.mat'])
load([path_to_save_analyses, 'Erosion_vs_elevation_part1.mat'])
load([path_to_save_analyses, 'Relief_scales_part1.mat'])

LRS_area_final(1:length(Models_part1)) = Total_LRS_area;
LRS_area_stage1(1:length(Models_part1)) = sum(LRS_area_distribution_stage1,2)';
Mean_ice_volume_stage1(1:length(Models_part1)) = Ice_vol_mean_stage1;
Mean_ice_volume_stage2(1:length(Models_part1)) = Ice_vol_mean_stage2;
for m=1:length(Models_part1)
    Temperature_time(:,m) = Temperature_evolution(:,m);
    ice_volume_time(:,m) = ice_volume_evolution(:,m);
    Mean_erosion_time(:,m) = Mean_erosion_rate_evolution(:,m);
    LRS_time1(:,m) = LRS_area_evolution(:,m);
%     Mean_water_discharge_time(:,m) = Mean_water_discharge_evolution(:,m);
end
T0_Initial(1:length(Models_part1)) = Temperature_evolution(1,:);
delta_Temp_glacialInterglacial(1:length(Models_part1)) = delta_max_Temp_array(:);
Total_cirque_area(1:length(Models_part1)) = cirques_density_time(end,:);
Erosion_vs_elevation(1:length(Models_part1),:) = Erosion_dist';
Mean_T0_interglacial(1:length(Models_part1))  = Mean_Temp_interglacials;
Mean_Erate_interglacial(1:length(Models_part1))  = Mean_Erate_interglacials;
local_reliefs(1:length(Models_part1)) = mean_local_relief;
local_reliefs_1km(1:length(Models_part1)) = median_local_relief_1km;
large_reliefs(1:length(Models_part1)) = mean_large_relief;
local_reliefs_init(1:length(Models_part1)) = mean_local_relief_init;
local_reliefs_init_1km(1:length(Models_part1)) = median_local_relief_init_1km;
large_reliefs_init(1:length(Models_part1)) = mean_large_relief_init;
erosion_reliefs(1:length(Models_part1)) = mean_erosion;

% Load result models part 2
load([path_to_save_analyses, 'Models_results_part2.mat'])
load([path_to_save_analyses, 'Erosion_timeSeries_part2.mat'])
load([path_to_save_analyses, 'LRS_area_evolution_part2.mat'])
load([path_to_save_analyses, 'Glacial_interglacials_part2.mat'])
load([path_to_save_analyses, 'Erosion_vs_elevation_part2.mat'])
load([path_to_save_analyses, 'Relief_scales_part2.mat'])
LRS_time2 = zeros(149,length(Models_part1));

indices = length(Models_part1) + 1 : length(Models_part1) + length(Models_part2)';
LRS_area_final(indices) = Total_LRS_area;
LRS_area_stage1(indices) = sum(LRS_area_distribution_stage1,2);
Mean_ice_volume_stage1(indices) = Ice_vol_mean_stage1;
Mean_ice_volume_stage2(indices) = Ice_vol_mean_stage2;
T0_Initial(indices) = Temperature_evolution(1,:);
delta_Temp_glacialInterglacial(indices) = delta_max_Temp_array(:);
Total_cirque_area(indices) = cirques_density_time(end,:);
Erosion_vs_elevation(indices,:) = Erosion_dist';
Mean_T0_interglacial(indices)  = Mean_Temp_interglacials;
Mean_Erate_interglacial(indices)  = Mean_Erate_interglacials;
local_reliefs(indices) = mean_local_relief;
local_reliefs_1km(indices) = median_local_relief_1km;
large_reliefs(indices) = mean_large_relief;
local_reliefs_init(indices) = mean_local_relief_init;
local_reliefs_init_1km(indices) = median_local_relief_init_1km;
large_reliefs_init(indices) = mean_large_relief_init;
erosion_reliefs(indices) = mean_erosion;
for m=1:length(Models_part2)
    Temperature_time(:,m+length(Models_part1)) = Temperature_evolution(:,m);
    ice_volume_time(:,m+length(Models_part1)) = ice_volume_evolution(:,m);
    Mean_erosion_time(:,m+length(Models_part1)) = Mean_erosion_rate_evolution(:,m);
    LRS_time2(:,m) = LRS_area_evolution(:,m);
%     Mean_water_discharge_time(:,m) = Mean_water_discharge_evolution(:,m);
end



%%
%--------------------------------------------------------------------------
%-------------------------- Section 3.2.1 ---------------------------------
%--------------------------------------------------------------------------

% Figure 5. Climatic control on the extent of ELRS after 3 Myr of
% glaciations

Models_to_plot = [1,2,3,16,17];
T0_init = T0_Initial(Models_to_plot);
size_font = 14;
MarkerSize = 100;

figure
set(gcf,'units','normalized','outerposition',[0 0 1 1])
% ELRS vs T0
ax1 = subplot(2,2,1);
hold on, box on, grid on;
scatter(T0_init,LRS_area_final(Models_to_plot),MarkerSize,Mean_ice_volume_stage2(Models_to_plot),'fill'), cb= colorbar(); colormap(flipud(cool));
hold off
axis square
ylabel(cb,'Mean ice volume (km^3)')
caxis([0 15])
xlabel('T_0 (°C)'), ylabel('LRS area (km²)')
set(gca,'fontname','arial','fontsize',size_font)
ylim([0 1600]), xlim([7, 12])

% ELRS vs Prec
Models_to_plot = [2, 8, 9, 18, 19];
Prec_rate = [1.0, 1.5, 0.75, 0.9, 1.25];
ax2 = subplot(2,2,2);
hold on, box on, grid on;
scatter(Prec_rate,LRS_area_final(Models_to_plot),MarkerSize,Mean_ice_volume_stage2(Models_to_plot),'fill','Marker','s'); cb= colorbar(); colormap(flipud(cool));
hold off
axis square
caxis([0 15])
ylabel(cb,'Mean ice volume (km^3)')
xlabel('Precipitation rate (mm/yr)'), ylabel('LRS area (km²)')
set(gca,'fontname','arial','fontsize',size_font)
ylim([0 1600]), xlim([0.7, 1.6])


% ELRS vs climate variability
Models_to_plot = [7, 5, 2, 6, 4];
ax3 = subplot(2,2,3);
hold on, box on, grid on;
scatter(1:5,LRS_area_final(Models_to_plot),MarkerSize,Mean_ice_volume_stage2(Models_to_plot),'fill'); cb= colorbar(); colormap(flipud(cool));
hold off
axis square
caxis([0 15])
ylabel(cb,'Mean ice volume (km^3)')
xlabel('Climate variability'), ylabel('LRS area (km²)')
ylim([0 1600]), xlim([0, 6])
set(gca,'fontname','arial','fontsize',size_font,'XTick',1:5,'XTickLabel',Legend_Models(Models_to_plot))

% ELRS vs mean ice volume
Models_to_plot = [1:9,16:19];
ax4 = subplot(2,2,4);
colormap(ax4,'parula')
hold on, box on, grid on;
scatter(Mean_ice_volume_stage1(Models_to_plot),LRS_area_stage1(Models_to_plot),MarkerSize,'fill','MarkerFaceColor',[0.5 0.5 0.5]);
scatter(Mean_ice_volume_stage2(Models_to_plot),LRS_area_final(Models_to_plot),100,delta_Temp_glacialInterglacial(Models_to_plot),'fill')
hold off
axis square
xlabel('Mean ice volume (km^3)'), ylabel('LRS area (km²)')
ylim([0 1600]), xlim([0, 15])
set(gca,'fontname','arial','fontsize',size_font)

% save figure
print([path_to_figures,'Figure5_mainText'],'-dsvg','-vector')



%%
%--------------------------------------------------------------------------
%---------------------------- Section 4.2 ---------------------------------
%--------------------------------------------------------------------------

% Figure 9. elevated low-relief surfaces and the occurence of the glacial shelter

Models_to_plot = [1:19];
minele = 0;
maxele = 1;
Nh = 30;
Norm_elevation_bins = linspace(minele,maxele,Nh+1);

figure
set(gcf,'units','normalized','outerposition',[0 0 1 1])
% Cumulative erosion vs normalized elevation
subplot(2,2,1)
hold on, box on, grid on,
cmap = flipud(colormap('cool'));
normalized2cmap = floor(LRS_area_final(Models_to_plot)./max(LRS_area_final(Models_to_plot)).*255)+1;
legendentry = cell(size(Models_to_plot));
hold on, grid on, box on
for m=1:length(Models_to_plot)
    p(m) = plot(Erosion_vs_elevation(Models_to_plot(m),:),Norm_elevation_bins,'-','Color',cmap(normalized2cmap(m),:));
    legendentry{m} = num2str(round(LRS_area_final(Models_to_plot(m)),1));
end
hold off
axis square
xlabel('Erosion (m)'),ylabel('Normalized elevation')
set(gca,'fontname','arial','fontsize',10)

% ELRS area vs Mean T0 during interglacial
subplot(2,2,2)
hold on, box on, grid on
scatter(Mean_T0_interglacial(Models_to_plot),LRS_area_final(Models_to_plot),100,Mean_Erate_interglacial(Models_to_plot),'fill')
cb = colorbar();colormap('parula')
caxis([0.15 0.4])
ylabel(cb,'Mean erosion rate during interglacials (mm/yr)')
hold off
xlim([5,11]), ylim([0 1600])
xlabel('Mean base-level temperature during interglacials (°C)'), ylabel('LRS area (km²)')
axis square
set(gca,'fontname','arial','fontsize',14)

% ELRS area vs valley-head cirques density
subplot(2,2,3)
hold on,box on, grid on
scatter(Total_cirque_area(Models_to_plot),LRS_area_final(Models_to_plot),100,Mean_ice_volume_stage2(Models_to_plot),'fill'), cb = colorbar;
hold off
ylabel('LRS area (km²)'), xlabel('Cirque area (km²)')
axis square
caxis([0 15])
ylim([0 1600]), xlim([0 70])
set(gca,'fontname','arial','fontsize',10)
ylabel(cb,'Mean ice volume (km^3)')

% Local vs large scale relief
subplot(2,2,4)
hold on, grid on, box on,
Relief_change_local = (local_reliefs-local_reliefs_init)./local_reliefs_init.*100;
Relief_change_large = (large_reliefs-large_reliefs_init)./large_reliefs_init.*100;
scatter(Relief_change_local(Models_to_plot),Relief_change_large(Models_to_plot),100,erosion_reliefs(Models_to_plot),'fill'), cb = colorbar();
plot([0 0],[0 60],'--','Color',[0.3 0.3 0.3],'LineWidth',1.5)
xlabel("Median 2 km relief change (%)"), ylabel('Median 20 km relief change (%)')
axis square
ylabel(cb,'Mean erosion (m)')
set(gca,'fontname','arial','fontsize',14)
hold off
xlim([-30 30]), ylim([0 60])

% save figure
print([path_to_figures,'Figure9_mainText'],'-dsvg','-vector')


%% Local relief for erosion efficacy variation models
Models_to_plot = [2,10,11];
kg = [1e-5, 5e-4, 8e-5];
figure
set(gcf,'units','normalized','outerposition',[0 0 1 1])
% vary kG
subplot(1,2,1)
hold on, grid on, box on,
cmap = colormap('cool');
scatter(kg, local_reliefs_1km(Models_to_plot),150,LRS_area_final(Models_to_plot),'fill'), cb = colorbar();
plot([0 0],[0 60],'--','Color',[0.3 0.3 0.3],'LineWidth',1.5)
xlabel("kg"), ylabel('Median 1 km relief (m)')
axis square
ylabel(cb,'LRS area (km2)')
set(gca,'fontname','arial','fontsize',14)
hold off
xlim([0 6e-4]),ylim([200 600])
caxis([0 800])

% vary Kf
Models_to_plot = [2,12,13];
kf = [1.1, 0.1, 0.6];
subplot(1,2,2)
cmap = colormap('cool');
hold on, grid on, box on,
scatter(kf, local_reliefs_1km(Models_to_plot),150,LRS_area_final(Models_to_plot),'s','fill'), cb = colorbar();
plot([0 0],[0 60],'--','Color',[0.3 0.3 0.3],'LineWidth',1.5)
xlabel("kF (m0.5)"), ylabel('Median 1 km relief (m)')
axis square
ylabel(cb,'LRS area (km2)')
set(gca,'fontname','arial','fontsize',14)
hold off
xlim([0 1.5]),ylim([200 600])
caxis([0 800])

print('./Figure3a_Local_relief_vs_erosion_efficacy','-dsvg','-vector')

