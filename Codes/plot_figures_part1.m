
% This code plot main figures from model results extracted by the code
% 'Extract_model_results_part1.m'.


% Plotting paramaters 
path_to_models = '../Models/';
path_to_results = '../Models_output';
path_to_figures = '../Figures/';
path_to_save_analyses = '../Model_analyses/';
path_to_codes = '../Codes/';
Models = {'Temp_-2C','Reference','Temp_+2C','Ampx2','41_kyr_cycles','100_kyr_cycles',...
    'No_cycles','Prec_x1-5','Prec_x075'};
number_of_Models = length(Models);
Colors_Models = {'c';[0,0,0];'b';'r';[0,255,0]./255;[0,153,0]./255;[0.5,0.92,1];
    [255,153,153]./255;[1,0.5,1]};
Legend_Models = {'Colder (-2°C)';'Reference';'Warmer (+2°C)';'Amp. x2';'41 kyrs';'100 kyrs';'No cycles';'Prx1.5';'Prx0.75'};
time1 = [1:85].*0.0205;
time2 = [89:149].*0.02;
Model_time = [time1,time2];

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

%% Load results analyses
load([path_to_save_analyses,'Models_results_part1.mat'])
load([path_to_save_analyses, 'Erosion_timeSeries_part1.mat'])
load([path_to_save_analyses, 'LRS_area_evolution_part1.mat'])


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Figures main text %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%--------------------------------------------------------------------------
%-------------------------- Section 2.1. ----------------------------------
%--------------------------------------------------------------------------

% Convert benthic 18&O into temperature
% This equation uses simple assumption that the observed signal
% is only attributed to temperature change alone, with the effects of
% salinity and ice volume change ignored.
% Epstein et al. (1953) estimated that a 18&O increase of 0.22 pour mille 
% is equivalent to a cooling of 1°C.

load lisiecki.mat; % 1Kyr intervals
T0 = 5;
miny1 = 2.5;
maxy1 = 5.5;
miny2 = -2;
maxy2 = 7;
timeWin = 5;

% Smooth signal over 5 kyr
temp = reshape(d18O, timeWin,[]);
M5 = mean(temp,1);
Age5 = Age(1:timeWin:end);
M5 = movmean(d18O,timeWin);
Age5 = Age;
Temperature2 = 16.5 - 4.3.*M5 + 0.14.*M5.^2;

% Simplify glacial cycles
maxtime = 3e6; % Model time [yr]
tempdrop = 3; % Background cooling [°C]
time = linspace(0,maxtime ,2000); %time vector in years
Temp = zeros(size(time));
shift = 0;
for i=1:length(time) 
    if time(i) < 1.7e6
        Temp(i) = 8.5+ -1.*sawtooth(shift+2.*pi.*time(i)./41e3,0.5)-(tempdrop/maxtime)*time(i);  %sea level temp 
    elseif time(i) >= 1.7e6 & time(i) < 2e6
        Temp(i) = 8.5+ -1.*(1+((time(i)-1.5e6)/0.5e6)*1).*sawtooth(shift+2.*pi.*time(i)./100e3,0.8)-(tempdrop/maxtime)*time(i);
    elseif time(i) >=2e6
        Temp(i) = 8.5 + -1.*2*sawtooth(shift+2*pi*time(i)/100e3,0.8)-(2.5/maxtime)*time(i);
    end
end   

%%%% Plot figure %%%%
figure
set(gcf,'units','normalized','outerposition',[0 0 1 1])
% Benthic delta 18O
subplot(2,1,1)
hold on; box on; grid on;
plot(Age5,M5,'k-')
xlabel('Age (Ma)')
ylabel(['Benthic \delta^1^8O (',char(8240),')'])
set(gca, 'Ydir','reverse')
set(gca, 'Xdir','reverse')
xlim([0 3])
hold off;
% Simplified climatic cycles
subplot(2,1,2)
hold on; box on; grid on;
plot(Age5,T0+Temperature2,'k-')
plot((maxtime-time)./1e6,Temp)
hold off
xlabel('Model time (Myr)')
ylabel('Temperature (°C)')
set(gca, 'Xdir','reverse')
xlim([0 3]), ylim([2 11])
xticks = get(gca,'XTick');
set(gca,'XTickLabel',3:-0.5:0)

% Save figure
print([path_to_figures, 'Figure1_part1'],'-dsvg','-vector')

%%%% Plot initial fluvial topography
cd([path_to_models,'Reference'])
show('file',0,'data','bed','noforeland','contours','ncon',10,'tlim',[0 1500], 'view',[90 90])
cd(['../',path_to_codes])




%%
%--------------------------------------------------------------------------
%-------------------------- Section 3.1. ----------------------------------
%--------------------------------------------------------------------------


% Figure 3. Formation and preservation of low-relief surfaces. (Reference model)

%%%%%%%%%%%%%%%%% Part 1: plot local relief distribution %%%%%%%%%%%%%%%%%%
Models_to_plot = [2];
size_font = 10;
range_local =1000; % 2 km window
fnr = 146; % Last output file

figure
set(gcf,'units','normalized','outerposition',[0 0 1 1])
cd([path_to_models,'Reference'])
subplot(3,1,1)
hypsometry('file',fnr,'relief','radius',range_local,'dlim',[0 1500],'nofig','noclose')
cd(['../',path_to_codes])

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
    [ELA_array,minElevGlacierExtent] = get_ELATime(1:85,interpFactor,Massb_threshold,IceThick_threshold);
    [Mean_ELA,index_ELA] = min(ELA_array(1:85)); 
    minElevGlacier = minElevGlacierExtent/1e3; 
    plot(ax,[-100,100],[Mean_ELA,Mean_ELA]./1e3,'b--')
    plot(ax,[-100,100],[minElevGlacier,minElevGlacier]./1e3,'r--')

    % Stage 2
    [ELA_array,minElevGlacierExtent] = get_ELATime(1:fnr,interpFactor,Massb_threshold,IceThick_threshold);
    [Mean_ELA,index_ELA] = min(ELA_array(1:fnr)); 
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
print([path_to_figures, 'Figure3_Part2_Hyspometry_LRS_vs_elevation_ReferenceModel'],'-dsvg','-vector')


%%%%%%%%%%% Part 3: Show topographies with LRS and cirques %%%%%%%%%%%%%%%%
figure 
set(gcf,'units','normalized','outerposition',[0 0 1 1])
t = tiledlayout(2,1);
nexttile
cd([path_to_models, Models{Models_to_plot(1)}])
show('file',fnr,'data','bed','tlim',[-800 1700], 'mask',0.01,'cirques','cirques','withLRS','cmap','vik','noforeland','view',[90 90],'noclose','nofig','nocbar','radius',5000)
set(gca,'Visible','off')
cd(['../',path_to_codes])

%%%%%%%%%%% Part 4: plot zoom erosion  %%%%%%%%%%%%%%%%
figure
t = tiledlayout(2,4);
cd([path_to_models, Models{Models_to_plot(1)}])
%%% Kg x50 %%%
% Sliding rate
nexttile
show('file',1:fnr,'data','bed','dlim',[0 100],'mask',0.001,'sliding','sliding','dlim',[0 10],'noforeland','zoom',[1,15,50,23,65,-1,2],'view',[90 90],'nofig','nocbar','noclose');%
set(gca,'Visible','off')
% Glacial erosion
nexttile
show('file',fnr,'data','abrasion_rate','rate',0,'dlim',[0 100],'dlim',[0 0.1],'noforeland','zoom',[1,15,50,23,65,-1,2],'view',[90 90],'nofig','nocbar','noclose','cmap','lapaz');%
set(gca,'Visible','off')
% Fluvial + hillslope erosion
nexttile
show('file',fnr,'data','fluvial_rate','rate',0,'dlim',[0 100],'dlim',[0 0.4],'noforeland','zoom',[1,15,50,23,65,-1,2],'view',[90 90],'nofig','nocbar','noclose','cmap','bamako');%
set(gca,'Visible','off')
% Mean erosion
nexttile
show('file',fnr,'data','erosion_rate','rate',0,'dlim',[0 100],'dlim',[0 0.5],'noforeland','zoom',[1,15,50,23,65,-1,2],'view',[90 90],'nofig','nocbar','noclose','cmap','lajolla');%
set(gca,'Visible','off')
cd(['../',path_to_codes])


%%
%--------------------------------------------------------------------------
%-------------------------- Section 3.2.1. --------------------------------
%--------------------------------------------------------------------------

% Figure 5. Formation and preservation of low-relief surfaces for climate
% variation.
%%%%%%%%%%%%%%%%% Part 1: plot local relief distribution %%%%%%%%%%%%%%%%%%
Models_to_plot = [3, 1, 7, 4];
size_font = 10;
range_local =1000; % 2 km window
fnr = 146; % Last output file

figure
set(gcf,'units','normalized','outerposition',[0 0 1 1])
cd([path_to_models,Models{Models_to_plot(1)}])
subplot(2,2,1)
hypsometry('file',fnr,'relief','radius',range_local,'dlim',[0 1500],'nofig','noclose')
cd(['../',path_to_codes])
subplot(2,2,2)
cd([path_to_models,Models{Models_to_plot(2)}])
hypsometry('file',fnr,'relief','radius',range_local,'dlim',[0 1500],'nofig','noclose')
cd(['../',path_to_codes])
subplot(2,2,3)
cd([path_to_models,Models{Models_to_plot(3)}])
hypsometry('file',fnr,'relief','radius',range_local,'dlim',[0 1500],'nofig','noclose')
cd(['../',path_to_codes])
subplot(2,2,4)
cd([path_to_models,Models{Models_to_plot(4)}])
hypsometry('file',fnr,'relief','radius',range_local,'dlim',[0 1500],'nofig','noclose')
cd(['../',path_to_codes])

% Save figure
print([path_to_figures,'Figure5_Part1_reliefDistributions'],'-dsvg','-vector')


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
    [ELA_array,minElevGlacierExtent] = get_ELATime(1:85,interpFactor,Massb_threshold,IceThick_threshold);
    [Mean_ELA,index_ELA] = min(ELA_array(1:85)); 
    minElevGlacier = minElevGlacierExtent/1e3; 
    plot(ax,[-100,100],[Mean_ELA,Mean_ELA]./1e3,'b--')
    plot(ax,[-100,100],[minElevGlacier,minElevGlacier]./1e3,'r--')

    % Stage 2
    [ELA_array,minElevGlacierExtent] = get_ELATime(1:fnr,interpFactor,Massb_threshold,IceThick_threshold);
    [Mean_ELA,index_ELA] = min(ELA_array(1:fnr)); 
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
print([path_to_figures, 'Figure3_Part2_Hyspometry_LRS_vs_elevation'],'-dsvg','-vector')


%%%%%%%%%%% Part 3: Show topographies with LRS and cirques %%%%%%%%%%%%%%%%
figure 
set(gcf,'units','normalized','outerposition',[0 0 1 1])
t = tiledlayout(4,1);
nexttile
cd([path_to_models, Models{Models_to_plot(1)}])
show('file',fnr,'data','bed','tlim',[-1000 1600], 'mask',0.01,'cirques','withLRS','cmap','vik','noforeland','view',[90 90],'noclose','nofig','nocbar','radius',1000)
set(gca,'Visible','off')
cd(['../',path_to_codes])

nexttile
cd([path_to_models, Models{Models_to_plot(2)}])
show('file',fnr,'data','bed','tlim',[-1000 1600], 'mask',0.01,'cirques','withLRS','cmap','vik','noforeland','view',[90 90],'noclose','nofig','nocbar','radius',1000)
set(gca,'Visible','off')
cd(['../',path_to_codes])

nexttile
cd([path_to_models, Models{Models_to_plot(3)}])
show('file',fnr,'data','bed','tlim',[-1000 1600], 'mask',0.01,'cirques','withLRS','cmap','vik','noforeland','view',[90 90],'noclose','nofig','nocbar','radius',1000)
set(gca,'Visible','off')
cd(['../',path_to_codes])

nexttile
cd([path_to_models, Models{Models_to_plot(4)}])
show('file',fnr,'data','bed','tlim',[-1000 1600], 'mask',0.01,'cirques','withLRS','cmap','vik','noforeland','view',[90 90],'noclose','nofig','nocbar','radius',1000)
set(gca,'Visible','off')
cd(['../',path_to_codes])


%%%%%%%%%%% Part 4: plot zoom erosion  %%%%%%%%%%%%%%%%
figure
t = tiledlayout(4,4);
cd([path_to_models, Models{Models_to_plot(1)}])
%%% Warmer model %%%
% Sliding rate
nexttile
show('file',1:fnr,'data','bed','dlim',[0 100],'mask',0.001,'sliding','sliding','dlim',[0 10],'noforeland','zoom',[1,15,50,23,65,-1,2],'view',[90 90],'nofig','nocbar','noclose');%
set(gca,'Visible','off')
% Glacial erosion
nexttile
show('file',fnr,'data','abrasion_rate','rate',0,'dlim',[0 100],'dlim',[0 0.1],'noforeland','zoom',[1,15,50,23,65,-1,2],'view',[90 90],'nofig','nocbar','noclose','cmap','lapaz');%
set(gca,'Visible','off')
% Fluvial + hillslope erosion
nexttile
show('file',fnr,'data','fluvial_rate','rate',0,'dlim',[0 100],'dlim',[0 0.4],'noforeland','zoom',[1,15,50,23,65,-1,2],'view',[90 90],'nofig','nocbar','noclose','cmap','bamako');%
set(gca,'Visible','off')
% Mean erosion
nexttile
show('file',fnr,'data','erosion_rate','rate',0,'dlim',[0 100],'dlim',[0 0.5],'noforeland','zoom',[1,15,50,23,65,-1,2],'view',[90 90],'nofig','nocbar','noclose','cmap','lajolla');%
set(gca,'Visible','off')
cd(['../',path_to_codes])

%%% Colder model %%%
cd([path_to_models, Models{Models_to_plot(2)}])
% Sliding rate
nexttile
show('file',1:fnr,'data','bed','dlim',[0 100],'mask',0.001,'sliding','sliding','dlim',[0 10],'noforeland','zoom',[1,15,50,23,65,-1,2],'view',[90 90],'nofig','nocbar','noclose');%
set(gca,'Visible','off')
% Glacial erosion
nexttile
show('file',fnr,'data','abrasion_rate','rate',0,'dlim',[0 100],'dlim',[0 0.1],'noforeland','zoom',[1,15,50,23,65,-1,2],'view',[90 90],'nofig','nocbar','noclose','cmap','lapaz');%
set(gca,'Visible','off')
% Fluvial + hillslope erosion
nexttile
show('file',fnr,'data','fluvial_rate','rate',0,'dlim',[0 100],'dlim',[0 0.4],'noforeland','zoom',[1,15,50,23,65,-1,2],'view',[90 90],'nofig','nocbar','noclose','cmap','bamako');%
set(gca,'Visible','off')
% Mean erosion
nexttile
show('file',fnr,'data','erosion_rate','rate',0,'dlim',[0 100],'dlim',[0 0.5],'noforeland','zoom',[1,15,50,23,65,-1,2],'view',[90 90],'nofig','nocbar','noclose','cmap','lajolla');%
set(gca,'Visible','off')
cd(['../',path_to_codes])

%%% No cycles model %%%
cd([path_to_models, Models{Models_to_plot(3)}])
% Sliding rate
nexttile
show('file',1:fnr,'data','bed','dlim',[0 100],'mask',0.001,'sliding','sliding','dlim',[0 10],'noforeland','zoom',[1,15,50,23,65,-1,2],'view',[90 90],'nofig','nocbar','noclose');%
set(gca,'Visible','off')
% Glacial erosion
nexttile
show('file',fnr,'data','abrasion_rate','rate',0,'dlim',[0 100],'dlim',[0 0.1],'noforeland','zoom',[1,15,50,23,65,-1,2],'view',[90 90],'nofig','nocbar','noclose','cmap','lapaz');%
set(gca,'Visible','off')
% Fluvial + hillslope erosion
nexttile
show('file',fnr,'data','fluvial_rate','rate',0,'dlim',[0 100],'dlim',[0 0.4],'noforeland','zoom',[1,15,50,23,65,-1,2],'view',[90 90],'nofig','nocbar','noclose','cmap','bamako');%
set(gca,'Visible','off')
% Mean erosion
nexttile
show('file',fnr,'data','erosion_rate','rate',0,'dlim',[0 100],'dlim',[0 0.5],'noforeland','zoom',[1,15,50,23,65,-1,2],'view',[90 90],'nofig','nocbar','noclose','cmap','lajolla');%
set(gca,'Visible','off')
cd(['../',path_to_codes])

%%% Amplitude x2 model %%%
cd([path_to_models, Models{Models_to_plot(4)}])
% Sliding rate
nexttile
show('file',1:fnr,'data','bed','dlim',[0 100],'mask',0.001,'sliding','sliding','dlim',[0 10],'noforeland','zoom',[1,15,50,23,65,-1,2],'view',[90 90],'nofig','nocbar','noclose');%
set(gca,'Visible','off')
% Glacial erosion
nexttile
show('file',fnr,'data','abrasion_rate','rate',0,'dlim',[0 100],'dlim',[0 0.1],'noforeland','zoom',[1,15,50,23,65,-1,2],'view',[90 90],'nofig','nocbar','noclose','cmap','lapaz');%
set(gca,'Visible','off')
% Fluvial + hillslope erosion
nexttile
show('file',fnr,'data','fluvial_rate','rate',0,'dlim',[0 100],'dlim',[0 0.4],'noforeland','zoom',[1,15,50,23,65,-1,2],'view',[90 90],'nofig','nocbar','noclose','cmap','bamako');%
set(gca,'Visible','off')
% Mean erosion
nexttile
show('file',fnr,'data','erosion_rate','rate',0,'dlim',[0 100],'dlim',[0 0.5],'noforeland','zoom',[1,15,50,23,65,-1,2],'view',[90 90],'nofig','nocbar','noclose','cmap','lajolla');%
set(gca,'Visible','off')
cd(['../',path_to_codes])

%%
%--------------------------------------------------------------------------
%-------------------------- Section 3.2.2 ---------------------------------
%--------------------------------------------------------------------------

% Figure 6. Time evolution of low-relief surfaces.

Model_ID = [1,5,8,2,6,4,3];

fig = figure;
set(fig,'defaultAxesColorOrder',[[0,0,0]; [0,0,0]]);
set(gcf,'units','normalized','outerposition',[0 0 1 1])
linewidth = 1.5;
% LRS production rate evolution - first 4 models (Part 1)
for m=1:4
    subplot(2,2,m)
    hold on, box on, grid on
    yyaxis right
    plot(time_from_time_series./1e6,Mean_erosion_rate_evolution(:,Model_ID(m)),'-','Color',[0.4,0.4,0.4],'LineWidth',1) % erosion rate
    plot(time_from_time_series./1e6,movmean(Mean_erosion_rate_evolution(:,Model_ID(m)),4000),'-','Color','k','LineWidth',1) % moving mean erosion rate
    ylabel('Erosion rate (mm/yr)'), ylim([0 0.6])

    yyaxis left
    dt = (Model_time(2:end)-Model_time(1:end-1))./2;
    plot(Model_time(1:end-1)+dt,(LRS_area_evolution(2:end,Model_ID(m))-LRS_area_evolution(1:end-1,Model_ID(m)))./(dt'.*2.*1e3),'-','Color','r','LineWidth',linewidth)
    plot([0 3],[0 0],'--','Color',[0.3 0.3 0.3])
    xlabel('Model time (Myr)');ylabel('Prod. rate of LRS (km²)')
    hold off

    set(gca,'fontsize',12,'fontname','arial')
    ylim([-2 2])
    title(Models{Model_ID(m)})
end
% Save figure
print([path_to_figures,'./Figure6_part1'],'-dsvg','-vector')

% LRS production rate evolution - Next 4 models (Part 2)
fig = figure;
set(fig,'defaultAxesColorOrder',[[0,0,0]; [0,0,0]]);
set(gcf,'units','normalized','outerposition',[0 0 1 1])
linewidth = 1.5;
for m=1:3
    k=m+4;
    subplot(2,2,m)
    hold on, box on, grid on
    yyaxis right
    plot(time_from_time_series./1e6,Mean_erosion_rate_evolution(:,Model_ID(k)),'-','Color',[0.4,0.4,0.4],'LineWidth',1) % erosion rate
    plot(time_from_time_series./1e6,movmean(Mean_erosion_rate_evolution(:,Model_ID(k)),4000),'-','Color','k','LineWidth',1) % moving mean erosion rate
    ylabel('Erosion rate (mm/yr)'), ylim([0 0.6])

    yyaxis left
    dt = (Model_time(2:end)-Model_time(1:end-1))./2;
    plot(Model_time(1:end-1)+dt,(LRS_area_evolution(2:end,Model_ID(k))-LRS_area_evolution(1:end-1,Model_ID(k)))./(dt'.*2.*1e3),'-','Color','r','LineWidth',linewidth)
    plot([0 3],[0 0],'--','Color',[0.3 0.3 0.3])
    xlabel('Model time (Myr)');ylabel('Prod. rate of LRS (km²)')
    hold off
    set(gca,'fontsize',12,'fontname','arial')
    ylim([-2 2])
    title(Models{Model_ID(k)})
end

% LRS area vs time
subplot(2,2,4)
Model_ID = [1,2,3,4,5,6,8];
hold on, box on, grid on;
for m=1:length(Model_ID)
    plot(Model_time,LRS_area_evolution(:,Model_ID(m)),'-','LineWidth',linewidth,'Color',Colors_Models{Model_ID(m)})
end
hold off
xlabel('Model time (Myr)');ylabel('LRS area (km²)')
set(gca,'fontsize',12,'fontname','arial')
ylim([0 1600])
legend(Legend_Models{Model_ID},'location','northwest')

% Save figure
print([path_to_figures,'./Figure6_part2'],'-dsvg','-vector')



%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Section 4.3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Figure 8. Model comparison with observations from Scandinavia.
% -------------------------------------------------------------------------

Model_ID = 1;
load([path_to_models,Models{Model_ID},'/input/SPM.mat'])
nx = 108; ny =236;
dx = SPM.mesh.dx; dy = SPM.mesh.dy;
x = linspace(0,468.75*(nx-1),108); 
y = linspace(0,468.75*(ny-1),236);
[X, Y] = meshgrid(x,y);
filetime = SPM.mesh.filetime;
include = SPM.data.include;
%%%%% Coordinates of samples %%%%%
load([path_to_models,Models{Model_ID},'/cursor_info_spatial.mat'])
% load(['./',Models{Model_ID},'/sample_cirques.mat'])
coord = cursor_info;
xcoord = zeros(length(coord),1);
ycoord = xcoord;
for i=1:length(coord)
    xcoord(i) = coord(i).Position(1).*1e3;
    ycoord(i) = coord(i).Position(2).*1e3;
end
%Find closest grid point from each coordinates
ind = zeros(size(xcoord));
for k=1:length(xcoord)
    ind(k) = find(abs(xcoord(k)-X) < dx/2 & abs(ycoord(k)-Y) < dy/2 );
end
index_matrix = zeros(size(X));
count = 1;
for i=1:size(X,1)
    for j=1:size(X,2)
        index_matrix(i,j) = count;
        count = count +1;
    end
end

%%% End Coordinates %%%%
nx = SPM.mesh.nx; ny =SPM.mesh.ny;
fnr_range = 1:number_of_output_files;
linecolors = lines(length(xcoord));
Ero_data = zeros(length(xcoord),length(fnr_range));
EroRate_data = Ero_data;
bed_data = Ero_data;
bed0_data = Ero_data;
Glacial_data = Ero_data;
sliding_data = Ero_data;
Fluvial_data = Ero_data;
isostasy_data = Ero_data;
relief_data = Ero_data;
LRS_data = Ero_data;
ELA_array = zeros(1,number_of_output_files);
std_array = zeros(1,number_of_output_files);
%Initial bed
bed0 = SPM.data.bed;
bed0 = bed0(include==1);
%Go through time
ero_prev = SPM.data.abrasion;
gla_prev = SPM.data.abrasion;
flu_prev = SPM.data.fluvial;
isos_prev = SPM.data.abrasion;
i = Model_ID;
for f=1:length(fnr_range)
    fnr = fnr_range(f);
    loaddata_models
    %Map LRS at sample location
    [LRS,relief] = map_LRS(xc,yc,bed,bslope,500,10000,include);

    erosion = abrasion+hillslope+fluvial;
    erosion = erosion(include==1);
    fluvial = fluvial(include==1);
    abrasion = abrasion(include==1);
    hillslope = hillslope(include==1);
    bed = bed(include==1);
    ReliefvsIcevsEro_part1 = isostasy(include==1);
    for j=1:length(xcoord)
        EroRate_data(j,f) = (erosion(ind(j))-ero_prev(ind(j)))./SPM.mesh.filetime.*1e3;
		Ero_data(j,f) = erosion(ind(j));
        Glacial_data(j,f) = (abrasion(ind(j))- gla_prev(ind(j)))./SPM.mesh.filetime.*1e3;
		Fluvial_data(j,f) = ((fluvial(ind(j))+hillslope(ind(j))) - flu_prev(ind(j)))./SPM.mesh.filetime.*1e3;
        sliding_data(j,f) = sliding(ind(j));
        bed_data(j,f) = bed(ind(j));
        bed0_data(j,f) = bed0(ind(j));
        isostasy_data(j,f) = isostasy(ind(j))+SPM.data.srate(50,50).*SPM.mesh.filetime.*f;
        LRS_data(j,f) = LRS(ind(j));
        relief_data(j,f) = relief(ind(j));
    end
	ero_prev = erosion;
    flu_prev = fluvial+hillslope;
    gla_prev = abrasion;
    isos_prev = isostasy+SPM.mesh.filetime.*SPM.data.srate(50,50);
    % get ELA
    index = find(massb(include==1) < 0.05 & massb(include==1) > -0.05 & ice(include==1) > 10);
    bed_include = bed; massb_include = massb(include==1); ice_include = ice(include==1);
    ELA_array(f) = mean(bed_include(index)+ice_include(index));
    if isnan(ELA_array(f)) %means ELA is above maximum topography
        ELA_array(f) = max(bed(:));
    end
    std_array(f) = std(bed_include(index)+ice_include(index));
end

% Plot sample location
figure
cd([path_to_models,Models{Model_ID}])
show('file',max(fnr_range),'data','bed','tlim',[-1000 1500],'mask2',0.177,'radius',10000,'noforeland','zoom',[1,15,50,23,65,-1,2],'view',[63 50],'nofig','nocbar','noclose','contours','ncon',10);%
hold on
for i=1:length(ind)
    scatter3(gca,X(ind(i))./1e3,Y(ind(i))./1e3,bed(ind(i))./1e3,30,linecolors(i,:),'filled','MarkerEdgeColor','k')
end
hold off
set(gca,'Visible','off'), light
cd(['../',path_to_codes])

figure
set(gcf,'units','normalized','outerposition',[0 0 1 1])
%Climate
subplot(6,4,[2,3])
tdata = load(['./',Models{Model_ID},'./output/tseries.dat']);
timeMa = tdata(:,1)./1e6;
Climate = tdata(:,10);
hold on, box on, grid on,
plot(timeMa,Climate,('r-'))
plot([1.7, 1.7], [2, 10],'--','LineWidth',0.5,'Color',[0.4,0.4,0.4])
hold off
ylim([1, 10])
ylabel('T_0 (°C)')
set(gca,'fontsize',9,'fontname','arial')

% LRS evolution
subplot(6,4,[6,7])
hold on, grid on ,box on,
yyaxis left
plot(Model_time(2:end)-(Model_time(2:end)-Model_time(1:end-1))./2,(LRS_area_evolution(2:end,Model_ID)-LRS_area_evolution(1:end-1,Model_ID))./20,'k-','LineWidth',0.5)
ylabel('LRS production rate (km²/kyr)')
ylim([-2 2])
yyaxis right
plot(Model_time,LRS_area_evolution(:,Model_ID),'k-','LineWidth',1)
ylabel('Total LRS area (km²)')
hold off
ylim([0 1600])
set(gca,'fontsize',9,'fontname','arial')
% Erosion
subplot(6,4,[18,19])
hold on,box on, grid on;
for i=1:length(xcoord)
    plot(Model_time,Ero_data(i,:),'color',linecolors(i,:))
    plot(Model_time,isostasy_data(i,:),'--','color',linecolors(i,:))
end
hold off
ylabel('Erosion (m)')
set(gca,'fontsize',9,'fontname','arial')

% erosion rate
subplot(6,4,[10,11])
plots = [];
hold on,box on, grid on;
for i=1:length(xcoord)
    plot(Model_time,Fluvial_data(i,:),'-','color',linecolors(i,:))
end
hold off
ylabel('E_F_H (mm/yr)')
legend(plots,'location','northeast')
set(gca,'fontsize',9,'fontname','arial')
ylimF = get(gca,'ylim');
ylim([0 0.9])

% glacial rate
subplot(6,4,[14,15])
plots = [];
hold on,box on, grid on;
for i=1:length(xcoord)
    plot(Model_time,Glacial_data(i,:),'-','color',linecolors(i,:))
end
ylabel('E_g (mm/yr)')
% ylim(ylimF)
ylim([0 1.0])
set(gca,'fontsize',9,'fontname','arial')

% Elevation
subplot(6,4,[22,23])
hold on,box on, grid on,
for i=1:length(xcoord)
    plot(Model_time,bed_data(i,:),'color',linecolors(i,:))
end
hold off
ylabel('Elevation (m)')
xlabel('Model time (Myr)')
% legend('location','northwest')
set(gca,'fontsize',9,'fontname','arial')
ylim([-1000 2000])



