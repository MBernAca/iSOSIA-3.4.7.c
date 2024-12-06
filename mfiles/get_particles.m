function get_particles

%This function compute the particles transfer times and show the transient
%arrival of particle within the frontal moraine.


clear all;
close all;
%Get information from particles
SPM=SPMload;
ny = SPM.mesh.ny; nx=SPM.mesh.nx;
dx = SPM.mesh.dx;
dy = SPM.mesh.dy;
L = 29.5e+3 ; H = 14.1e+3;
nx = SPM.mesh.nx; ny = SPM.mesh.ny;
[Xc, Yc] = meshgrid([1:nx]*dx-dx/2,[1:ny]*dy-dy/2);
Xkm = Xc ./1000; Ykm = Yc./1000;
load('include.mat'); load('cmap.mat'); load('cmap_jet.mat');
fnr_vector = [11:95];

if SPM.hwprop.kmin == 0.0041 %Model M16
    xg = [21.150,24.650,24.650,21.150].*1e3; yg = [5.550,5.550,4.250,4.250].*1e3;
    if SPM.mesh.dopulse == 0
        disp('The total number of particles is high, we advice to run this function for the model M16 !')
    end
elseif SPM.hwprop.kmin == 0.0040 %Model M17
    xg = [23.550,25.650,25.650,23.550].*1e3; yg = [5.550,5.550,4.350,4.350].*1e3; 
    if SPM.mesh.dopulse == 0
        disp('The total number of particles is high, we advice to run this function for the model M17 !')
    end
else
    xg = [22.550,26.650,26.650,22.550].*1e3; yg = [5.550,5.550,4.350,4.350].*1e3;
    if SPM.mesh.dopulse == 0
        disp('The total number of particles is high, we advice to run this function for the model Tied_UniformPart !')
    end
end


%Time to have 90% of particles in the morraine
Tot_particles_vector = zeros(size(fnr_vector));
nb_particles_vector = zeros(size(fnr_vector));
hillslope_type_vector = zeros(size(fnr_vector));
glacial_type_vector = zeros(size(fnr_vector));
fnr = fnr_vector(end);
loaddata;loadparticles; inactivep = sum(inparticles(:));
inparticles_end = sum(inactivep(:));
nb_hillslope_end = length(typep(typep==-1));
nb_glacial_end = length(typep(typep==1));
pfound_end = find(xp >= xg(1)-dx./1e3/2 & xp <= xg(2)+dx./1e3/2 & yp <= yg(1)+dy./1e3/2 & yp >= yg(3)-dy./1e3/2);
fnr=fnr_vector(1);
loadparticles;
glacial_cells = length(ice(ice>5)); hillslope_cells = length(polyg_Tied(polyg_Tied==1)) - glacial_cells;
nb_hillslope =hillslope_cells - (hillslope_cells - nb_hillslope_end);
nb_glacial = glacial_cells - (glacial_cells - nb_glacial_end);
Tot_particles = sum(particles(:)) - inparticles_end;
%For time_transfert
Ttransfert_part = [];
Vtransfert = [];
bxtransfert_part = [];
bytransfert_part = [];

for i=1:length(fnr_vector)
    fnr = fnr_vector(i);
    loadparticles; loaddata;
    pfound = find(xp >= xg(1)-dx./1e3/2 & xp <= xg(2)+dx./1e3/2 & yp <= yg(1)+dy./1e3/2 & yp >= yg(3)-dy./1e3/2);
    nb_particles_vector(i) = numel(pfound);
    Total_type = typep(pfound);
    hillslope_type_vector(i) = length(Total_type(Total_type==-1));
    glacial_type_vector(i) = length(Total_type(Total_type==1));
    
    %for time transfer
    bx_pfound = bx(pfound); by_pfound = by(pfound); dl_pfound = dl(pfound); type_pfound= typep(pfound);
    pmoraine = [];
    for ip=1:length(bxtransfert_part)
        pmoraine(ip) = find(bx_pfound==bxtransfert_part(ip) & by_pfound==bytransfert_part(ip)); %find particles previously stored
    end
    bx_pfound(pmoraine) = []; by_pfound(pmoraine) = []; dl_pfound(pmoraine) = [];type_pfound(pmoraine) = [];
    bxtransfert_part = [bxtransfert_part, bx_pfound'];
    bytransfert_part = [bytransfert_part, by_pfound'];
    Ttransfert_part = [Ttransfert_part, zeros(size(by_pfound))'+ SPM.mesh.filetime*fnr-1000];
    Vtransfert = [Vtransfert,dl_pfound'./(SPM.mesh.filetime*fnr-1000)];
   
            
end
percent_particles = nb_particles_vector./Tot_particles .*100;
percent_hillslope = hillslope_type_vector./nb_particles_vector .*100; %Percent of the total particles in FM
percent_glacial = glacial_type_vector./nb_particles_vector .*100;
percent_hillslope2 = hillslope_type_vector./nb_hillslope.*100; %Percent of the total hillslope vs Total hillslope particle
percent_glacial2 = glacial_type_vector./nb_glacial.*100;

figure
hold on;
time_100 = fnr_vector.*SPM.mesh.filetime-1000;
save('UniformPart_100.mat','percent_particles','percent_hillslope','percent_glacial','percent_hillslope2','percent_glacial2','time_100');
plot(fnr_vector.*SPM.mesh.filetime-1000,percent_particles,'k','LineWidth',2);
plot(fnr_vector.*SPM.mesh.filetime-1000,percent_hillslope,'r','LineWidth',2);
plot(fnr_vector.*SPM.mesh.filetime-1000,percent_glacial,'b','LineWidth',2);
plot(fnr_vector.*SPM.mesh.filetime-1000,percent_glacial2,'b--','LineWidth',1);
plot(fnr_vector.*SPM.mesh.filetime-1000,percent_hillslope2,'r--','LineWidth',1);
hold off;
xlabel('Time (yr)'); ylabel('Amount of particles (%)');
ylim([0 100]);
legend('N_p_,_F_M/N_p','N_p_h_,_F_M/N_p_,_F_M','N_p_g_,_F_M/N_p_,_F_M','N_p_g_,_F_M/N_p_g','N_p_h_,_F_M/N_p_h')
set(gca,'FontSize',8,'FontName','Arial')
set(gcf,'Units','centimeters','Position',[0,0,12,10])

%------------------- Time transfert --------------------------------------
fnr = fnr_vector(end);
loaddata; loadparticles;

inmoraine = zeros(size(bytransfert_part));
for ip=1:length(bytransfert_part)
    inmoraine(ip) = find(bx==bxtransfert_part(ip) & by==bytransfert_part(ip));
end
bx_outFM = bx; bx_outFM(inmoraine) = [];
by_outFM = by; by_outFM(inmoraine) = [];
bx_FM = bx(inmoraine); by_FM = by(inmoraine);
    
figure
ncon =30;
elev = bed./1000; 
cblim = [0 5];
caxis(cblim)
cdata_M = get_color(elev,colormap(cmap_jet),cblim,1,1,elev,ice);
hp=surf(Xkm,Ykm,elev,'Cdata',cdata_M); hold on; shading interp;
[cdata_N, mindata,maxdata] = get_color(log10(Ttransfert_part),colormap(cmap_jet),cblim,0,0,elev,0);
scatter3(bx_outFM./1000,by_outFM./1000,zeros(size(bx_outFM))+4e3,1,'ko','filled');
scatter3(bx_FM./1000,by_FM./1000,zeros(size(bx_FM))+4e3,1,'fill','cdata',cdata_N); cb = colorbar('position',[0.85,0.55,0.02,0.4]);
vc = linspace(min(elev(:)),max(elev(:)),20);
contour3(Xkm,Ykm,elev,vc,'k','LineWidth',0.05,'Tag','hContour');
xlabel('Distance (km)'); ylabel('Distance (km)');
material([.4,.4,.4,2,0.1]); 
hold off;
set(gca,'dataaspectratio',[1,1,1/1]);
set(hp,'facelighting','gouraud','edgelighting','gouraud');
set(gca,'ambientlightcolor',[.9,.8,.9]);
caxis([0 5])
vt = [0,0.25,0.5,0.75,1.0];
set(cb,'ylim',[4,5],'ytick',vt+4,'yticklabel',10.^(mindata+vt*(maxdata-mindata)));
set(get(cb,'ylabel'),'string','Average transfer time (yr)');
light('position',[100,0,500]*1e3);
view(90, 90)

figure
ncon =30;
elev = bed./1000; 
cblim = [0 5];
caxis(cblim)
cdata_M = get_color(elev,colormap(cmap_jet),cblim,1,1,elev,ice);
hp=surf(Xkm,Ykm,elev,'Cdata',cdata_M); hold on; shading interp;
[cdata_N, mindata,maxdata] = get_color(log(Vtransfert),colormap(cmap_jet),cblim,0,0,elev,0);
scatter3(bx_outFM./1000,by_outFM./1000,zeros(size(bx_outFM))+4e3,1,'ko','filled');
scatter3(bx_FM./1000,by_FM./1000,zeros(size(bx_FM))+4e3,1,'fill','cdata',cdata_N); cb = colorbar('position',[0.85,0.55,0.02,0.4]);
vc = linspace(min(elev(:)),max(elev(:)),20);
contour3(Xkm,Ykm,elev,vc,'k','LineWidth',0.05,'Tag','hContour');
xlabel('Distance (km)'); ylabel('Distance (km)');
material([.4,.4,.4,2,0.1]); 
hold off;
set(gca,'dataaspectratio',[1,1,1/1]);
set(hp,'facelighting','gouraud','edgelighting','gouraud');
set(gca,'ambientlightcolor',[.9,.8,.9]);
caxis([0 5])
vt = [0,0.25,0.5,0.75,1.0];
set(cb,'ylim',[4,5],'ytick',vt+4,'yticklabel',exp(mindata+vt*(maxdata-mindata)));
set(get(cb,'ylabel'),'string','Average particle speed (m/yr)');
light('position',[100,0,500]*1e3);
view(90, 90)




    