
function template()

close all;

%mesh dimensions for iSOSIA 900 m resolution
L = 50e3; H = 100e3;
nx =128; ny = 256;
dy = H/ny; dx = L/nx;
tfac = 100;

%load default SPM structure
SPM = setup(L,H,nx,ny);
SPM.mesh.hmin = min(SPM.mesh.dx,SPM.mesh.dy);
SPM.mesh.doice = 0;
SPM.mesh.dosliding = 0;  %sliding flag
SPM.mesh.doglacialerosion = 0;
SPM.mesh.doglacialhydrology = 0;
SPM.mesh.dosediment = 1;
SPM.mesh.doavalance = 0;
SPM.mesh.coldbased = 0;  %do not include coldbased ice
SPM.mesh.doglacialsedi = 0; %advect sediment volumes in glaciers
SPM.mesh.doparticles = 0;  
SPM.mesh.doisostasy = 1; %isostasy flag
SPM.mesh.dolandslide = 0;  %landslide flag

SPM.mesh.nmoni = 100;
SPM.mesh.maxtime = 60e6/tfac;  %length of simulation
SPM.mesh.filetime =2e6/tfac;  %time between outputs
SPM.mesh.maxdt = 1;    %maximum time step
SPM.mesh.ct = 0.1;       %time step scaling factor
SPM.mesh.slidingmode = 3; %0:weertman, 1:budd , 2:Schoof, 3: cavity
SPM.mesh.maxb = 1;     %maximum bed gradient
SPM.mesh.maxs = 0.7;     %maximum ice surface gradient
SPM.mesh.hydromode = 3;  %1: steady state, 2: transient Iverson, 3: cavity/channel

SPM.mesh.dohillslope = 1; %hillslope processes
SPM.mesh.dohillslopeerosion = 1; %bedrock erosion on hillslopes
SPM.mesh.dofluvial =1;
SPM.mesh.doperiglacial = 0;
SPM.mesh.Flaw = 4;		   %Fluvial erosion law: 1 : Stream Power law -> E = Ke*Qw^m*Slope^1
						   %			 2 : E = Ke*coverfac*Qw^0.5*S^1
						   %			 3 : E = Ke*coverfac*(tau-tauc)
SPM.mesh.steadybed =0;


%****** ice properties *************
SPM.iprop.ifac = 0.05;       %Relaxation constant (velocity)
SPM.iprop.sfac = 0.05;        %Relaxation constant (stress)
SPM.iprop.vbfac = 0.01;      %Relaxation constant (sliding)
SPM.iprop.maxitt_v = 400;    %Max iteration (velocity)
SPM.iprop.maxitt_s = 100;    %Max iteration (stress)
SPM.iprop.C = 0.25;          %Schoof sliding parameter
SPM.iprop.Cs = 0.5e-2;          %Schoof sliding parameter
SPM.iprop.L0 = 1.0e-4;       %Schoof sliding parameter
SPM.iprop.minbeta = 0.02;    %Schoof sliding parameter
SPM.iprop.gamma0 = 1.0e-2;  %Viscosity regelation

%glacial erosion
SPM.iprop.Kq = 0.0*tfac;        %pre-multiplication factor for quarrying
SPM.iprop.Ka = 1e-4*tfac;         %bedrock abrasion rate
SPM.iprop.ap = 1.0;          %abrasion sliding power
SPM.iprop.minefac = 1.0e-6;  %min abrasion constant

%**** hydrology properties ******
SPM.hwprop.a2w = 0.91;     %fraction of ablation to glacial water 
SPM.hwprop.kh = 0.1506;    %run code htest.m to find these
SPM.hwprop.h0 = 0.1;  
SPM.hwprop.ds = 0.25;
SPM.hwprop.ss0 = 0.3;
SPM.hwprop.kmin = 0.001;    %min transmicity for cavities 
SPM.hwprop.alpha = 0.7;    %cavity area scaling 
SPM.hwprop.Ls = 4;         %Cavity spacing
SPM.hwprop.Lc = 400;       %channel spacing
SPM.hwprop.S0 = 0.1;       %min cavity length
SPM.hwprop.minqw = 1.0e-4; %minimum water flux
SPM.hwprop.dss = 1.0e6*tfac; %sediment carrying capacity
                           
%**** mass balance parameters *********
SPM.mprop.mtype = 2; %linear
SPM.mprop.avaslope = 0.5; %max stable slope for avalances
SPM.mprop.avacurv = -0.0001;  %min stable curvature for avalances (wind)
SPM.mprop.lrate = 0.006;   %lapse rate
SPM.mprop.mPDD = 0.4e-3; %2.0
SPM.mprop.dTa = 10;
SPM.mprop.Ldebris = 0.1;  %ablation reduction beneath debris cover
SPM.mprop.maxacc =1.0;  %maximum accumulation rate
SPM.mprop.scc =0.1;  %maximum accumulation rate
                          
%**** hillslope processes *********
SPM.hprop.Ks = 5*tfac;       %hillslope diffusivity
SPM.hprop.sc = 0.8;        %critical slope
SPM.hprop.Ke = 0.2*tfac;       %Hillslope erosion coefficient
SPM.hprop.kt = 1e-4;       %time coefficient for landslide
SPM.hprop.Nc = 100;         %number of random landslide cells
C = 60e3;               %Hillslope cohesion 6000 kg/m/s (or Pa)
SPM.hprop.gamma = 4*C/(2800*10);       %hillslope landslide strength

%**** Particles properties ********
SPM.parprop.npmax = 1e6; %max number of particles
SPM.parprop.minice = 10.0; %minimum ice thickness for particle formation
SPM.parprop.minsedi = 0.4; %minimum sediment thickness for particle formation
SPM.parprop.maxsedi = 0.5; %max sediment per particle
SPM.parprop.maxpm = 100; %max particle number in margin cell
SPM.parprop.maxp = 200; %max particle number in ordinary cell
SPM.parprop.minpm = 0; %minimum particle number at margin 
SPM.parprop.minsedim = 0; %mimum sediment pickup at margin
                          %The two latter parameters are used to
                          %form sediment along margins without
                          %production of sediment in the landscape

%**** fluvial properties *********
SPM.fprop.pr = 2;        %precipitation rate m/yr
SPM.fprop.rho_s = 2400;    %density of sediment grains
SPM.fprop.Dg = 0.01;       %average grain size
SPM.fprop.kw = 10.0;       %channel width scaling factor
SPM.fprop.tau_c = 0.01;    %critical shear stress
SPM.fprop.Kt = 0.6e-4*tfac;       %transport capacity scaling
SPM.fprop.Ke = 1.1*tfac;       %bedrock erosion rate
SPM.fprop.m = 0.5;         %discharge exponent (Stream Power law)

%***** Flexural isostasy ********
SPM.flexprop.Te = 10.0e3;  %Elastic thickness
SPM.flexprop.rho_r = 2750; %Density of eroded bedrock
SPM.flexprop.rho_s = 2300; %Density of sediment
SPM.flexprop.rho_a = 3250; %Density of asthenosphere
SPM.flexprop.ncall = 1;   %Call isostasy every ncall timestep


%******** periglacial data **********
load ./input/periglacial_input.mat;
SPM.pgprop.nHs = length(Hsv);
SPM.pgprop.nT0 = length(T0v);
SPM.pgprop.Hsv = Hsv;
SPM.pgprop.T0v = T0v;
SPM.pgprop.Ci = Ci;
SPM.pgprop.Tr = Tr;
SPM.pgprop.rho_b = 2700.0;
SPM.pgprop.rho_s = 2400.0;
SPM.pgprop.Ke = 10e-3;
SPM.pgprop.Kt = 1.0;         %Transport scaling
SPM.pgprop.maxsedi = 10.0;   %Max sediment for periglacial erosion
SPM.pgprop.maxice = 10.0;    %Max ice for periglacial erosion
SPM.pgprop.minslope = 0.5;   %minimum slope for sediment free cells
SPM.pgprop.minsedi = 1.0;    %expected minimum sediment cover on
                           %slopes smaller than minslope
                 
%Sea level temperature through time
time = linspace(0,SPM.mesh.maxtime,100); %time vector in years
T0 = zeros(size(time)) + 9.5;
figure
plot(time,T0)

SPM.mprop.Temp = [time(:)';T0(:)'];   %t;T data

bed = zeros(SPM.mesh.ny,SPM.mesh.nx);
xshift = 1;
yshift = 1;
for i=1+yshift:ny-yshift
    for j=1+xshift:nx-xshift
          bed(i,j) = 50;
    end
end

%make fractal noise for the critical slope
x=linspace(0,L,nx);
y=linspace(0,H,ny);
nu=0.8; %roughness
Lx=1000; %wave length
seed=1; %same random numbers
karmanM=karman2d(x,y,nu,Lx,seed);
karmanM=2*(karmanM-min(karmanM(:)))/(max(karmanM(:))-min(karmanM(:)))-1;%normilized Von Karman
phi = SPM.hprop.sc + 0.2*karmanM;
bed = bed +phi;
bed(:,1) = 0;
bed(1,:) = 0;
bed(:,end) = 0;
bed(end,:) = 0;
%Mesh grid
x = linspace(0,SPM.mesh.L,SPM.mesh.nx);
y = linspace(0,SPM.mesh.H,SPM.mesh.ny);
[X,Y] = meshgrid(x,y);

%load mask to limit accumulation
include = zeros(size(bed));
facy = 20;
facx = 10;
%save input data on grid 
SPM.data.bed = bed;
SPM.data.fixflag = zeros(size(bed));
SPM.data.fixflag(:,1) = 1;
SPM.data.fixflag(1,:) = 1;
SPM.data.fixflag(:,end) = 1;
SPM.data.fixflag(end,:) = 1;
SPM.data.ice = zeros(size(bed));
SPM.data.sedi = zeros(size(bed));
SPM.data.precip = ones(size(bed));
SPM.data.include = ones(size(bed));
SPM.data.phi = 0.5*ones(SPM.mesh.ny,SPM.mesh.nx);
SPM.data.Kq = SPM.iprop.Kq.*ones(size(bed));
SPM.data.Ka = SPM.iprop.Ka.*ones(size(bed));
SPM.data.mrate = ones(size(bed));
SPM.data.srate = zeros(size(bed));

for i=1+yshift:ny-yshift
    for j=1+xshift:nx-xshift
          SPM.data.srate(i,j) = 0.1e-3*tfac;
    end
end

%mass balance function (positive degree day)
[SPM,h,h2] = massbalance(SPM,1);

%plot hypsometry with snowlines
elevation = linspace(-500,5000,200);
count = histcounts(bed(bed>0),elevation); figure
plot([count,0],elevation,'k'); hold on;
plot([0,max(count)],[h, h],'b');
plot([0,max(count)],[h2, h2],'b--');
hold off
xlabel('Count'), ylabel('Elevation (m)')

%plot input topography
figure
surf(X,Y,SPM.data.bed); view(0,90),shading interp
set(gca,'dataaspectratio',[1,1,1]);
figure
pcolor(X,Y,SPM.data.include), shading flat
set(gca,'dataaspectratio',[1,1,1]);

%write input files
write(SPM);

