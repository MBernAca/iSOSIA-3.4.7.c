% Check river profile
close all
SPM=SPMload;
nx = SPM.mesh.nx; ny = SPM.mesh.ny;
st = load('./status.dat');
latestfnr = st(3);
fnr = latestfnr;
loaddata

L = SPM.mesh.L;
H = SPM.mesh.H;
dx = SPM.mesh.dx;
dy = SPM.mesh.dy;
include = SPM.data.include;
[X,Y] = meshgrid([1:nx-2]*dx-dx/2,[1:ny-2]*dy-dy/2);

% keep DEM in include area
bed = bed(2:end-1,2:end-1);

DEM = GRIDobj(xc(2:end-1,2:end-1),yc(2:end-1,2:end-1),bed);
figure
imagesc(DEM)
set(gca,'DataAspectRatio',[1,1,1])

% Flow accumulation
FD = FLOWobj(DEM,"preprocess",'fill');
A = flowacc(FD);
figure
subplot(1,2,1)
imageschs(DEM,dilate(sqrt(A),ones(2)),'colormap','flowcolor','caxis',[0 50])

% Stream
S = STREAMobj(FD,A>0.5e2);
S2 = streamorder(FD,A>0.5e2);
S3 = klargestconncomps(S,1);
St = trunk(S3);
s4 = streamorder(S3);
z = imposemin(S,DEM);
s = streamorder(S);

%Drainage basins
D = drainagebasins(FD,S2,3);
subplot(1,2,2)
imageschs(DEM,shufflelabel(D))

figure
subplot(1,2,1)
imageschs(DEM,[],'colormap',[1 1 1]);
hold on
plotc(S3,DEM);
colormap(jet);
hold off

subplot(1,2,2)
plotdz(S3,DEM,'color',s4,'LineWidth',1,'cbarlabel','Stream order')
