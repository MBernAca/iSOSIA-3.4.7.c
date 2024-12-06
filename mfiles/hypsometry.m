function hypsometry(varargin)


% hypsometry - shows hypsometry for SPM model
%
%   DLE 11/3 2013
%
% last update: 18/04/2022 by Maxime

% set(gcf,'units','centimeters','position',[10,10,20,25]);
set(gcf,'color','w');

%load SPM structure
SPM = SPMload;
nx = SPM.mesh.nx;
ny = SPM.mesh.ny;
Nc = length(SPM.data.include(SPM.data.include==1));
% load('../../Cli03_bedSlope/bslopeInit.mat');
% bslope0 = load('./bslopeInit.mat');

%default values
fnr = 0;
Nh = 50;
doclose = 1;
dlim = 0;
showSlope = 0;
dorelief = 0;
slopeElev = 0;
showErosion = 0;
nc = 500;
minbed = -1000;
maxbed = 2000;
normalized = 0;
GOT = 0;
hypso = 0;
reliefElev = 0;
reliefGOT = 0;
fig = 1;
use_init_topo = 0;

%Input variables
if nargin > 0,   
  lArgin = varargin;
  while length(lArgin) >= 1,
    type = lArgin{1};
    lArgin = lArgin(2:end);
    switch type,
     case{'file'}
      fnr = lArgin{1};
      lArgin = lArgin(2:end);
	 case{'dlim'}
	  dlim = 1;
      drange = lArgin{1};
      lArgin = lArgin(2:end);
	 case{'hypsometry'}
		hypso = 1;
	 case{'noclose'}
      doclose = 0;
	 case{'slope'};
	   showSlope = 1;
	 case{'slopeelev'};
	   slopeElev = 1;
	 case{'relief'};
	    dorelief = 1;
	 case{'radius'};
		nc = lArgin{1};
		lArgin = lArgin(2:end);
	 case{'erosion'};
		showErosion = 1;
	 case{'normalize'};
		normalized = 1;
	 case{'GOT'};
	    GOT = 1;
	 case{'nofig'};
	    fig = 0;
     case{'reliefElev'};
	    reliefElev = 1;
	 case{'reliefGOT'};
	    reliefGOT = 1;
	 case{'use_init_topo'};
		use_init_topo = 1;
	 case{'Nh'}
	    Nh = lArgin{1};
		lArgin = lArgin(2:end);
    end;
  end;
end;

%Determine file number
if fnr ~= 0,
  st = load('./status.dat');
  latestfnr = st(3);
  if fnr == 'latest',
    fnr = latestfnr;
    disp(['fnr = ',num2str(fnr)]);
  elseif fnr ,
    fnr = fnr;
  else,
    ['Invalid file number']
    return
  end;
end;
include = SPM.data.include;
fnr_stored = fnr;
fnr = 1;
loaddata;
% bslope0 = bslope0(include==1);
bslope0 = bslope(include==1);
Elev0 = bed(include==1);
fnr = fnr_stored;
%initial bed elevation
bed0 = SPM.data.bed;
abrasion0 = SPM.data.abrasion;

if doclose, 
    close all; 
else
	if fig
	disp('yes')
		figure
	end
end 

if fnr == 0,
  
  bed = bed0;
  dx = 60e3/128; dy = 120e3/256;
  x = dx/2:dx:dx*128;
  y = dy/2:dy:dy*256;
  [xc,yc] = meshgrid(x,y);
  % relief0 = localrelief(bed0,nc);
  % relief = localrelief(bed,nc);
  DEM = GRIDobj(xc,yc,bed);
  relief = localtopography(DEM,nc);
  relief = flipud(relief.Z);
  relief0 = relief;
  Meanbed = localtopography(DEM,10000,'type','mean');
  Meanbed = flipud(Meanbed.Z);
  abrasion = SPM.data.abrasion;
  got = SPM.data.abrasion;
  Meanbed0 = Meanbed;

else,

  %load data
  fid = fopen(['./output/output',num2str(fnr),'.dat']);
  xc = fread(fid,[ny,nx],'double');
  yc = fread(fid,[ny,nx],'double');
  bed = fread(fid,[ny,nx],'double');
  dx = 60e3/128; dy = 120e3/256;
  x = dx/2:dx:dx*128;
  y = dy/2:dy:dy*256;
  [xc,yc] = meshgrid(x,y);
  % relief0 = localrelief(bed0,nc);
  % relief = localrelief(bed,nc);
  DEM = GRIDobj(xc,yc,bed);
  DEM0 = GRIDobj(xc,yc,bed0);
  relief = localtopography(DEM,nc);
  relief = flipud(relief.Z);
  relief0 = localtopography(DEM0,nc);
  relief0 = flipud(relief0.Z);
  relief0 = relief0(include==1);
  Meanbed = localtopography(DEM,10000,'type','mean');
  Meanbed = flipud(Meanbed.Z);
  Meanbed0 = localtopography(DEM0,10000,'type','mean');
  Meanbed0 = flipud(Meanbed0.Z);
  fclose(fid);
  fnr_saved = fnr;
  loaddata;
  
  % handle 2-stage modelling Cli03_kaPrec4c
  got_temp = got;
  loaddata;
  got_temp = got_temp + got;
  fnr = fnr_saved;
  loaddata;
  got = got_temp;

end;

bed = bed(include==1);
bed0 = bed0(include==1);
relief = relief(include==1);
abrasion = abrasion(include==1);
abrasion0 = abrasion0(include==1);
fluvial = fluvial(include==1);
hillslope = hillslope(include==1);
got = got(include==1);
Meanbed = Meanbed(include==1);
Meanbed0 = Meanbed0(include==1);
bedTrue = bed;

if use_init_topo
	bed = bed0;
end
if normalized,
	minbed0 = min(bed0(:));
	maxbed0 = max(bed0(:));
	bed0 = (bed0-minbed0)./(maxbed0-minbed0);
	minbed = min(bed(:));
	maxbed = max(bed(:));
	bed = bed - minbed;
	bed = (bed-min(bed(:)))./(max(bed(:)));
end

if dlim, 
	minele = drange(1);
	maxele = drange(2);
else,
	maxele = max(bed(:));
	minele = min(bed(:));
end

if hypso,
	minele = minbed;%min([bslope(include==1);bslope0(:)]);
	maxele = maxbed,%max([bslope(include==1);bslope0(:)]);
	edges = linspace(minele,maxele,Nh+1);
	N = histc(bed(:),edges).*dx.*dy.*1e-6;
	set(gca,'xlim',[0,500],'ylim',[minele,maxele]);

elseif showSlope,
	minele = minbed;%min([bslope(include==1);bslope0(:)]);
	maxele = maxbed,%max([bslope(include==1);bslope0(:)]);
	edges = linspace(minele,maxele,Nh+1);
	N0 = histc(bslope0(:),edges)/Nc;
	N = histc(bslope(include==1),edges)/Nc;
	set(gca,'xlim',[0,.15],'ylim',[0,1]);
elseif showErosion,
	minele = minbed;%min([bslope(include==1);bslope0(:)]);
	maxele = maxbed,%max([bslope(include==1);bslope0(:)]);
	edges = linspace(minele,maxele,Nh+1);
	[~,~,elev_bin] = histcounts(bed(:),edges);
	Ero = zeros(size(edges)); EroStd = zeros(size(edges));
	Flu = zeros(size(edges)); FluStd = zeros(size(edges));
	Gla = zeros(size(edges)); GlaStd = zeros(size(edges));
	Hill = zeros(size(edges)); HillStd = zeros(size(edges));
	gotM = zeros(size(edges)); gotStd = zeros(size(edges));
	erosion = fluvial+abrasion+hillslope;
	for i=1:length(edges)
		index = find(elev_bin==i);
		Ero(i) = mean(erosion(index));
		EroStd(i) = std(erosion(index));
		Flu(i) = mean(fluvial(index));
		FluStd(i) = std(fluvial(index));
		Gla(i) = mean(abrasion(index));
		GlaStd(i) = std(abrasion(index));
		Hill(i) = mean(hillslope(index));
		HillStd(i) = std(hillslope(index));
		gotM(i) = mean(got(index));
		gotStd(i) = std(got(index));
	end
	
elseif GOT,
	minele = minele;
	maxele = maxele;
	edges = linspace(minele,maxele,Nh+1);
	[~,~,elev_bin] = histcounts(bed(:),edges);
	slopes = bslope(include==1);
	N = zeros(size(edges));
	Ero = zeros(size(edges));
	Ncs = 0;
	for i=1:length(edges)
		Relc = 200;
		index = find(elev_bin==i);
		% count slopes < 10° for the bin elevation
		slope_found = slopes(index);
		relief_found = relief(index);
		index_slope = find(slope_found < 0.177 & relief_found< Relc & bedTrue(index) > Meanbed(index));
		N(i) = length(index_slope);
		% Quantify erosion
		temp2 = got(index); 
		Ero(i) = mean(temp2);
	end
	N = N;
	EroGla2 = Ero./1e3;
	%set(gca,'xlim',[minele,maxele]);
	
elseif slopeElev %Find slope < 10° along elevation
	minele = minele;
	maxele = maxele;
	edges = linspace(minele,maxele,Nh+1);
	[~,~,elev_bin] = histcounts(bed(:),edges);
	[~,~,elev_bin0] = histcounts(bed0(:),edges);
	slopes = bslope(include==1);
	slopes0 = bslope0(:);
	N = zeros(size(edges));
	N0 = zeros(size(edges));
	Ne = zeros(size(edges)); %for elevation
	EroFlu = zeros(size(edges));
	EroGla = zeros(size(edges));
	EroHill = zeros(size(edges));
	N_flu = zeros(size(edges));
	N_gla = zeros(size(edges));
	N_hill = zeros(size(edges));
	% Ncs = [];
	Ncs = 0;
	Nc0s = 0;
	for i=1:length(edges)
		Relc = 200;
		index = find(elev_bin==i);
		% count slopes < 10° for the bin elevation
		slope_found = slopes(index);
		relief_found = relief(index);
		index_slope = find(slope_found < 0.177 & relief_found< Relc & bedTrue(index) > Meanbed(index));
		N(i) = length(index_slope)*dx*dy*1e-6;
		Ncs = Ncs + length(index_slope);
		Ne(i) = length(index)*dx*dy*1e-6;
		% Quantify erosion
		temp = fluvial(index); temp = temp(index_slope);
		temp2 = abrasion(index); temp2 = temp2(index_slope);
		temp3 = hillslope(index); temp3 = temp3(index_slope);
        total_erosion = sum(temp)+sum(temp2)+sum(temp3);
        N_flu(i) = sum(temp)/total_erosion; % contribution of fluvial erosion
        N_gla(i) = sum(temp2)/total_erosion;  % contribution of glacial erosion
        N_hill(i) = sum(temp3)/total_erosion;  % contribution of hillslope erosion
		EroFlu(i) = sum(temp);%./(sum(temp)+sum(temp2));
		EroGla(i) = sum(temp2);%./(sum(temp)+sum(temp2));;
		EroHill(i) = sum(temp3);
		%Initial topo
		index = find(elev_bin0==i);
		% count slopes < 10° for the bin elevation for initial topo
		slope_found = slopes0(index);
		relief_found = relief0(index);
		index_slope = find(slope_found < 0.177 & relief_found< Relc & bed(index) > Meanbed0(index));
		N0(i) = length(index)*dx*dy*1e-6;
		Nc0s = Nc0s + length(index_slope);
		
	end

	EroFlu2 = EroFlu./sum(EroFlu+EroGla+EroHill);
	EroGla2 = EroGla./sum(EroFlu+EroGla+EroHill);
	EroHill2 = EroHill./sum(EroFlu+EroGla+EroHill);
	% set(gca,'xlim',[-0.15,0.15],'ylim',[minele,maxele]);
elseif dorelief
	edges = linspace(minele,maxele,Nh+1);
	N0 = histc(relief0(:),edges)/Nc;
	N = histc(relief(:),edges)/Nc;
	set(gca,'xlim',[minele./1e3,maxele./1e3],'ylim',[0,0.25]);
elseif reliefElev | reliefGOT,
	minele = minele;%min([bslope(include==1);bslope0(:)]);
	maxele = maxele,%max([bslope(include==1);bslope0(:)]);
	edges = linspace(minele,maxele,Nh+1);
	[~,~,elev_bin] = histcounts(bed(:),edges);
	[~,~,elev_bin0] = histcounts(bed0(:),edges);
	Rel = zeros(size(edges)); relStd = zeros(size(edges));
	Rel0 = zeros(size(edges)); rel0Std = zeros(size(edges));
	gotM = zeros(size(edges)); gotStd = zeros(size(edges));
	for i=1:length(edges)
		index = find(elev_bin==i);
		index0 = find(elev_bin0==i);
		Rel(i) = mean(relief(index));
		relStd(i) = std(relief(index));
		Rel0(i) = mean(relief0(index0));
		rel0Std(i) = std(relief0(index0));
		gotM(i) = mean(got(index));
		gotStd(i) = std(got(index));
	end
else,
	edges = linspace(minele,maxele,Nh+1);
	N0 = histc(bed0(:),edges)/Nc;
	N = histc(bed(:),edges)/Nc;
	set(gca,'xlim',[0,.1],'ylim',[minele,maxele]);
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Ploting part %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on; box on;
if hypso,
	for i=1:Nh,
	 h2 = line([0,N(i),N(i),0],[edges(i),edges(i),edges(i+1),edges(i+1)],'color',[0.3,0.3,0.3],'linewidth',1);
	end;
	xlabel('Area (km²)'), ylabel('Elevation (m)')
	axis square
	
elseif dorelief,
	edges = edges./1e3;
	for i=1:Nh,
	  patch([edges(i),edges(i),edges(i+1),edges(i+1)],[0,N(i),N(i),0],[178,255,102]./255);
	end;
	for i=1:Nh,
	 line([edges(i),edges(i),edges(i+1),edges(i+1)],[0,N0(i),N0(i),0],'color','k','linewidth',1.5);
	end;
	xlabel('Relief (km)')
	ylabel('Frequency')
	set(gca,'fontsize',9,'fontname','arial')
	
elseif showSlope,	
	for i=1:Nh,
	  patch([0,N(i),N(i),0],[edges(i),edges(i),edges(i+1),edges(i+1)],[1,.8,0]);
	end;
	for i=1:Nh,
	 line([0,N0(i),N0(i),0],[edges(i),edges(i),edges(i+1),edges(i+1)],'color','b','linewidth',1.5);
	end;

	xlabel('Frequency')
	ylabel('bed gradient')
	
elseif showErosion,	
	hold off
	t = tiledlayout(1,1);
	ax1 = axes(t);
	hold on
	errorbar(Ero./1e3,edges./1e3,EroStd./1e3,'horizontal','o','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor',[0.3,0.3,0.3],'Color',[0.3,0.3,0.3])
	axis square
	xlabel('Erosion (km)'),xlim([0 1.5])
	ylabel('Elevation (m)')
	errorbar(Gla./1e3,edges./1e3,GlaStd./1e3,'horizontal','o','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor',[102,178,255]./255,'Color',[102,178,255]./255)
	axis square
	xlabel('Erosion (km)'),xlim([0 1.5])
	ylabel('Elevation (m)')
	errorbar(Flu./1e3,edges./1e3,FluStd./1e3,'horizontal','o','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor',[102,255,102]./255,'Color',[102,255,102]./255)
	axis square
	xlabel('Erosion (km)'),xlim([0 1.5]),ylim([-1 2])
	ylabel('Elevation (m)')
	errorbar(Hill./1e3,edges./1e3,HillStd./1e3,'horizontal','o','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor',[1,153/255,51/255],'Color',[1,153/255,51/255])
	axis square
	xlabel('Erosion (km)'),xlim([0 1.5]),ylim([-1 2])
	ylabel('Elevation (m)')
	hold off
	set(ax1,'fontsize',9,'fontname','arial')
	legend(ax1,'Total erosion','Glacial erosion','Fluvial erosion','Hillslope erosion','location','northeast')

elseif GOT,
	minele = 0;%min([bed(:);bed0(:)]);
	maxele = 2e3;%max([bed(:);bed0(:)]);
	edges = linspace(minele,maxele,50+1);
	max(EroGla2(:))
	[~,~,elev_bin] = histcounts(EroGla2(:),edges);
	N2 = zeros(size(elev_bin));
	NC = 0;
	for i=1:length(edges)
		index = find(elev_bin==i);
		N2(i) = sum(N(index))./sum(N(:));
	end

	for i=1:50,
	  patch([edges(i),edges(i),edges(i+1),edges(i+1)],[0,N2(i),N2(i),0],[135/255,206/255,1]);
	end;
	% plot(EroGla2,N,'ko')
	xlabel('Mean glacial occupation time (kyr)')
	ylabel('LRS frequency')
	
elseif slopeElev
	t = tiledlayout(1,1);
	ax1 = axes(t);
	for i=1:Nh,
        Flu = N_flu(i)*N(i);
        Gla = N_gla(i)*N(i);
        hill = N_hill(i)*N(i);
        hf = patch(ax1,[0,Flu,Flu,0],[edges(i),edges(i),edges(i+1),edges(i+1)],[188,240,99]./255,'EdgeColor',[0.3,0.3,0.3],'LineWidth',0.5);
        hg = patch(ax1,[Flu,Flu+Gla,Flu+Gla,Flu],[edges(i),edges(i),edges(i+1),edges(i+1)],[102,178,255]./255,'EdgeColor',[0.3,0.3,0.3],'LineWidth',0.5);
        hh = patch(ax1,[Flu+Gla,Flu+Gla+hill,Flu+Gla+hill,Flu+Gla],[edges(i),edges(i),edges(i+1),edges(i+1)],[255,178,102]./255,'EdgeColor',[0.3,0.3,0.3],'LineWidth',0.5);
    end;
	% Surfaces area
    for i=1:Nh,
        he = patch(ax1,[0,-Ne(i)./500.*100,-Ne(i)./500.*100,0],[edges(i),edges(i),edges(i+1),edges(i+1)],[0.6,0.6,0.6], 'EdgeColor',[0.3,0.3,0.3],'linewidth',0.5);
    end;
	for i=1:Nh,
        hh = patch(ax1,[0,-N0(i)./500.*100,-N0(i)./500.*100,0],[edges(i),edges(i),edges(i+1),edges(i+1)],[0.6,0.6,0.6],'FaceColor','none', 'EdgeColor',[0,0,0],'linewidth',0.5);
    end;
	axis square
    xlim([-100 100])
    ax1.XTick = [-100, -50, 0, 50, 100];
    xlabel('Area (km²)')
    ax1.XTickLabel = {'500','250','0','50','100'};
    ylabel('Elevation (km)')
    set(ax1,'fontsize',8,'fontname','arial')
    hold off
elseif reliefElev,	
	t = tiledlayout(1,1);
	ax1 = axes(t);
	hold on;
	e1=errorbar(Rel0./1e3,edges,rel0Std./1e3,'horizontal','o','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor',[0.4,0.4,0.4],'Color','k');
	e2=errorbar(Rel./1e3,edges,relStd./1e3,'horizontal','o','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor','r','Color','r');
	hold off
	axis square
	xlabel('relief (km)'),xlim([0 2])
	if normalized,
		ylabel('Normalized elevation')
	else
		ylabel('Elevation (m)')
	end
	% ax2 = axes(t);
	% e3=errorbar(gotM./1e6,edges./1e3,gotStd./1e6,'horizontal','o','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor',[51/255,153/255,255/255],'Color',[51/255,153/255,255/255]);
	% ax2.XAxisLocation = 'top';ax2.YAxisLocation = 'right';
	% ax2.YTick = [];
	% ax2.Color = 'none'; ax1.Box = 'off'; ax2.Box = 'off';
	% axis square
	% xlabel('Glacial occupation time (Myr)')
	set(ax1,'fontsize',9,'fontname','arial')
	% set(ax2,'fontsize',9,'fontname','arial')
	legend([e1,e2],'init. relief','final relief','got')
	ylim([minele maxele])
	xlim([0 2.5])
elseif reliefGOT,	
    yneg = relStd./1e3; ypos = relStd./1e3;
	xneg = gotStd./1e6; xpos = gotStd./1e6;
	errorbar(gotM./1e6,Rel./1e3,yneg,ypos,xneg,xpos,'o','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor','r','Color','r');
	axis square
	xlabel('Glacial occupation time (Myr)')
	ylabel('Relief (km)')
	set(gca,'fontsize',9,'fontname','arial')
	
else
	for i=1:Nh,
	  patch([0,N(i),N(i),0],[edges(i),edges(i),edges(i+1),edges(i+1)],[1,.8,0]);
	end;
	for i=1:Nh,
	 line([0,N0(i),N0(i),0],[edges(i),edges(i),edges(i+1),edges(i+1)],'color','b','linewidth',1.5);
	end;

	xlabel('Frequency')
	ylabel('Elevation (m)')
	legend({'Present','Initial'})
end




% print -dtiff -r100 hypso60.tif