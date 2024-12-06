function show(varargin)

%
% show function for iSOSIA 3.4.7b
%
% DLE 21/2 2013
%


%load SPM structure
SPM = SPMload;
nx = SPM.mesh.nx;
ny = SPM.mesh.ny;
dx = SPM.mesh.dx;
dy = SPM.mesh.dy;
cell_area = dx*dy;
include = SPM.data.include;

%default values
fnr = 0;
showtime = 0;
contours = 0;
datatype = ['bed'];
ncon = 50;
withice = 0;
withsedi = 0;
showparticles = 0;
dlim = 0;
tlim = 0;
jpg = 0;
inview = [0,90];
hx = 1;
doclose = 1;
lowdata = -1;
cbar = 1;
ELA = 0; %Maxime
normalized = 0; %Maxime - normalized any data to its maximal value
loga = 0;
nc = 4;
combine = 0; %Combined tseries from two iSOSIA models
noforeland = 0;
mask =0;
mask2 = 0;
mask3 = 0;
Relc = 200;
erosionRate = 0;
anim = 0;
nofig = 0;
zoomin = 0;
hp =0;
localTopo = 'mean';
Radius = 10000;
cmapId = 'jet';
reversed = 0;
higherResolution = 0;
useTTB = 0;
EroRatio = 0;
show_ero_ratio = 0;
with_LRS = 0;
zlight = [100,0,500];
show_average = 0;
use_pcolor = 0;
drange = [0 100];
datatypep = 'bed';
datatypep2 = 'bed';
relief_cirque_threshold = 500;
abrasion_threshold = 5;
curvature_threshold = 0.001;

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
     case{'time'}
      showtime = 1;
     case{'data'}
      datatype = lArgin{1};
      lArgin = lArgin(2:end);
     case{'withLRS'}
         with_LRS = 1;
     case{'contours'}
      contours = 1;
     case{'ncon'}
      ncon = lArgin{1};
      lArgin = lArgin(2:end);
     case{'dlim'}
      dlim = 1;
      drange = lArgin{1};
      lArgin = lArgin(2:end);
     case{'tlim'}
      tlim = 1;
      trange = lArgin{1};
      lArgin = lArgin(2:end);
	 case{'cmap'}
      cmapId = lArgin{1};
      lArgin = lArgin(2:end);
	 case{'reverse'}
      reversed = 1;
     case{'mindata'}
      lowdata = lArgin{1};
      lArgin = lArgin(2:end);
     case{'TTB'}
      useTTB = 1;
     case{'pcolor'}
      use_pcolor = 1;
     case{'view'}
      inview = lArgin{1};
      lArgin = lArgin(2:end);
     case{'hx'}
      hx = lArgin{1};
      lArgin = lArgin(2:end);
     case{'withice'}
      withice = 1;
     case{'withsedi'}
      withsedi = 1;
	 case{'nc'}
      nc = lArgin{1};
      lArgin = lArgin(2:end);
     case{'relc'}
      Relc = lArgin{1};
      lArgin = lArgin(2:end);
	 case{'combine'}
	   combine = 1;
	   Models = lArgin{1};
	   lArgin = lArgin(2:end);
     case{'eroRatio'}
	   EroRatio = 1;
     case{'zlight'}
	   zlight =  lArgin{1};
      lArgin = lArgin(2:end);
     case{'particles'}
      showparticles = 1;
	  datatypep = lArgin{1};
      lArgin = lArgin(2:end);
	 case{'ELA'}
      ELA = 1;
	 case{'localTopo'}
	   localTopo = lArgin{1};
	   lArgin = lArgin(2:end);
	 case{'radius'}
	   Radius = lArgin{1};
	   lArgin = lArgin(2:end);
     case{'interpolate'}
       higherResolution = 1;
	   interpFactor = lArgin{1};
	   lArgin = lArgin(2:end);
	 case{'normalize'}
      normalized = 1;
	 case{'noforeland'}
		noforeland = 1;
     case{'nocbar'}
      cbar = 0;
     case{'noclose'}
      doclose = 0;
	 case{'mask'}
	  mask = 1;
	  Sc = lArgin{1};
	  lArgin = lArgin(2:end);
	  datatypep = lArgin{1};
	  lArgin = lArgin(2:end);
	 case{'mask2'}
	  mask2 = 1;
	  Sc = lArgin{1};
	  lArgin = lArgin(2:end);
     case{'mask3'}
	  mask3 = 1;
	  Sc = lArgin{1};
	  lArgin = lArgin(2:end);
	 case{'rate'}
	  erosionRate = 1;
	  fnr_array = lArgin{1};
	  lArgin = lArgin(2:end);
     case{'zoom'}
      zoomin = 1;
	  zoomfac = lArgin{1};
	  lArgin = lArgin(2:end);
     case{'show_ratio'}
      show_ero_ratio = 1;
     case{'jpg'}
      jpg = 1;
	 case{'log'}
      loga = 1;
	 case{'anim'}
		anim = 1;
	  case{'nofig'}
		nofig = 1;
    end;
  end;
end;

%Determine file number
if length(fnr) > 1 % Then we consider we want to plot average value
    fnr_array = fnr;
    fnr = fnr_array(1);
    show_average = 1;
else
    if fnr ~= 0,
      st = load('./status.dat');
      latestfnr = st(3);
      if fnr == 'latest',
        fnr = latestfnr;
        disp(['fnr = ',num2str(fnr)]);
      else
        fnr = fnr;
      end;
    end;
end

if doclose, 
    close all; 
elseif anim == 0 & nofig == 0
    figure;
	set(gcf,'units','normalized','outerposition',[0 0 1 1])
end;


%Time series
if showtime,
  
  %load data
  if combine,
	dim = 1;
	tdata1 = load(['../',Models,'/output/tseries.dat']);
	tfac = 100;
	tdata2 = load(['./output/tseries.dat']);
	tdata2(:,1) = tdata1(end,1).*tfac + tdata2(:,1);
	time = cat(dim,tdata1(:,1).*tfac,tdata2(:,1));
	tdata = zeros(size(time,1),40);
	tdata(:,1) = time;
	for i=2:size(tdata1,2)
		if i == 19 || i == 20 || i == 21 || i == 28 || i == 36 || i == 35 || i == 33 || i == 37 ||i == 38 || i == 39 || i == 40
			tdata(:,i) = cat(dim,tdata1(:,i)./tfac,tdata2(:,i));
		else
			tdata(:,i) = cat(dim,tdata1(:,i),tdata2(:,i));
		end
		
	end
  else
	tdata = load('./output/tseries.dat');
  end
  
  subplot(3,3,1);
  hold on; grid on; box on;
  plot(tdata(:,1),tdata(:,3).*cell_area.*1e-6,'-k');
  %plot(tdata(:,1),tdata(:,15),'-r');
  %plot(tdata(:,1),tdata(:,25),'-b');
  plot(tdata(:,1),tdata(:,26).*cell_area,'-b'); % water head hw
  xlabel('time'); ylabel('H'); title('Mean ice Volume (km^3)');
  hold off;
  
  subplot(3,3,2);
  hold on; grid on; box on;
  plot(tdata(:,1),tdata(:,4),'-k'); %Velocity itt
  plot(tdata(:,1),tdata(:,13),'-r'); %Stress itt
  plot(tdata(:,1),tdata(:,24),'-b');
  xlabel('time'); ylabel('Nitt'); title('# of velocity itterations');
  hold off;

  subplot(3,3,3);
  hold on; grid on; box on;
  plot(tdata(:,1),tdata(:,2),'-k');
  plot(tdata(:,1),tdata(:,11),'-b');
  xlabel('time'); ylabel('dt'); title('Time step');
  hold off;

  subplot(3,3,4);
  hold on; grid on; box on;
  plot(tdata(:,1),tdata(:,5),'-k');
  xlabel('time'); ylabel('speed'); title('Max average speed');
  hold off;

  subplot(3,3,5);
  hold on; grid on; box on;
  plot(tdata(:,1),tdata(:,6).*1e3,'-k');
% fluvial_ero = tdata(:,36).*(nx*ny)/length(include(include==1));
% plot(tdata(:,1),(tdata(:,19)+tdata(:,21)+fluvial_ero).*1e3,'-k');
  plot(tdata(:,1),tdata(:,9).*1e3,'-b'); %Quarrying
  %plot(tdata(:,1),tdata(:,15),'-r');
  plot(tdata(:,1),tdata(:,19).*1e3,'-g'); %Hillslope
  plot(tdata(:,1),tdata(:,20).*1e3,'-y'); %landslide
  plot(tdata(:,1),tdata(:,21).*1e3,'-r'); %Abrasion
  plot(tdata(:,1),tdata(:,28).*1e3,'-m'); %Weathering
  if SPM.mesh.doice == 1, 
	plot(tdata(:,1),tdata(:,36).*1e3,'-c');
  else
	plot(tdata(:,1),tdata(:,36).*1e3,'-c');
  end %Fluvial temporaly correct for include 
  xlabel('time'); ylabel('erosion rate (mm/yr)'); title('erosion rate');
  hold off;
  % ylim([0, 1.5])

  
  subplot(3,3,6);
  hold on; grid on; box on;
  %plot(tdata(:,1),tdata(:,22),'-k');
  % plot(tdata(:,1),tdata(:,33),'-r');
  plot(tdata(:,1),tdata(:,35).*1e3,'-b');
  %plot(tdata(:,1),(tdata(:,35)+SPM.data.srate(100,100)).*1e3,'-g');
  xlabel('time'); ylabel('Right BC'); title('uplift rate (mm/yr)');
  hold off;
  ylim([-50 50])

  subplot(3,3,7);
  hold on; grid on; box on;
  plot(tdata(:,1)./1e6,tdata(:,10),'-r');
  xlabel('time (Myr)'); ylabel('Sealevel T'); title('Sealevel T');
  % ylim([3.5 11])
  hold off;

  
  subplot(3,3,8);
  hold on; grid on; box on;
  plot(tdata(:,1),tdata(:,16),'-k');
  plot(tdata(:,1),tdata(:,17),'-b');
  plot(tdata(:,1),tdata(:,18),'-r');
  xlabel('time'); ylabel('Elevation (m)'); 
  title('min, max, and mean elevaion');
  hold off;
%   ylim([-1000, 2000])

  subplot(3,3,9);
  hold on; grid on; box on;
  plot(tdata(:,1),tdata(:,37),'-k'); %left
  plot(tdata(:,1),tdata(:,38),'-g'); %right
  plot(tdata(:,1),tdata(:,39),'-b'); %bottom
  plot(tdata(:,1),tdata(:,40),'-r'); %left
  plot(tdata(:,1),tdata(:,33),'-m'); %left
  xlabel('time'); ylabel('Sediment flux'); 
  title('Sediment');
  hold off;

  
else,
  
  %Retrieve data
  nx = SPM.mesh.nx;
  ny = SPM.mesh.ny;
  L = SPM.mesh.L;
  H = SPM.mesh.H;
  dx = SPM.mesh.dx;
  dy = SPM.mesh.dy;
  include = SPM.data.include;
  fixflag = SPM.data.fixflag;
  bed_Init = SPM.data.bed;
  [Xc,Yc] = meshgrid([1:nx]*dx-dx/2,[1:ny]*dy-dy/2);
  if show_average, % Compute average value
        [Xc,Yc] = meshgrid([1:nx]*dx-dx/2,[1:ny]*dy-dy/2);
        sliding_average = zeros(size(SPM.data.bed));
        ice_average = zeros(size(SPM.data.bed));
        for i=1:length(fnr_array)
            fnr = fnr_array(i);
            loaddata
            sliding_average = sliding_average + sliding;
            ice_average = ice_average+ice;
        end
        sliding = sliding_average./length(fnr_array);
        ice = ice_average./length(fnr_array);
        cirques = zeros(size(SPM.data.bed));
        erosion_prev = zeros(size(SPM.data.bed));
        erosion_rate = zeros(size(SPM.data.bed));
        erosion_diff = zeros(size(SPM.data.bed));
        got_diff = zeros(size(SPM.data.bed));
        fluvial_diff = zeros(size(SPM.data.bed));
        abrasion_diff = zeros(size(SPM.data.bed));
        hillslope_diff = zeros(size(SPM.data.bed));
        bcurv = zeros(size(SPM.data.bed));
        fluvial_rate= zeros(size(SPM.data.bed));
        abrasion_rate = zeros(size(SPM.data.bed));
        hillslope_rate = zeros(size(SPM.data.bed));
        surface_uplift_prev = zeros(size(SPM.data.bed));
        relief = zeros(size(SPM.data.bed));
        Meanbed = zeros(size(SPM.data.bed));

  else   
      if fnr == -1 | fnr == 0
        
        bed = SPM.data.bed;
        ice = SPM.data.ice;
        bslope = get_bed_curvature(SPM.data.bed,dx,dy);
        if abs(abs(dx)-abs(dy))>1e-9
		    x = dx/2:dx:dx*nx;
		    y = dx/2:dx:dx*ny;
		    [Xc,Yc] = meshgrid(x,y);
        else
            [Xc,Yc] = meshgrid([1:nx]*dx-dx/2,[1:ny]*dy-dy/2);
        end
        DEM = GRIDobj(Xc,Yc,bed);
	    %Compute local relief
	    relief = localtopography(DEM,Radius);
        relief = flipud(relief.Z);
	    sediment = SPM.data.sedi;
	    %make mesh for plotting
        abrasion = zeros(size(SPM.data.bed));
	    quarrying = zeros(size(SPM.data.bed));
	    fluvial = zeros(size(SPM.data.bed));
	    hillslope = zeros(size(SPM.data.bed));
	    periglacial = zeros(size(SPM.data.bed));
	    isostasy = zeros(size(SPM.data.bed));
        erosion_rate = zeros(size(SPM.data.bed));
        abrasion_rate = zeros(size(SPM.data.bed));
        fluvial_rate = zeros(size(SPM.data.bed));
        hillslope_rate = zeros(size(SPM.data.bed));
        Tb = zeros(size(SPM.data.bed));
        SLf = zeros(size(SPM.data.bed));
        te = zeros(size(SPM.data.bed));
        ts = zeros(size(SPM.data.bed));
        icemargin = zeros(size(SPM.data.bed));
        bQw = zeros(size(SPM.data.bed));
	    curv = zeros(size(SPM.data.bed));
	    deformation = zeros(size(SPM.data.bed));
	    Meanbed = zeros(size(SPM.data.bed));
	    got = zeros(size(SPM.data.bed));
	    sliding = zeros(size(SPM.data.bed));
        cirques = zeros(size(SPM.data.bed));
        erosion_prev = include.*0;
        surface_uplift_prev = include.*0;
        erosion_diff = include.*0;
        fluvial_diff = include.*0;
        abrasion_diff = include.*0;
        hillslope_diff = include.*0;
        bcurv = include.*0;
        got_diff = include.*0;
      else,
        loaddata;
        %res = sqrt(vxres(:,2:end).^2+vyres(2:end,:).^2);
        %xres = vxres(:,2:end);
        %yres = vyres(2:end,:);
	    %If ELA - computes the mean ELA position and its 1 sigma deviation
	    %ELA is defined where effective mass balance is roughly null
	    if ELA
	      index = find(massb(SPM.data.include==1) < 0.0001 & massb(SPM.data.include==1) > -0.0001 & ice(SPM.data.include==1) > 10);
	      bed_include = bed(SPM.data.include==1); massb_include = massb(SPM.data.include==1); ice_include = ice(SPM.data.include==1);
	      meanELA = mean(bed_include(index)+ice_include(index));
	      stdELA = std(bed_include(index)+ice_include(index));
	      disp(['mean ELA = ',num2str(meanELA), ' m'])
	      disp(['one-sigma std = ', num2str(stdELA), ' m'])
	    end
	    %make mesh for plotting
	    [Xc,Yc] = meshgrid([1:nx]*dx-dx/2,[1:ny]*dy-dy/2);
    
	    %color data
	    icevelo = deformation + sliding; %Mean ice velocity
    
	    %Compute erosion rate
	    fnr_stored = fnr;
	    if erosionRate
		    fnr =  fnr_array
	    elseif fnr < 500 & erosionRate ==0
	        fnr = fnr -1;
	    else
		    fnr = 0;
	    end
	    if fnr <= 0;
		    abrasion = zeros(size(SPM.data.bed));
		    quarrying = zeros(size(SPM.data.bed));
		    fluvial = zeros(size(SPM.data.bed));
		    hillslope = zeros(size(SPM.data.bed));
		    periglacial = zeros(size(SPM.data.bed));
		    isostasy = zeros(size(SPM.data.bed));
		    erosion_rate = zeros(size(SPM.data.bed));
            abrasion_rate = zeros(size(SPM.data.bed));
            fluvial_rate = zeros(size(SPM.data.bed));
		    Meanbed = zeros(size(SPM.data.bed));
		    bed = SPM.data.bed;
            erosion_prev = include.*0;
            bcurv =  include.*0;
	    else 
		    loaddata;
        end
	    abrasion_prev = abrasion;
	    hillslope_prev = hillslope;
	    fluvial_prev = fluvial;
	    erosion_prev = abrasion+quarrying+hillslope+fluvial+periglacial;
        got_prev = got;
	    isostasy_prev = isostasy;
        surface_uplift_prev = (SPM.data.srate(floor(ny/2),floor(nx/2))*(fnr*SPM.mesh.filetime) + isostasy) - erosion_prev;
	    bed_prev = bed;
	    fnr_diff = abs(fnr-fnr_stored);
	    fnr = fnr_stored;
	    loaddata;
	    abrasion_rate = (abrasion-abrasion_prev)./(SPM.mesh.filetime.*fnr_diff).*1e3; %(mm/yr)
	    hillslope_rate = (hillslope-hillslope_prev)./(SPM.mesh.filetime.*fnr_diff).*1e3;%(mm/yr)
	    fluvial_rate = (fluvial-fluvial_prev)./(SPM.mesh.filetime.*fnr_diff).*1e3;%(mm/yr)
        flu_hill_rate = (fluvial+hillslope-fluvial_prev-hillslope_prev)./(SPM.mesh.filetime.*fnr_diff).*1e3;
	    %Evaluate topographic changes
	    Topo_init = bed-isostasy;
	    elevdiff = bed - bed_Init;
	    %Correct for uplift
	    uplift = SPM.data.srate(floor(ny/2),floor(nx/2));
	    filetime = SPM.mesh.filetime;
	    maxtime = filetime*fnr;
	    topo = Topo_init ;%- (uplift*maxtime);
	    uplift_rate = ((isostasy-isostasy_prev).*1e3)./(SPM.mesh.filetime) + SPM.data.srate(floor(ny/2),floor(nx/2));
	    % %Total erosion
	    erosion = abrasion+quarrying+hillslope+fluvial+periglacial;
	    erosion_rate = (erosion-erosion_prev)./(SPM.mesh.filetime.*fnr_diff).*1e3;
        erosion_diff = (erosion-erosion_prev);
        got_diff = (got-got_prev);
        abrasion_diff = (abrasion-abrasion_prev);
        fluvial_diff = (fluvial-fluvial_prev);
        hillslope_diff = (hillslope-hillslope_prev);
	    bQw = bQw*(3600*24*365.25);
	    %Compute local relief
	    if abs(abs(dx)-abs(dy))>1e-9
		    x = dx/2:dx:dx*nx;
		    y = dx/2:dx:dx*ny;
		    [Xc,Yc] = meshgrid(x,y);
	    end
	    DEM = GRIDobj(Xc,Yc,bed);
        bcurv = get_bed_curvature(bed,dx,dy);
	    relief = localtopography(DEM,Radius);
	    relief = flipud(relief.Z);
	    Meanbed = localtopography(DEM,Radius,'type',localTopo);
	    Meanbed = flipud(Meanbed.Z);
        cirques = zeros(size(bed));
        [cirque_density,~,~,index_cirques] = get_cirques_density(bcurv,relief,bed,curvature_threshold,relief_cirque_threshold,abrasion,fluvial,abrasion_threshold,xc,yc,0);
        cirques(index_cirques) = 1;
        cirques = flipud(cirques); % indices are flipped
      end
  end
  
  % Show higher resolution model?
  if higherResolution,
    x = linspace(SPM.mesh.L/(SPM.mesh.nx*interpFactor)/2,SPM.mesh.L,SPM.mesh.nx*interpFactor);
    y = linspace(SPM.mesh.H/(SPM.mesh.ny*interpFactor)/2,SPM.mesh.H,SPM.mesh.ny*interpFactor);
    [Xc,Yc] = meshgrid(x,y);
    massb = interp2(xc,yc,massb,Xc,Yc);
    include = interp2(xc,yc,include,Xc,Yc); include(include>0)= 1;
    ice = interp2(xc,yc,ice,Xc,Yc);
    bed = interp2(xc,yc,bed,Xc,Yc);
  end

  %figure
  set(gca,'position',[0.1,0.05,0.7,1],'visible','off');
  set(gca,'visible','on');
  hold on; grid on; box on; %axis vis3d;
  if exist(['../../spm-3.4.7.b/mfiles/colormap/',cmapId,'.mat'])
        load(['../../spm-3.4.7.b/mfiles/colormap/',cmapId,'.mat']);
  else
    load(['../spm-3.4.7.b/mfiles/colormap/',cmapId,'.mat']);
  end
  colormap(map); %colormap
  caxis([0,5]);
  cmap = colormap;
  
  if noforeland,
	if mask | mask2
        if abs(abs(dx)-abs(dy))>1e-9
		    x = dx/2:dx:dx*nx;
		    y = dx/2:dx:dx*ny;
		    [Xc,Yc] = meshgrid(x,y);
	    end
	    DEM = GRIDobj(Xc,Yc,bed);
		relief = localtopography(DEM,500);
		relief = flipud(relief.Z);
		%relief = localrelief(bed,1);
		relief = relief(include==1);
		Meanbed = localtopography(DEM,Radius,'type','mean');
		Meanbed = flipud(Meanbed.Z);
	else
		relief = relief(include==1);
	end
	nx = 108; ny =236;
	x = linspace(0,468.75*(nx-1),108); 
	y = linspace(0,468.75*(ny-1),236);
	[Xc, Yc] = meshgrid(x,y);
	relief = reshape(relief,ny,nx);
	Meanbed = Meanbed(include==1);
	Meanbed = reshape(Meanbed,ny,nx);
	bed = bed(include==1);
	bed = reshape(bed,ny,nx);
    index = bed>=Meanbed;
    indexBelow = bed<Meanbed;
    BelowMeanBed = Meanbed;BelowMeanBed(index) = 0; BelowMeanBed(indexBelow) = 1; 
    AboveMeanBed = Meanbed;AboveMeanBed(index) = 1; AboveMeanBed(indexBelow) = 0; 
%     qw = qw(include==1);
%     qw = reshape(qw,ny,nx);
	bed = reshape(bed,ny,nx);
	ice = ice(include==1);
    ice = reshape(ice,ny,nx);
    cirques = cirques(include==1);
    cirques = reshape(cirques,ny,nx);
	curv = curv(include==1);
    curv = reshape(curv,ny,nx);
	deformation = deformation(include==1);
    deformation = reshape(deformation,ny,nx);
    Tb = Tb(include==1);
    Tb = reshape(Tb,ny,nx);
	bslope = bslope(include==1);
	bslope = reshape(bslope,ny,nx);
	sediment = sediment(include==1);
	abrasion = abrasion(include==1);
	fluvial = fluvial(include==1);
	hillslope = hillslope(include==1);
	periglacial = periglacial(include==1);
	quarrying = quarrying(include==1);
	sediment = reshape(sediment,ny,nx);
	abrasion = reshape(abrasion,ny,nx);
    hillslope = reshape(hillslope,ny,nx);
    SLf = SLf(include==1);
    SLf = reshape(SLf,ny,nx);
    bQw = bQw(include==1);
    bQw = reshape(bQw,ny,nx);
    te = te(include==1);
    te = reshape(te,ny,nx);
    ts = ts(include==1);
    ts = reshape(ts,ny,nx);
    icemargin = icemargin(include==1);
    icemargin = reshape(icemargin,ny,nx);
	fluvial = reshape(fluvial,ny,nx)+hillslope; % Add hillslope !
	periglacial = reshape(periglacial,ny,nx);
	quarrying = reshape(quarrying,ny,nx);
	erosion = abrasion+fluvial+hillslope+quarrying+periglacial;
    erosion = reshape(erosion,ny,nx);
    erosion_prev = erosion_prev(include==1);
    erosion_prev = reshape(erosion_prev,ny,nx);
    erosion_rate = erosion_rate(include==1);
    erosion_rate = reshape(erosion_rate,ny,nx);
    erosion_diff = erosion_diff(include==1);
    erosion_diff = reshape(erosion_diff,ny,nx);
    hillslope_diff = hillslope_diff(include==1);
    hillslope_diff = reshape(hillslope_diff,ny,nx);
    fluvial_diff = fluvial_diff(include==1);
    fluvial_diff = reshape(fluvial_diff,ny,nx)+hillslope_diff; % Add hillslope !
    abrasion_diff = abrasion_diff(include==1);
    abrasion_diff = reshape(abrasion_diff,ny,nx);
    hillslope_rate = hillslope_rate(include==1);
    hillslope_rate = reshape(hillslope_rate,ny,nx);
    fluvial_rate = fluvial_rate(include==1);
    fluvial_rate = reshape(fluvial_rate,ny,nx)+hillslope_rate; % add hillslope
    abrasion_rate = abrasion_rate(include==1);
    abrasion_rate = reshape(abrasion_rate,ny,nx);
    bcurv = bcurv(include==1);
    bcurv = reshape(bcurv,ny,nx);
    Ratio_EfEg = fluvial_rate./abrasion_rate;
	got = got(include==1);
    got= reshape(got,ny,nx)./1e6;
    got_diff = got_diff(include==1);
    got_diff = reshape(got_diff,ny,nx)./1e6;
	sliding = sliding(include==1);
    sliding= reshape(sliding,ny,nx);
	isostasy = isostasy(include==1);
    isostasy= reshape(isostasy,ny,nx);
    surface_uplift_prev = surface_uplift_prev(include==1);
    surface_uplift_prev = reshape(surface_uplift_prev,ny,nx);
    surface_uplift = ((SPM.data.srate(floor(ny/2),floor(nx/2))*(fnr*SPM.mesh.filetime) + isostasy) - erosion) - surface_uplift_prev;
    flatness_index = erosion./(SPM.data.srate(floor(ny/2),floor(nx/2))*(fnr*SPM.mesh.filetime) + isostasy);
    % Erosion ratio
    if EroRatio
        time149 = (149-87)*SPM.mesh.filetime;
        time87 = 87*SPM.mesh.filetime;
        GlaFlu_ratio = abrasion./(fluvial+1e-6);
        erosion_ratio = (erosion-erosion_prev+1e-6);
        indexes = bed >= Meanbed;
        indexesBelow = bed < Meanbed;
        if show_ero_ratio
            erosion_ratio_Meanbed = (mean(fluvial_rate(indexes) + abrasion_rate(indexes) + 1e-6)./length(indexes)) ./ (mean(fluvial_rate(indexesBelow) + abrasion_rate(indexesBelow) + 1e-6)./length(indexesBelow))
        end
    end
	% Get streams using topotoolbox
	  valleyCenter = zeros(size(bed));
	  DEM = GRIDobj(Xc./1e3,Yc./1e3,bed);
	  FD = FLOWobj(DEM);
	  % calculate flow accumulation
	  A = flowacc(FD);
	  % Note that flowacc returns the number of cells draining
	  % in a cell. Here we choose a minimum drainage area of 10000 cells.
	  W = A>20;
	  % create an instance of STREAMobj
	  S = STREAMobj(FD,W);
	  valleyCenter(S.IXgrid) = 1;
	  valleyCenter = flipud(valleyCenter);
	  reliefCenter = relief;
	  reliefCenter(valleyCenter==0) = 0;
  end
  if tlim, 
      mintopo = trange(1);
      maxtopo = trange(2);
  else,
      mintopo = min(bed(:));
      maxtopo = max(bed(:));
  end;

  eval(['data = ',datatype,';']);
  % data = data(:,2:end);
  %if not topography

    %Show logarythm values
  if loga == 1
	data = log10(data+1e-16);
  end
  
  if ~strcmp(datatype,'bed');
    
      %I = find(abs(data(:)) > lowdata);
		
      if dlim, 
          if numel(drange) == 2
              mindata = drange(1);
              maxdata = drange(2);
              I = find(data(:) > mindata);
          else
            min_pos = drange(3); max_pos = drange(4);
            min_neg = drange(2); max_neg = drange(1);
            index_pos = find(data>=0);
            index_neg = find(data< 0);
          end
	  elseif normalized,
		  mindata = 0; maxdata = 1;
		  data = data./max(data(:));
		  I = find(abs(data(:)) > 0.000001*(max(data(:))-min(data(:))));
      else,
          mindata = min(data(:)); maxdata = max(data(:));
          I = find(abs(data(:)) > 0.000001*(max(data(:))-min(data(:))));
      end; 
      
	  %Create a mask of data (works with dlim)
	  if mask,
	    % low-relief surface = bslope<10° and area of ~870 m2
		index = find(data>Sc);
	    data2 = nan(size(data));
		data2(index) = data(index);
		data = data2;
	  elseif mask2,
		index = find(bslope<Sc & relief < Relc & bed > Meanbed);
	    data2 = zeros(size(data));
		data2(index) = 1;
        data2 = bwmorph(data2,'clean');
        index_LRS = find(data2==1);

      elseif mask3,
        fnr = 87;
        loaddata;
		index = find(bslope<Sc & relief < Relc & bed > Meanbed);
	    data2 = zeros(size(data));
		data2(index) = 1;
        fnr = 149;
        loaddata;
		index = find(bslope<Sc & relief < Relc & bed > Meanbed);
        data2(index) = data2(index) + 1;
	  end
		 
      if numel(drange) == 4
            data2 = zeros(size(data));
            data2(index_pos) = ((data(index_pos)./max_pos) + 1) .* 0.5;
            data2(index_neg) = (1- (abs(data(index_neg)) ./ abs(max_neg))).*0.5;
            data2(data2 < 0) = 0;
            data2(data2 > 1) = 1;
            cdata = data2;
            mindata = drange(1);
            maxdata = drange(4);
      else
        cdata =  (data-mindata)/(maxdata-mindata+1e-16); 
      end
	  
      I = find(cdata(:)>1); cdata(I) = 1;
      I = find(cdata(:) <= 0); cdata(I) = 0;
	  if reversed,
		cdata = 1-cdata;
	  end
	  cdata = cdata + 4;
      zdata = bed;%+sediment+ice;
	  if withice,
		zdata=bed+ice;
	  end
      %make Cdata_M
      %Data
      xi = linspace(0,5,length(cmap(:,1)));
      rc = interp1(xi,cmap(:,1),cdata);
      gc = interp1(xi,cmap(:,2),cdata);
      bc = interp1(xi,cmap(:,3),cdata);
      cdata_M(:,:,1) = rc;
      cdata_M(:,:,2) = gc;
      cdata_M(:,:,3) = bc;
      if mask2,
			rc(index) = 0;
			gc(index) = 0;
			bc(index) = 1;
			cdata_P(:,:,1) = rc;
			cdata_P(:,:,2) = gc;
			cdata_P(:,:,3) = bc;
	  end
      
  else

        %bed topography
        zdata = bed;
        cdata = 0.99*(zdata-mintopo)/(maxtopo-mintopo);
        if normalized,
          if tlim
              maxtopo = maxtopo;
          else
              maxtopo = 1;
          end
          zdata = zdata + abs(min(zdata(:)));
		  zdata = zdata./max(zdata(:));
          cdata = 0.99*(zdata)/(maxtopo);
        end
        I = find(cdata(:) > 0.99); cdata(I) = 0.99;

        %make Cdata_M
        %Elevation
        xi = linspace(0,5,length(cmap(:,1)));
        rc = interp1(xi,cmap(:,1),cdata);
        gc = interp1(xi,cmap(:,2),cdata);
        bc = interp1(xi,cmap(:,3),cdata);
        cdata_M(:,:,1) = rc;
        cdata_M(:,:,2) = gc;
        cdata_M(:,:,3) = bc;
		
        if mask,
			eval(['datap = ',datatypep,';']);
			% find data > sc
			index = find(datap>Sc);
			data2 = zeros(size(bed));
			data2(index) = datap(index);
			datap = data2;
			if dlim, 
			    mindata = drange(1);
			    maxdata = drange(2);
			    Ip = find(datap(:) > mindata);
			else,
				mindata = min(datap(:)); maxdata = max(datap(:));
				Ip = find(abs(datap(:)) > 0.000001*(max(datap(:))-min(datap(:))));
			end; 
			cdatap =  (datap-mindata)/(maxdata-mindata+1e-16); 
			cdatap(cdatap>0) = cdatap(cdatap>0) + 4;
			if reversed,
				cdatap(cdatap>4) = abs(max(cdatap(cdatap>4))-cdatap(cdatap>4))+4;
            end
            cdatap(cdatap>5) = 5;
            cdatap(cdatap==0) = cdata(cdatap==0);
            xi = linspace(0,5,length(cmap(:,1)));
			rc = interp1(xi,cmap(:,1),cdatap);
			gc = interp1(xi,cmap(:,2),cdatap);
			bc = interp1(xi,cmap(:,3),cdatap);
            % Add LRS
            if with_LRS
                index = find(bslope<0.177 & relief < Relc & bed > Meanbed);
			    data2 = zeros(size(bed));
			    data2(index) = 1;
                data2 = bwmorph(data2,'clean');
                data2 = double(data2);
                index = find(data2==1);
                data2(data2==0) = nan;
                rc(index) = 1;
			    gc(index) = 0;
			    bc(index) = 0;
            end
			cdata_P(:,:,1) = rc;
			cdata_P(:,:,2) = gc;
			cdata_P(:,:,3) = bc;
		elseif mask2,
			index = find(bslope<Sc & relief < Relc & bed > Meanbed);
			data2 = zeros(size(bed));
			data2(index) = 1;
            data2 = bwmorph(data2,'clean');
            data2 = double(data2);
            index = find(data2==1);
            data2(data2==0) = nan;
            rc(index) = 1;
			gc(index) = 0;
			bc(index) = 0;
			cdata_P(:,:,1) = rc;
			cdata_P(:,:,2) = gc;
			cdata_P(:,:,3) = bc;
        
        end

        if (withice),
            %Add ice
            ci = 1.1;
            rci = interp1(xi,cmap(:,1),ci);
            gci = interp1(xi,cmap(:,2),ci);
            bci = interp1(xi,cmap(:,3),ci);
            fi = .5*(erf((ice-20)/20)+1);
  
            cdata_M(:,:,1) = (1-fi).*cdata_M(:,:,1)+fi.*rci;
            cdata_M(:,:,2) = (1-fi).*cdata_M(:,:,2)+fi.*gci;
            cdata_M(:,:,3) = (1-fi).*cdata_M(:,:,3)+fi.*bci;
  
            zdata = zdata + ice;       
        
        end;
        
        if (withsedi);
                        
            %surface sediment blending
            cs = 2.7;
            rcs = 0.0;
            gcs = 0.0;
            bcs = 0.0;
            
            mins = 0.05;
            ds = 0.05;
    
            fs = .5*(erf((ssedi-mins)/ds)+1);
            I = find(ssedi < 0.01); fs(I) = 0;
  
            cdata_M(:,:,1) = (1-fs).*cdata_M(:,:,1)+fs.*rcs;
            cdata_M(:,:,2) = (1-fs).*cdata_M(:,:,2)+fs.*gcs;
            cdata_M(:,:,3) = (1-fs).*cdata_M(:,:,3)+fs.*bcs;
  
            zdata = zdata + sediment;
            
        end;
        
    
        
  
  end;

        
  %plot surface
  if anim
   clear clf
   end
  
   % Do I want a cropped image ?
  if useTTB % Use topoToolBox 3D visualization?
    [cmap,zlimits] = ttcmap(DEM,'cmap','mby');
    azi = 160;
    RGB = imageschs(DEM,[],'colormap',cmap);
    C = castshadow(DEM,'azi',azi,'alt',10);
    C = C.Z;
    clr=cmap;
    minval = min(C(:));
    maxval = max(C(:));
    valuesatclr = linspace(minval,...
                       maxval,...
                       size(clr,1))';
    RGBmat = interp1(valuesatclr,clr,C(:)); 
    RGBmat = reshape(RGBmat,DEM.size(1),DEM.size(2),3);
    RGBmat = RGBmat*255;
    
    a2     = (C-minval)./(maxval-minval)*1;
    
    figure('WindowState','maximized')
    MASK = DEM; MASK.Z(isnan(MASK.Z)) = 0; MASK.Z(MASK.Z~=0) = 1;
    DEMc = crop(DEM,MASK,0);
    h = surf(DEMc,'block',true);colormap(cmap);cb=colorbar;
    h.CData = RGB;
    axis off
    exaggerate(gca,1)
    view([-90 90])
    
  else %use iSOSIA visualization
      if use_pcolor
           hp = image(gca,Xc./1e3,Yc./1e3,cdata_M);
           shading interp
      else
          if mask2 | mask,
	        hc=surf(gca,Xc./1e3,Yc./1e3,(zdata+50)./1e3,'Cdata',cdata_P,'facealpha',0.8); shading interp
	        set(hc,'facelighting','gouraud','edgelighting','gouraud');
	        set(gca,'ambientlightcolor',[.9,.8,.9],'Clipping','off');
	        material([.4,.4,.4,2,0.1]);
            %hcc = contour3(Xc./1e3,Yc./1e3,(zdata+data2)./1e3,1,'k')
          else
            hp = surf(gca,Xc./1e3,Yc./1e3,zdata./1e3,'Cdata',cdata_M);
            set(hp,'facelighting','gouraud','edgelighting','gouraud');
            set(gca,'ambientlightcolor',[.9,.8,.9],'Clipping','off');
            material([.4,.4,.4,2,0.1]);
          end
          if strcmp(datatypep,'cirques')
              indices = find(cirques==1);
              scatter3(Xc(indices)./1e3,Yc(indices)./1e3,zdata(indices)./1e3,5,'b','fill')
          end

          %axis equal; 
          set(gca,'dataaspectratio',[1,1,1/hx]);
          shading interp,
          
          view(inview);  
      end
      xlim([0 max(Xc(:))./1e3]);ylim([0 max(Yc(:))./1e3]);

      if zoomin
          set(gca,'Clipping','on')
          % zoomfac is a array of size 7
          % 1: zoom factor, 2-3: xlim, 4-5: ylim, 6-7: zlim
          zoom(zoomfac(1))
          axis([zoomfac(2) zoomfac(3) zoomfac(4) zoomfac(5)])
          zlim([zoomfac(6) zoomfac(7)])
      end
  end
  
  if showparticles,
            
      clear xp;
      clear yp;
      clear zp;
      loadparticles;
      
      eval(['datap = ',datatypep,';']);
	  
	  mindatap = min(datap(:)); maxdatap = max(datap(:));
      
      max(zpr)
      maxdatap
	  disp(['Nb particles =  ', num2str(length(zpr))])
      
      scatter3(xp./1e3,yp./1e3,zpr./1e3,5,datap,'filled');

      
  end;
  
  
  %keyboard
   
  if contours, 
      
      if tlim, vc = linspace(trange(1),trange(2),ncon);
      elseif withice, vc = linspace(min(bed(:)+sediment(:)+ice(:)),max(bed(:)+sediment(:)+ice(:)),ncon);
      else vc = linspace(min(bed(:)+sediment(:)),max(bed(:)+sediment(:)),ncon);
      end;
  
      if withice,
          contour3(Xc./1e3,Yc./1e3,(bed+sediment+ice+10)./1e3,vc./1e3,'color',[0.3,0.3,0.3]); 
      else, 
          contour3(Xc./1e3,Yc./1e3,(bed+1)./1e3,vc./1e3,'color',[0.4,0.4,0.4],'Linewidth',0.5);
      end;
  end;
  
  if anim & fnr<=1
	  if cbar, 
		  vt = [0,0.25,0.5,0.75,1.0];
		  cb = colorbar('position',[0.85,0.1,0.02,0.4]);
		  set(cb,'ylim',[0,1],'ytick',vt,'yticklabel',mintopo+vt*(maxtopo-mintopo));
		  set(get(cb,'ylabel'),'string',['elevation (m)']);
	  end;
	  delete(cb)
	  if ~strcmp(datatype,'bed') & cbar,
		vt = [0,0.25,0.5,0.75,1.0];
		cb = colorbar('position',[0.85,0.55,0.02,0.4]);
		set(cb,'ylim',[4,5],'ytick',vt+4,'yticklabel',mindata+vt*(maxdata-mindata));
		set(get(cb,'ylabel'),'string',datatype);
	  end;
	  if mask,
		  vt = [0,0.25,0.5,0.75,1.0];
		  cb = colorbar('position',[0.85,0.55,0.02,0.4]);
		  set(cb,'ylim',[4,5],'ytick',vt+4,'yticklabel',mindata+vt*(maxdata-mindata));
		  set(get(cb,'ylabel'),'string',datatypep);

	  end
  elseif anim == 0
    if cbar & useTTB==0, 
		  vt = [0,0.25,0.5,0.75,1.0];
		  cb = colorbar('position',[0.85,0.1,0.02,0.4]);
		  set(cb,'ylim',[0,1],'ytick',vt,'yticklabel',mintopo+vt*(maxtopo-mintopo));
		  set(get(cb,'ylabel'),'string',['elevation (m)']);
	  end;

	  if ~strcmp(datatype,'bed') && cbar == 1 | mask & useTTB==0,
		delete(cb)
		vt = [0,0.25,0.5,0.75,1.0];
		cb = colorbar('position',[0.85,0.55,0.02,0.4]);
		set(cb,'ylim',[4,5],'ytick',vt+4,'yticklabel',mindata+vt*(maxdata-mindata));
		if mask,
			set(get(cb,'ylabel'),'string',datatypep);
		else
			set(get(cb,'ylabel'),'string',datatype);
		end
	end;
  
  end
  
  if showparticles,
		vt = [0,0.25,0.5,0.75,1.0];
		cb = colorbar('position',[0.85,0.55,0.02,0.4]);
		set(cb,'ylim',[4,5],'ytick',vt+4,'yticklabel',mindatap+vt*(maxdatap-mindatap));
		set(get(cb,'ylabel'),'string',datatypep);
	end;
 
  
  if anim==0 & useTTB == 0     
	l = light('position',[100 0, 500].*1e3);
  end
%   lightangle(l,zlight(1),zlight(2)); 
%   lightangle(l,inview(1)+5,inview(2)+5); 
    
  hold off;
  
  if anim
    clear gca
	clear cb
	drawnow
  end
  
  if jpg,
    pause(.2);
    
    eval(['print -djpeg90 -r100 ./flic/flic',sprintf('%04d',fnr),'.jpg']);
    
  end;

  if SPM.mesh.docelldata,
    load input/nodelist.mat;
    hold on;
    for i=1:length(nnx),
      plot3(Xc(nny(i),nnx(i)),Yc(nny(i),nnx(i)),3000* ...
            ones(size(nnx)),'.r','markersize',30);
      end;
  end;

end;

