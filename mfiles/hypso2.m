function hypso2(fnr)

%close all;
figure
% set(gcf,'units','centimeters','paperunits','centimeters');
% set(gcf,'position',[20,20,20,12],'paperposition',[5,10,10,6]);

%colormap graded to white
colormap('default');
map = colormap;
rd = map(1,1); gd = map(1,2); bd = map(1,3);
xx = linspace(0,1,5);
rx = interp1([0,1],[1,rd],xx);
gx = interp1([0,1],[1,gd],xx);
bx = interp1([0,1],[1,bd],xx);
map2 = [rx(:),gx(:),bx(:);map];
colormap(map2);

SPM=SPMload;
nx=SPM.mesh.nx;
ny=SPM.mesh.ny;
loaddata;
%sedi = sediment(:,4:(end-3));
%I = find(sedi(:) < 0.1); sedi(I) = 0.1;
R = localrelief(bed,4);

ed = linspace(-1000,1000,80);
[N,BIN] = histc(bed(:),ed);

I = find(R(:) < 100);
[Nf,BINf] = histc(bed(I),ed);

Nf = Nf./(N+1);
N = N/length(bed(:));

R = R/600;
edr = linspace(0,1,40);

hold on; box on;

line([-500,2000],[0,0],'color','k');
for i=1:length(ed)-1,
  I = find(BIN == i);
  if length(I) > 10,
    bsa = R(I);
    
    [Nr,BINr] = histc(bsa(:),edr);
    Nr = Nr(:)/length(bsa);
    
    for j=2:(length(edr)-1),
        if Nr(j) > 0.01,
            fill([ed(i),ed(i+1),ed(i+1),ed(i)],[edr(j),edr(j),edr(j+1),edr(j+1)],[Nr(j),Nr(j),Nr(j),Nr(j)],'linestyle','none');
        
        end;
    end;
    
    dx = 0.2*(ed(i+1)-ed(i));
    %xx = [ed(i),ed(i+1)];
    %yy = edr(:);
    %imagesc(xx,yy,[Nr(:),Nr(:)]);
    
    
    patch([ed(i)+dx,ed(i)+dx,ed(i+1)-dx,ed(i+1)-dx],[0,Nf(i),Nf(i),0]-1,[.2,.5,.3],'linewidth',.5);
  
  end;
  
end;
%line(ed(2:end),.4*avslope);
axis([-1000,1000,-1,1]);
set(gca,'xtick',[-500:500:2000],'ytick',[-1:0.25:1]);
%set(gca,'xticklabel',[],'yticklabel',[]);

%eval(['print -depsc -loose hypso',num2str(fnr),'.eps']);

