function hypso(fnr)

%close all;
figure
set(gcf,'units','centimeters','paperunits','centimeters');
%set(gcf,'position',[20,20,20,12],'paperposition',[5,10,10,6]);

SPM=SPMload;
nx=SPM.mesh.nx;
ny=SPM.mesh.ny;
loaddata;
%sedi = sediment(:,4:(end-3));
%I = find(sedi(:) < 0.1); sedi(I) = 0.1;
R = localrelief(bed,4);

ed = linspace(-1000,1000,40);
[N,BIN] = histc(bed(:),ed);

I = find(R(:) < 100);
[Nf,BINf] = histc(bed(I),ed);

Nf = Nf./(N+1)
N = N/length(bed(:));

R = R/600;


hold on; box on;

line([-500,2000],[0,0],'color','k');
for i=1:length(ed)-1,
  I = find(BIN == i);
  if length(I) > 10,
    bsa = R(I);
    Q2 = median(bsa);
    II = find(bsa < Q2); Q1 = median(bsa(II));
    II = find(bsa > Q2); Q3 = median(bsa(II));
    dx = 0.2*(ed(i+1)-ed(i));
    patch([ed(i)+dx,ed(i)+dx,ed(i+1)-dx,ed(i+1)-dx],[Q1,Q3,Q3,Q1], ...
          [.7,.7,.2],'linewidth',.5);
    line([ed(i)+dx,ed(i+1)-dx],[Q2,Q2],'color','k','linewidth',1);
    line([ed(i)+dx,ed(i+1)-dx],[min(bsa),min(bsa)],'color','k','linewidth',.5);
    line([ed(i)+dx,ed(i+1)-dx],[max(bsa),max(bsa)],'color','k','linewidth',.5);
    line(.5*(ed(i)+ed(i+1))*[1,1],[min(bsa),max(bsa)],'color','k','linewidth',.5);
    avslope(i) = mean(bslope(I));
    
    patch([ed(i)+dx,ed(i)+dx,ed(i+1)-dx,ed(i+1)-dx],[0,Nf(i),Nf(i),0]-1,[.2,.5,.3],'linewidth',.5);
  
  end;
  
end;
%line(ed(2:end),.4*avslope);
axis([-1000,1000,-1,1]);
set(gca,'xtick',[-500:500:2000],'ytick',[-1:0.25:1]);
%set(gca,'xticklabel',[],'yticklabel',[]);

%eval(['print -depsc -loose hypso',num2str(fnr),'.eps']);

