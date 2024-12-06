% Create animation 

SPM = SPMload;
nx = SPM.mesh.nx;
ny = SPM.mesh.ny;
nfiles = 149;
M(nfiles) = struct('cdata',[],'colormap',[]);
M_LR(nfiles) = struct('cdata',[],'colormap',[]);
M_slope(nfiles) = struct('cdata',[],'colormap',[]);
M_bed(nfiles) = struct('cdata',[],'colormap',[]);
M_cirques(nfiles) = struct('cdata',[],'colormap',[]);

%%

%Ice bed evolution
f = figure;
l = light('position',[100,0,500]*1e3);
lightangle(l,60, 90);
for i=1:nfiles
    title(['Time = ',num2str(SPM.mesh.filetime*i/1e6),' Myr'])
    show('file',i,'data','bed','mask',0.01,'ice','tlim',[-500 1500],'dlim',[0 500],'contours','ncon',8,'anim','noclose','view',[90 90],'noforeland')
%     set(gca,'Visible','off')
    cdata = print('-RGBImage','-r300');
    F = im2frame(cdata);
    M(i) = F;
    pause(0.1)
    delete(findobj('Type','surface'))
	delete(findobj('Type','line'))
	delete(findobj('Type','contour'))
end
save('IceBed_movie.mat',"M")


% water
M_water(nfiles) = struct('cdata',[],'colormap',[]);
f = figure;
l = light('position',[100,0,500]*1e3);
lightangle(l,60, 90);
for i=1:nfiles
    title(['Time = ',num2str(SPM.mesh.filetime*i/1e6),' Myr'])
    show('file',i,'data','qw','withice','dlim',[5 9],'log','contours','ncon',8,'anim','noclose','view',[90 90],'noforeland')
%     set(gca,'Visible','off')
    cdata = print('-RGBImage','-r300');
    F = im2frame(cdata);
    M_water(i) = F;
    pause(0.1)
    delete(findobj('Type','surface'))
	delete(findobj('Type','line'))
	delete(findobj('Type','contour'))
end
save('water_movie.mat',"M")

% bed evolution
f = figure;
l = light('position',[100,0,500]*1e3);
lightangle(l,60, 90);
for i=1:nfiles
    title(['Time = ',num2str(SPM.mesh.filetime*i/1e6),' Myr'])
    show('file',i,'data','bed','tlim',[-1000 1500],'contours','ncon',8,'anim','noclose','view',[63 40],'noforeland')
    set(gca,'Visible','off')
    cdata = print('-RGBImage','-r300');
    F = im2frame(cdata);
    M_bed(i) = F;
    pause(0.1)
    delete(findobj('Type','surface'))
	delete(findobj('Type','line'))
	delete(findobj('Type','contour'))
end
save('Bed_movie.mat',"M_bed")

%Low-relief surface evolution
f = figure;
l = light('position',[100,0,500]*1e3);
lightangle(l,60, 90);
for i=1:nfiles
    show('file',i,'data','bed','tlim',[-1000 1500],'contours','ncon',8,'mask2',0.177,'anim','noclose','view',[90 90],'noforeland','radius',10000)
%     set(gca,'Visible','off')
    title(['Time = ',num2str(SPM.mesh.filetime*i/1e6),' Myr'])
    cdata = print('-RGBImage','-r300');
    F = im2frame(cdata);
    M_LR(i) = F;
    pause(0.2)
    delete(findobj('Type','surface'))
	delete(findobj('Type','line'))
	delete(findobj('Type','contour'))
end
save('LR_movie_3Dview.mat',"M_LR")

%%
%cirques
f = figure;
l = light('position',[100,0,500]*1e3);
lightangle(l,60, 90);
for i=1:nfiles
    show('file',i,'data','bed','tlim',[-1000 1500],'contours','ncon',8,'mask',0.1,'cirques','withLRS','cmap','vik','anim','noclose','view',[90 90],'noforeland','radius',5000)
    set(gca,'XTickLabel',[],'YTickLabel',[])
    title(['Time = ',num2str(SPM.mesh.filetime*i/1e6),' Myr'])
    cdata = print('-RGBImage','-r300');
    F = im2frame(cdata);
    M_cirques(i) = F;
    pause(0.2)
    delete(findobj('Type','surface'))
	delete(findobj('Type','line'))
	delete(findobj('Type','contour'))
end
save('cirques_movie_3Dview.mat',"M_cirque")

%%
M_erate(nfiles) = struct('cdata',[],'colormap',[]);
%erate evolution
f = figure;
l=light('position',[100,0,500]*1e3);
lightangle(l,60, 90);
% Initial
show('file',-1,'data','erosion_rate','tlim',[-500 1500],'contours','ncon',10,'dlim',[0 0.5],'anim','noclose','view',[-90 90],'noforeland')
cdata = print('-RGBImage','-r150');
F = im2frame(cdata);
M_erate(1) = F;
pause(0.1)
delete(findobj('Type','surface'))
delete(findobj('Type','Line'))
delete(findobj('Type','contour'))
for i=1:nfiles
    title(['Time = ',num2str(SPM.mesh.filetime*i/1e6),' Myr'])
    show('file',i,'data','erosion_rate','rate',i-1,'dlim',[0 0.5],'tlim',[-500 1500],'contours','ncon',8,'anim','cmap','roma','noclose','view',[-90 90],'noforeland')
    cdata = print('-RGBImage','-r150');
    F = im2frame(cdata);
    M_erate(i+1) = F;
    pause(0.1)
    delete(findobj('Type','surface'))
	delete(findobj('Type','Line'))
	delete(findobj('Type','contour'))
end
save('Erate_movie.mat',"M_erate")


M_abrasion(nfiles+1) = struct('cdata',[],'colormap',[]);
% abrasion evolution
f = figure;
l=light('position',[100,0,500]*1e3);
lightangle(l,60, 90);
% Initial
show('file',-1,'data','abrasion_rate','tlim',[-500 1500],'contours','ncon',10,'dlim',[0 0.5],'anim','noclose','cmap','roma','view',[-90 90],'noforeland')
cdata = print('-RGBImage','-r150');
F = im2frame(cdata);
M_abrasion(1) = F;
pause(0.1)
delete(findobj('Type','surface'))
delete(findobj('Type','Line'))
delete(findobj('Type','contour'))
for i=1:nfiles
    title(['Time = ',num2str(SPM.mesh.filetime*i/1e6),' Myr'])
    show('file',i,'data','abrasion_rate','rate',i-1,'dlim',[0 0.5],'tlim',[-500 1500],'contours','ncon',8,'anim','cmap','roma','noclose','view',[-90 90],'noforeland')
    cdata = print('-RGBImage','-r150');
    F = im2frame(cdata);
    M_abrasion(i+1) = F;
    pause(0.1)
    delete(findobj('Type','surface'))
	delete(findobj('Type','Line'))
	delete(findobj('Type','contour'))
end
save('abrasion_movie.mat',"M_abrasion")

M_fluvial(nfiles+1) = struct('cdata',[],'colormap',[]);
% fluvial evolution
f = figure;
l=light('position',[100,0,500]*1e3);
lightangle(l,60, 90);
% Initial
show('file',-1,'data','fluvial_rate','tlim',[-500 1500],'contours','ncon',10,'dlim',[0 0.5],'anim','noclose','cmap','bamako','view',[-90 90],'noforeland')
cdata = print('-RGBImage','-r150');
F = im2frame(cdata);
M_fluvial(1) = F;
pause(0.1)
delete(findobj('Type','surface'))
delete(findobj('Type','Line'))
delete(findobj('Type','contour'))
for i=1:nfiles
    title(['Time = ',num2str(SPM.mesh.filetime*i/1e6),' Myr'])
    show('file',i,'data','fluvial_rate','rate',i-1,'dlim',[0 0.5],'tlim',[-500 1500],'contours','ncon',8,'anim','cmap','bamako','noclose','view',[-90 90],'noforeland')
    cdata = print('-RGBImage','-r150');
    F = im2frame(cdata);
    M_fluvial(i+1) = F;
    pause(0.1)
    delete(findobj('Type','surface'))
	delete(findobj('Type','Line'))
	delete(findobj('Type','contour'))
end
save('fluvial_movie.mat',"M_fluvial")

M_sliding(nfiles) = struct('cdata',[],'colormap',[]);
%sliding evolution
f = figure;
l=light('position',[100,0,500]*1e3);
lightangle(l,60, 90);
for i=1:nfiles
    title(['Time = ',num2str(SPM.mesh.filetime*i/1e6),' Myr'])
    if i == 1
        show('file',-1,'data','bed','tlim',[-1000 1500],'mask',0,'sliding','contours','ncon',10,'dlim',[0 50],'anim','noclose','view',[-90 90],'noforeland','cmap','cool')
    else
        show('file',i,'data','bed','dlim',[0 50],'tlim',[-1000 1500],'mask',0,'sliding','contours','ncon',8,'anim','noclose','view',[-90 90],'noforeland','cmap','cool')
    end
    cdata = print('-RGBImage','-r150');
    F = im2frame(cdata);
    M_sliding(i) = F;
    pause(0.1)
    delete(findobj('Type','surface'))
	delete(findobj('Type','Line'))
	delete(findobj('Type','contour'))
end
save('sliding_movie.mat',"M_sliding")

M_deformation(nfiles) = struct('cdata',[],'colormap',[]);
%sliding evolution
f = figure;
l=light('position',[100,0,500]*1e3);
lightangle(l,60, 90);
for i=1:nfiles
    title(['Time = ',num2str(SPM.mesh.filetime*i/1e6),' Myr'])
    if i == 1
        show('file',-1,'data','bed','tlim',[-1000 1500],'mask',0,'deformation','contours','ncon',10,'dlim',[0 200],'anim','noclose','view',[-90 90],'noforeland','cmap','cool')
    else
        show('file',i,'data','bed','dlim',[0 200],'tlim',[-1000 1500],'mask',0,'deformation','contours','ncon',8,'anim','noclose','view',[-90 90],'noforeland','cmap','cool')
    end
    cdata = print('-RGBImage','-r150');
    F = im2frame(cdata);
    M_deformation(i) = F;
    pause(0.1)
    delete(findobj('Type','surface'))
	delete(findobj('Type','Line'))
	delete(findobj('Type','contour'))
end
save('IceDeformation_movie.mat',"M_deformation")

%bed slope evolution
f = figure;
l=light('position',[100,0,500]*1e3);
lightangle(l,60, 90);
for i=1:nfiles
    title(['Time = ',num2str(SPM.mesh.filetime*i/1e6),' Myr'])
    show('file',i,'data','bslope','dlim',[0 0.177],'contours','ncon',8,'anim','cmap','lajolla','noclose','view',[90 90],'noforeland')
    %set(gca,'Visible','off')
    cdata = print('-RGBImage','-r150');
    F = im2frame(cdata);
    M_slope(i) = F;
    pause(0.1)
    delete(findobj('Type','surface'))
	delete(findobj('Type','Line'))
	delete(findobj('Type','contour'))
end
save('slopemovie.mat',"M_slope")

%Frequency elev LR evolution
M_LRfreq(nfiles) = struct('cdata',[],'colormap',[]);
f = figure;
for i=1:nfiles
    hypsometry('file',i,'slopeelev','dlim',[-1000 2000],'noclose')
    title(['Time = ',num2str(SPM.mesh.filetime*i/1e6),' Myr'])
    cdata = print('-RGBImage','-r150');
    F = im2frame(cdata);
    M_LRfreq(i) = F;
    pause(0.1)
	delete(findobj('Type','Line'))
    delete(findobj('Type','Patch'))
end
% save('LRfreq_movie_norm.mat',"M_LRfreq")
save('LRfreq_movie.mat',"M_LRfreq")


%%
%% Anim on same fig
load IceBed_movie.mat
load LR_movie.mat
load Erate_movie.mat
load LRfreq_movie.mat

%%
% Save cirque evolution
v = VideoWriter('cirques_movie','MPEG-4');
v.FrameRate = 2;
open(v)
writeVideo(v,M_cirques)
close(v)

%%
% Save abrasion evolution
v = VideoWriter('abrasion_movie','MPEG-4');
v.FrameRate = 2;
open(v)
writeVideo(v,M_abrasion)
close(v)

% Save fluvial evolution
v = VideoWriter('fluvial_movie','MPEG-4');
v.FrameRate = 2;
open(v)
writeVideo(v,M_fluvial)
close(v)

% Save LR frequency evolution
v = VideoWriter('LRfreq_movie','MPEG-4');
v.FrameRate = 2;
open(v)
writeVideo(v,M_LRfreq)
close(v)

% Save erate evolution
v = VideoWriter('erate_movie','MPEG-4');
v.FrameRate = 2;
open(v)
writeVideo(v,M_erate)
close(v)

% Save bed evolution
v = VideoWriter('bed_movie','MPEG-4');
v.FrameRate = 2;
open(v)
writeVideo(v,M_bed)
close(v)

% Save sliding evolution
v = VideoWriter('sliding_movie','MPEG-4');
v.FrameRate = 2;
open(v)
writeVideo(v,M_sliding)
close(v)

% Save slope evolution
v = VideoWriter('slope_movie','MPEG-4');
v.FrameRate = 2;
open(v)
writeVideo(v,M_slope)
close(v)

% Save bed evolution
v = VideoWriter('bed_IceThickness_movie','MPEG-4');
v.FrameRate = 2;
open(v)
writeVideo(v,M)
close(v)

% Save LR evolution
v = VideoWriter('LR_movie','MPEG-4');
v.FrameRate = 2;
open(v)
writeVideo(v,M_LR)
close(v)

%plot mean uplift rate over time
figure
hold on;
for i=1:nfiles
    subplot(1,2,2)
    RGB= frame2im(M_erate(i));
    I1 = imshow(RGB);
    set(I1,'Tag',['Imageb',num2str(i)]')
    zoom(1)
    drawnow
    

    subplot(1,2,1)
    RGB2 = frame2im(M_LR(i));
    I2 = imshow(RGB2);
    set(I2,'Tag',['Imageb',num2str(i)]')
    zoom(1)
    drawnow

    pause(0.2)
    if i > 1
        delete(findobj('Tag',['Image',num2str(i-1)]))
        delete(findobj('Tag',['Imageb',num2str(i-1)]))
    end
    
end
hold off

% plot climate over time
tdata = load('./output/tseries.dat');
figure
plot(tdata(:,1),tdata(:,10),'r')
