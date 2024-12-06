function [SPM,h3,h4]=massbalance(SPM,fg)

%computes and plots massbalance models
%
% SPM is SPM structure
% fg = 0 : no graphics
% fg = 1 : make graphics
%

%close all;

%load SPM data
mPDD = SPM.mprop.mPDD; %Positive degree day melt rate constant
dTa = SPM.mprop.dTa; %annual temperature variation amplitude
precip = 1.0;%acc rates are nomalised to precipitation rate = 1 m/y
             %local acc rates must therefore be scaled by
             %SPM.data.precip inside isosia
Tsl = SPM.mprop.Temp(2,1); %Sea level temperture
Tsl2 = min(SPM.mprop.Temp(2,:));
lrate = SPM.mprop.lrate; %lapse rate

%input 
nt = 365; %days per year
nd = 200; %number of temperature points

%Temperature range
T0 = linspace(-2*dTa,2*dTa,nd);

%Time vector - one year
time = [1:365];

%compute temperature matrix
Temp = T0(:)*ones(1,nt) + ones(nd,1)*dTa*sin(2*pi*time(:)'/365);
% figure
% hold on;
% % for p=1:40:size(Temp,1),
% 	% plot(time,Temp(p,:))
% % end
% plot(time,Temp(1,:))
% plot(time,Temp(size(Temp,1),:))
% hold off;
% ylabel('Temperature (°C)')
% xlabel('Time (days)')

%Compute positive degree day index and number of days with frost
PDD = zeros(nd,1);
nFd = zeros(nd,1);
nFw = zeros(nd,1); %fraction of the year where rain happened - Maxime

for i=1:nd,
    I = find(Temp(i,:) > 0);
    PDD(i) = sum(Temp(i,I));
    nFd(i) = nt - length(I); 
	nFw(i) = length(I)/nt; %Maxime
	mean_temp = mean(Temp(i,:));
	if nFd(i) < 0,
		nFd(i) = 0;
	end
end;

%melt rate for varying T0
melt = mPDD*PDD;
mean(melt)/(max(melt)-min(melt))
%accumulation rate - number of days with frost times average precip rate
macc = precip*nFd/nt;
% macc = precip*nFd;

accgrad = mean(macc)/(max(macc)-min(macc));
%massbalance
mrate = macc - melt;

%Associated elevation range - using first temperature input
h = (Tsl - T0)/lrate;
h5 = (Tsl2-T0)/lrate;
%Evaluate for the minimum sea-level temperature
h2 = (min(SPM.mprop.Temp(2,:)) - T0)/lrate;

%ELA over time
Tsl_array = SPM.mprop.Temp(2,:); %Sea level temperture
h_array = zeros(size(T0));
h5_array = zeros(size(Tsl_array));
for t=1:length(Tsl_array)
	h_array = (Tsl_array(t) - T0)/lrate;
	[~,ind] = min(abs(mrate));
	h5_array(t) = h_array(ind);
end

%Snowline altitude
[~,ind] = min(abs(mrate));
minff = min(abs(mrate));
h3 = h(ind);
h4 = h2(ind);
%save massbalance data to SPM structure
SPM.mprop.Mrate_T = [T0(:)';macc(:)';melt(:)';nFw(:)']; %Maxime

if fg,
	figure
	subplot(1,2,1);
	hold on; box on; grid on;
	li(1)=line(T0,-melt,'color','r');
	li(2)=line(T0,macc,'color','b');
	li(3)=line(T0,mrate,'color','k');
	legend(li,'Melt','Accumulation','Mass balance','Location','SouthWest');
	xlabel('Temperature (C)');
	ylabel('Mass balance (m/y)');

	subplot(1,2,2);
	hold on; box on; grid on;
	li(1)=line(-melt,h,'color','r');
	li(2)=line(macc,h,'color','b');
	li(3)=line(mrate,h,'color','k');
	li(4)=line(mrate,h2,'color','k','linestyle','--');
	li(5) = line(macc,h5,'color','b','linestyle','--');
	legend(li,'Melt','Accumulation','Mass balance','LGM Mb','maccCold','Location','NorthWest');
	ylabel('Elevation (m)');
	xlabel('Melt rate (m/y)');
	
	figure
	li(1)=line(T0,nFw,'color','r');
	li(2)=line(T0,nFd/nt,'color','b');
	legend(li,'rain','snow','Location','NorthWest');
	xlabel('Mean annual temperature (C)');
	ylabel('theta');
	axis square
	
	figure
	li(1)=line(SPM.mprop.Temp(1,:)./1e6,h5_array,'color','r');
	legend(li,'ELA','Location','NorthWest');
	xlabel('Time (Myr)');
	ylabel('Elevation (m)');
	axis square
	index = find(SPM.mprop.Temp(1,:)./1e6 < 1.7);
	disp(['Mean ELA stage 1 = ',num2str(mean(h5_array(index))), ' +/- ', num2str(std(h5_array(index)))])
	disp(['Min ELA stage 1 = ',num2str(min(h5_array(index))), ' m '])
	index = find(SPM.mprop.Temp(1,:)./1e6 > 1.7);
	disp(['Mean ELA stage 2 = ',num2str(mean(h5_array(index))), ' +/- ', num2str(std(h5_array(index)))])
	disp(['Min ELA stage 2 = ',num2str(min(h5_array(index))), ' m '])
	save('massb_param.mat','macc','melt','mrate')
end
