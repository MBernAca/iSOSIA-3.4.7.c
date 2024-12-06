% Get uplift rate series

st = load('./status.dat');
latestfnr = st(3);
SPM = SPMload;
nx = SPM.mesh.nx;
ny = SPM.mesh.ny;
filetime = SPM.mesh.filetime; %in yrs
tectonic_upliftRates = SPM.data.srate(50,50);

Uplift_rate = zeros(latestfnr,1);
isostasy_rate = zeros(latestfnr,1);
time = zeros(latestfnr,1);

isostasy_prev = zeros(size(SPM.data.srate));
for i=1:latestfnr
    fnr = i;
    time(i) = filetime*i;
    loaddata;
    isostasy_rate(i) = mean(mean((isostasy - isostasy_prev)))./filetime;
    uplift_rate(i) = isostasy_rate(i) + tectonic_upliftRates;
    isostasy_prev = isostasy;
end

figure
hold on, grid on;
plot(time,isostasy_rate,'b-')
plot(time,uplift_rate,'r-')
hold off;
