samplePoints= 2048; %%%%CHOOSE VALUE%%%%
Ndetectors = 1;
%Choose to apply electric filter. Put a high value if none
maxDepth = 1.12; %%%%CHOOSE VALUE%%%%
%Remove the photon with very high likelihood
L_filter = 1; %%%%CHOOSE VALUE%%%%
%Compress image for refractive index
n_cor = 1;

%Generate photons
sig = zeros(samplePoints,Ndetectors);
pL = (0:0.001:1)*2*3*1e-3;
W = 1;
L = 1;

lambda_start = 850e-9; %%%%CHOOSE VALUE%%%%
lambda_stop = 950e-9; %%%%CHOOSE VALUE%%%%
k_start = 2*pi/lambda_stop;
k_stop = 2*pi/lambda_start;
k = linspace(k_stop,k_start,samplePoints)';


phase = k*pL;

% pLc = sqrt(W.*L);
% test = pLc/2;
% test2 = test(test<0.13);
% [y,z] = hist(test2,100);
% figure
% plot(db(y))



s = exp(1i.*phase).*sqrt(W.*L);
s_sample = sum(s(:,:),2);
ref_amp = max(mean(abs(s_sample),1));
s_ref = exp(1i.*k*0).*ref_amp;

%Make the interference with the reference arm and calculate the intensity
% I = (abs(s_sample + ref_amp)).^2 - (abs(s_sample - ref_amp) ).^2;
rawAline = (abs(s_sample + s_ref)).^2 - (abs(s_sample - s_ref) ).^2;


OCT = abs(fft(rawAline.*hann(length(k)))).^2;
OCT = OCT(1:floor(length(OCT(:,1))/2)+1,:);
OCT(2:end-1) = 2*OCT(2:end-1);

z = pi/(k(1)-k(2))*(0:(length(k)/2))/length(k)./n_cor;
plot(z,OCT)


%% Create the data of the photons

nAvg = 1000;

W = 1;
L = 1;

samplePoints= 2048; %%%%CHOOSE VALUE%%%%
Ndetectors = 1;
%Choose to apply electric filter. Put a high value if none
maxDepth = 1.12; %%%%CHOOSE VALUE%%%%
%Remove the photon with very high likelihood
L_filter = 1; %%%%CHOOSE VALUE%%%%
%Compress image for refractive index
n_cor = 1;

lambda_start = 850e-9; %%%%CHOOSE VALUE%%%%
lambda_stop = 950e-9; %%%%CHOOSE VALUE%%%%
k_start = 2*pi/lambda_stop;
k_stop = 2*pi/lambda_start;
k = linspace(k_stop,k_start,samplePoints)';


for n = 1:nAvg
    r = exprnd(1e-3*0.5*ones(1,10000));
    W = 1/3*log(rand(1,10000)/(1/(exp(3)-1))+1);
    L = 10*rand(1,10000);
    pL = r;
    phase = k*pL;
    
    s = exp(1i.*phase).*sqrt(W.*L);
    s_sample = sum(s(:,:),2);
    ref_amp = max(mean(abs(s_sample),1));
    s_ref = exp(1i.*k*0).*ref_amp;

    rawAline = (abs(s_sample + s_ref)).^2 - (abs(s_sample - s_ref) ).^2;
    
    OCT = abs(fft(rawAline.*hann(length(k)))).^2;
    OCT = OCT(1:floor(length(OCT(:,1))/2)+1,:);
    OCT(2:end-1) = 2*OCT(2:end-1);

    z = pi/(k(1)-k(2))*(0:(length(k)/2))/length(k)./n_cor;
    OCTt(:,n) = OCT;
    n
    
    
end

dx = z(2)-z(1);
OCTt = OCTt./sum(OCTt*dx*2,1);
OCT = mean(OCTt,2);

%r = exprnd(1e-3*0.5*ones(1,10000));

% add weigth and likelihood
[counts edges] = histcounts(r,100);
for n = 1:100
    sel = and((r > edges(n)),(r < edges(n+1)));
    counts(n) = sum(W(sel).*L(sel));
end

dx = edges(2)-edges(1);
counts = counts./sum(counts*dx);

plot(z*2,OCT,edges(1:end-1)+dx/2,counts)

%%
r = exprnd(1e-3*0.5*ones(1,10000));
[counts cbins] = hist(r,100);
dx = cbins(2)-cbins(1);
counts = counts./sum(counts*dx);
plot(cbins,counts)
pL = r;

dx = z(2)-z(1);
OCT = OCT./sum(OCT*dx);
plot(z,OCT,cbins,counts)

